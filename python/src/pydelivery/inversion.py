"""Trace-by-trace inversion driver and result containers."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import ModelConfig, apply_prior_trace, load_config
from .fluids import FLUID_CHARS, FluidState, enumerate_fluid_states
from .model import StateModel
from .modes import Mode, find_modes, weigh_modes
from .sampler import run_chains
from .seismic import SeismicModel, build_seismic_model, load_wavelet
from .su import Trace, read_su, write_su

# Block order of the Java "realisations.su" format (format version 22).
SU_FIELDS = (
    "d", "LFIV", "NG", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s", "S_h",
    "rho_h", "fluid_type", "vp_eff", "rho_eff", "Z_eff", "R_near", "R_far",
    "seismic_chisquare", "t_fc", "d_fc", "R_near_fc", "R_far_fc", "S_oil",
    "S_lsg", "S_gas", "vs_eff", "Zs_eff", "t",
)
SU_FORMAT_VERSION = 22


@dataclass
class TraceResult:
    """Posterior samples for one trace.

    ``samples[name]`` has shape (n_samples, n_layers) for layer properties,
    and (n_samples,) for 't_base', 'seismic_chisquare' and 'state'.
    Times are in seconds.
    """

    samples: dict[str, np.ndarray]
    states: list[FluidState]
    modes: list[Mode]
    rw_acceptance: float = 0.0
    jump_acceptance: float = 0.0
    header: dict = field(default_factory=dict)

    @property
    def n_samples(self) -> int:
        return len(self.samples["t_base"])

    def fluid_state_probabilities(self) -> dict[str, float]:
        """Posterior probability of each fluid combination (e.g. 'bgb')."""
        counts = np.bincount(self.samples["state"], minlength=len(self.states))
        return {s.label: c / self.n_samples for s, c in zip(self.states, counts)}

    def fluid_probabilities(self) -> np.ndarray:
        """(n_layers, 4) posterior probability of brine/oil/lsg/gas per layer."""
        ft = self.samples["fluid_type"].astype(int)
        return np.stack([(ft == k).mean(axis=0) for k in range(4)], axis=1)

    def pinchout_probability(self) -> np.ndarray:
        """Posterior probability that each layer has zero thickness."""
        t = np.column_stack([self.samples["t"], self.samples["t_base"]])
        return (np.diff(t, axis=1) <= 0).mean(axis=0)

    def thickness(self) -> np.ndarray:
        """Depth thickness of each layer per sample, (n_samples, n_layers)."""
        t = np.column_stack([self.samples["t"], self.samples["t_base"]])
        return 0.5 * np.diff(t, axis=1) * self.samples["vp_eff"]

    def quantile(self, name: str, q) -> np.ndarray:
        x = self.thickness() if name == "thickness" else self.samples[name]
        return np.quantile(x, q, axis=0)

    def mean(self, name: str) -> np.ndarray:
        x = self.thickness() if name == "thickness" else self.samples[name]
        return x.mean(axis=0)

    def std(self, name: str) -> np.ndarray:
        x = self.thickness() if name == "thickness" else self.samples[name]
        return x.std(axis=0)

    def to_su_traces(self) -> list[Trace]:
        """Realisations in the layout written by the Java ``delivery``
        (times in ms, one trace per realisation, numbered in 'mark')."""
        s = self.samples
        n, nl = s["t"].shape
        zeros = np.zeros((n, nl))
        blocks = []
        for name in SU_FIELDS:
            if name == "t":
                blocks.append(1e3 * s["t"])
            elif name == "seismic_chisquare":
                blocks.append(np.repeat(s["seismic_chisquare"][:, None], nl, axis=1))
            elif name in s:
                blocks.append(s[name])
            else:
                blocks.append(zeros)
        data = np.concatenate(blocks + [1e3 * s["t_base"][:, None]], axis=1)
        traces = []
        for i in range(n):
            h = dict(self.header)
            h.update(mark=(i + 1) % 32768, duse=SU_FORMAT_VERSION)
            traces.append(Trace(h, data[i]))
        return traces


def _collect(modes: list[Mode], mode_idx: np.ndarray, ms: list[np.ndarray],
             cfg: ModelConfig, seismic: SeismicModel | None) -> dict[str, np.ndarray]:
    n, nl = len(mode_idx), cfg.n_layers
    names = ("NG", "LFIV", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s", "S_h",
             "rho_h", "vp_eff", "vs_eff", "rho_eff", "R_near", "R_far", "t", "d",
             "fluid_type", "S_oil", "S_lsg", "S_gas")
    out = {k: np.zeros((n, nl)) for k in names}
    out.update(t_base=np.zeros(n), seismic_chisquare=np.zeros(n),
               state=np.zeros(n, dtype=int), mode=mode_idx.copy())
    for j, mo in enumerate(modes):
        idx = np.flatnonzero(mode_idx == j)
        if not len(idx):
            continue
        model = mo.model
        x = np.stack([ms[i] for i in idx])
        p = model.physical(x)
        perm = np.array([lay.is_permeable for lay in cfg.layers])
        for k in ("NG", "LFIV", "vp_m", "vs_m", "rho_m", "vp_eff", "vs_eff", "rho_eff"):
            out[k][idx] = p[k]
        for k in ("vp_s", "vs_s", "phi_s"):
            out[k][idx] = np.where(perm, p[k], 0.0)
        types = np.array(model.state.types)
        out["fluid_type"][idx] = types
        out["S_h"][idx] = np.where(model.has_hc, p["S_h"], 0.0)
        out["rho_h"][idx] = np.where(model.has_hc, p["rho_h"], np.where(perm, p["rho_b"], 0.0))
        for kind, code in (("S_oil", 1), ("S_lsg", 2), ("S_gas", 3)):
            out[kind][idx] = np.where(types == code, p["S_h"], 0.0)
        out["t"][idx] = p["t"][:, :nl]
        out["t_base"][idx] = p["t"][:, nl]
        # depths hung from the top layer's depth using effective velocities
        dz = 0.5 * np.diff(p["t"], axis=1) * p["vp_eff"]
        d0 = p["d"][:, :1]
        out["d"][idx] = d0 + np.concatenate([np.zeros((len(idx), 1)),
                                             np.cumsum(dz[:, :-1], axis=1)], axis=1)
        rts = model.reflection_tables(p)
        if rts:
            out["R_near"][idx] = rts[0]
            if len(rts) > 1:
                out["R_far"][idx] = rts[1]
            ev = model.evaluate(x)
            out["seismic_chisquare"][idx] = 2.0 * ev.seismic_nll / seismic.effective_n
        out["state"][idx] = mo.state_index
    out["Z_eff"] = out["vp_eff"] * out["rho_eff"]
    out["Zs_eff"] = out["vs_eff"] * out["rho_eff"]
    return out


def invert_trace(cfg: ModelConfig, seismic: SeismicModel | None, n_samples: int = 100,
                 rng: np.random.Generator | None = None, *, n_chains: int = 64,
                 burn_in: int = 1000, thin: int = 40, n_random_starts: int = 0,
                 ignore_seismic: bool = False, use_isopachs: bool = True,
                 header: dict | None = None) -> TraceResult:
    """Invert one trace.

    ``cfg`` carries the priors for this trace location and ``seismic`` its
    likelihood set-up (``None`` or ``ignore_seismic=True`` samples the prior).
    """
    rng = rng or np.random.default_rng()
    states = enumerate_fluid_states(cfg)
    modes: list[Mode] = []
    for k, st in enumerate(states):
        model = StateModel(cfg, st, seismic, use_isopachs)
        modes.extend(find_modes(model, k, rng, n_random_starts, ignore_seismic))
    if not modes:
        raise RuntimeError("no posterior modes found (rock physics invalid at all starts?)")
    modes = weigh_modes(modes, [s.probability for s in states])
    chain = run_chains(modes, n_samples, rng, n_chains=n_chains, burn_in=burn_in,
                       thin=thin, ignore_seismic=ignore_seismic)
    samples = _collect(modes, chain.mode, chain.m, cfg, seismic)
    return TraceResult(samples, states, modes, chain.rw_acceptance, chain.jump_acceptance,
                       dict(header or {}))


@dataclass
class Inversion:
    """A whole inversion job: model file, seismic stacks, wavelets, prior traces."""

    cfg: ModelConfig
    stacks: list[list[Trace]]
    wavelets: list
    prior_traces: list[Trace] | None = None

    @classmethod
    def from_files(cls, xml_path, prior_traces=None) -> "Inversion":
        cfg = load_config(xml_path)
        stacks = [read_su(cfg.path(s.filename)) for s in cfg.stacks]
        wavelets = [load_wavelet(cfg.path(s.wavelet_filename)) for s in cfg.stacks]
        pri = read_su(prior_traces) if prior_traces is not None else None
        if cfg.varying and pri is None:
            raise ValueError("model has areally varying priors: prior traces are required")
        if pri is not None and len(pri) != len(stacks[0]):
            raise ValueError("prior traces and seismic traces differ in number")
        return cls(cfg, stacks, wavelets, pri)

    @property
    def n_traces(self) -> int:
        return len(self.stacks[0])

    def trace_setup(self, i: int) -> tuple[ModelConfig, SeismicModel]:
        cfg = self.cfg if self.prior_traces is None else apply_prior_trace(
            self.cfg, self.prior_traces[i].data)
        seis = build_seismic_model(cfg, [s[i] for s in self.stacks], self.wavelets)
        return cfg, seis

    def run(self, n_samples: int = 100, traces=None, seed: int | None = None,
            **kwargs) -> list[TraceResult]:
        rng = np.random.default_rng(seed)
        results = []
        for i in (range(self.n_traces) if traces is None else traces):
            cfg, seis = self.trace_setup(i)
            header = dict(self.stacks[0][i].header)
            results.append(invert_trace(cfg, seis, n_samples, rng, header=header, **kwargs))
        return results


def write_realisations(path: str | Path, results: list[TraceResult]) -> None:
    """Write all realisations as one SU file in the Java Delivery layout."""
    write_su(path, (tr for r in results for tr in r.to_su_traces()))
