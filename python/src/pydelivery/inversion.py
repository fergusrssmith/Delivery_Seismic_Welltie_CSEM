"""Trace-by-trace inversion driver and result containers."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .config import ModelConfig, apply_prior_trace, load_config
from .fluids import FluidState, enumerate_fluid_states
from .io import TraceSet, read_prior_table, read_traces, write_results
from .model import StateModel
from .modes import Mode, find_modes, weigh_modes
from .sampler import run_chains
from .seismic import SeismicModel, build_seismic_model, load_wavelet

# per-layer outputs taken straight from StateModel.physical()
LAYER_OUTPUTS = ("t", "NG", "LFIV", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s",
                 "NG_hard", "vp_hard", "vs_hard", "rho_hard", "float", "vp_b", "rho_b",
                 "vp_eff", "vs_eff", "rho_eff", "S_h", "fluid_type", "z")


@dataclass
class TraceResult:
    """Posterior samples for one trace.

    ``samples[name]`` has shape (n_samples, n_layers) for layer properties
    and (n_samples,) for per-sample scalars ('t_base', 'z_base',
    'seismic_chisquare', 'state', 'A', 'B'). Times are seconds; depths are in
    the XML's units. 'z' is the depth of each layer top, hung from the
    reference layer; 'S_oil', 'S_lsg', 'S_gas' are thickness-weighted
    saturations (they differ from 'S_h' only in contacts mode).
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
        """(n_layers, 4) posterior probability of brine/oil/lsg/gas per layer
        (in contacts mode: of the dominant fluid of each layer)."""
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

    def _get(self, name):
        return self.thickness() if name == "thickness" else self.samples[name]

    def quantile(self, name: str, q) -> np.ndarray:
        return np.quantile(self._get(name), q, axis=0)

    def mean(self, name: str) -> np.ndarray:
        return self._get(name).mean(axis=0)

    def std(self, name: str) -> np.ndarray:
        return self._get(name).std(axis=0)


def _collect(modes: list[Mode], mode_idx: np.ndarray, ms: list[np.ndarray],
             cfg: ModelConfig, seismic: SeismicModel | None) -> dict[str, np.ndarray]:
    n, nl = len(mode_idx), cfg.n_layers
    n_stacks = len(seismic.stacks) if seismic is not None else 0
    out: dict[str, np.ndarray] = {}
    for j, mo in enumerate(modes):
        idx = np.flatnonzero(mode_idx == j)
        if not len(idx):
            continue
        model = mo.model
        x = np.stack([ms[i] for i in idx])
        p = model.physical(x)
        vals = {k: p[k] for k in LAYER_OUTPUTS}
        vals["t"] = p["t"][:, :nl]
        vals["t_base"] = p["t"][:, nl]
        vals["z_base"] = p["z_base"]
        vals["A"], vals["B"] = p["A"], p["B"]
        for h in ("oil", "lsg", "gas"):
            vals[f"S_{h}"] = p[f"S_eff_{h}"]
        if p["contact_depth"].shape[1]:
            vals["contact_depth"] = p["contact_depth"]
        if model.contacts_mode:
            vals["n_sublayers"] = p["n_sublayers"]
        for k in range(n_stacks):
            top, internal = model.layer_reflections(p, k)
            vals[f"R_{k}"] = top
            if model.contacts_mode:
                vals[f"R_contact_{k}"] = internal
                vals["t_contact"] = p["sub_t0"][..., 1:] * (p["sub_t1"][..., 1:] > p["sub_t0"][..., 1:])
        if seismic is not None:
            ev = model.evaluate(x)
            vals["seismic_chisquare"] = 2.0 * ev.seismic_nll / seismic.effective_n
        vals["state"] = np.full(len(idx), mo.state_index)
        vals["mode"] = np.full(len(idx), j)
        for k, v in vals.items():
            if k not in out:
                out[k] = np.zeros((n,) + v.shape[1:], dtype=v.dtype)
            out[k][idx] = v
    out["Z_eff"] = out["vp_eff"] * out["rho_eff"]
    out["Zs_eff"] = out["vs_eff"] * out["rho_eff"]
    return out


def invert_trace(cfg: ModelConfig, seismic: SeismicModel | None, n_samples: int = 100,
                 rng: np.random.Generator | None = None, *, n_chains: int = 64,
                 burn_in: int = 1000, thin: int = 40, n_random_starts: int | None = None,
                 ignore_seismic: bool = False, use_isopachs: bool = True,
                 use_stacking: bool = True, header: dict | None = None) -> TraceResult:
    """Invert one trace.

    ``cfg`` carries the priors for this trace location and ``seismic`` its
    likelihood set-up (``None`` or ``ignore_seismic=True`` samples the prior).
    """
    rng = rng or np.random.default_rng()
    if n_random_starts is None:
        # contact crossings make the posterior lumpier: search more widely
        n_random_starts = 8 if cfg.contacts_mode else 0
    states = enumerate_fluid_states(cfg)
    modes: list[Mode] = []
    for k, st in enumerate(states):
        model = StateModel(cfg, st, seismic, use_isopachs, use_stacking)
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
    """A whole inversion job: model file, seismic stacks, wavelets, prior table."""

    cfg: ModelConfig
    stacks: list[TraceSet]
    wavelets: list
    priors: dict[str, np.ndarray] | None = None

    @classmethod
    def from_files(cls, xml_path, priors=None) -> "Inversion":
        """Load a job. Seismic and wavelet files named in the XML may be
        .h5, .npz, .sgy or legacy .su; ``priors`` is an optional prior table
        (.h5/.npz named arrays, or a Java-format .su/.sgy block file)."""
        cfg = load_config(xml_path)
        stacks = [read_traces(cfg.path(s.filename)) for s in cfg.stacks]
        wavelets = [load_wavelet(cfg.path(s.wavelet_filename)) for s in cfg.stacks]
        table = read_prior_table(priors, cfg) if priors is not None else None
        if any(len(s) != len(stacks[0]) for s in stacks):
            raise ValueError("stacks have different numbers of traces")
        if table is not None:
            for k, v in table.items():
                if len(v) != len(stacks[0]):
                    raise ValueError(f"prior table entry '{k}' has {len(v)} traces, "
                                     f"seismic has {len(stacks[0])}")
        return cls(cfg, stacks, wavelets, table)

    @property
    def n_traces(self) -> int:
        return len(self.stacks[0])

    def trace_setup(self, i: int) -> tuple[ModelConfig, SeismicModel]:
        cfg = self.cfg if self.priors is None else apply_prior_trace(self.cfg, self.priors, i)
        seis = build_seismic_model(cfg, [s[i] for s in self.stacks], self.wavelets)
        return cfg, seis

    def run(self, n_samples: int = 100, traces=None, seed: int | None = None,
            **kwargs) -> list[TraceResult]:
        rng = np.random.default_rng(seed)
        results = []
        for i in (range(self.n_traces) if traces is None else traces):
            cfg, seis = self.trace_setup(i)
            header = self.stacks[0][i].header
            results.append(invert_trace(cfg, seis, n_samples, rng, header=header, **kwargs))
        return results

    def save(self, path, results: list[TraceResult]) -> None:
        """Write results to .h5 or .npz."""
        write_results(path, results, self.cfg)


__all__ = ["Inversion", "TraceResult", "invert_trace", "write_results"]
