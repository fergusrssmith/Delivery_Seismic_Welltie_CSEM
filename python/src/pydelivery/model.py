"""The layered earth model for one trace and one fluid state.

A :class:`StateModel` knows, for a fixed choice of fluid in every layer,

* which quantities are free parameters (the model vector ``m``),
* the Gaussian prior on ``m``, written as a list of linear residuals,
* how to map an *unconstrained* ``m`` to a physical model: values are clipped
  to their physical range and layer times are put in order, so that a layer
  whose base would lie above its top is **pinched out** (zero thickness), and
* the residuals of that physical model against the seismic and any isopach
  constraints.

The negative log posterior density is ``0.5 * sum(residuals**2)`` (plus
constants).  Clipping and time-ordering are what let a continuous Gaussian
prior put finite probability on discrete geological outcomes such as a
pinched-out layer or NG = 0.

Everything is vectorised over a leading batch axis so that many candidate
models are evaluated in one call.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ModelConfig
from .fluids import BRINE, FLUID_NAMES, FluidState
from .rockphysics import effective_properties, reflection_coefficient
from .seismic import SeismicModel

# per-layer property tables, in the order they are stored
LAYER_PROPS = (
    "d", "LFIV", "NG", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s",
    "vp_b", "rho_b", "S_h", "vp_h", "rho_h",
)
UNIT_RANGE = ("NG", "phi_s", "S_h")  # clipped to [0, 1]; all others to >= 0

Key = tuple[str, int]  # (parameter name, layer index); layer -1 for globals


@dataclass
class Evaluation:
    """Result of evaluating a batch of model vectors."""

    residuals: np.ndarray  # (batch, n_res): prior, then seismic, then isopach
    valid: np.ndarray  # (batch,) False where rock physics has no solution
    n_prior: int
    n_seismic: int

    @property
    def nll(self) -> np.ndarray:
        """Negative log posterior (up to constants); +inf where invalid."""
        v = 0.5 * np.sum(self.residuals**2, axis=1)
        return np.where(self.valid, v, np.inf)

    @property
    def seismic_nll(self) -> np.ndarray:
        s = self.residuals[:, self.n_prior : self.n_prior + self.n_seismic]
        return 0.5 * np.sum(s**2, axis=1)


class StateModel:
    def __init__(self, cfg: ModelConfig, state: FluidState, seismic: SeismicModel | None,
                 use_isopachs: bool = True):
        self.cfg = cfg
        self.state = state
        self.seismic = seismic
        nl = self.nl = cfg.n_layers
        types = np.array(state.types)
        self.has_hc = types != BRINE

        # ---- default (prior-mean or fixed) value of every table entry ------
        dflt = {name: np.zeros(nl) for name in LAYER_PROPS}
        self.K_g = np.ones(nl)
        self.rho_g = np.ones(nl)
        self.K_b_ref = np.ones(nl)
        self.rho_b_ref = np.ones(nl)
        self.keys: list[Key] = []  # free parameters, in model-vector order
        self.slots: list[tuple[int, str, int]] = []  # (index in m, table, layer)
        index: dict[Key, int] = {}
        factors: list[tuple[dict[int, float], float, float]] = []  # coeffs, const, sigma

        def free(key: Key, table: str | None = None, layers=None) -> int:
            if key not in index:
                index[key] = len(self.keys)
                self.keys.append(key)
                if table is not None:
                    for lyr in layers if layers is not None else [key[1]]:
                        self.slots.append((index[key], table, lyr))
            return index[key]

        def factor(terms: list[tuple[Key, float]], const: float, sigma: float,
                   fixed: dict[Key, float]) -> None:
            """Add residual (sum_i a_i*x_i - const)/sigma; fixed x go to const."""
            coeffs: dict[int, float] = {}
            for key, a in terms:
                if key in index:
                    coeffs[index[key]] = coeffs.get(index[key], 0.0) + a
                else:
                    const -= a * fixed[key]
            factors.append((coeffs, const, sigma))

        # times first (they dominate the problem)
        self.time_sigma = np.array([lay.sigma_time for lay in cfg.layers] + [cfg.sigma_t_base])
        self.time_default = np.array([lay.time for lay in cfg.layers] + [cfg.t_base])
        self.time_index = np.full(nl + 1, -1)
        fixed: dict[Key, float] = {}
        for i in range(nl + 1):
            if self.time_sigma[i] > 0:
                self.time_index[i] = free(("t", i))
                factor([(("t", i), 1.0)], self.time_default[i], self.time_sigma[i], fixed)
        # times are fixed in order of increasing uncertainty, top-down on ties
        self.time_order = np.argsort(self.time_sigma, kind="stable")

        # depth and LFIV (possibly shared between layers)
        for name, attr in (("d", "depth"), ("LFIV", "lfiv")):
            for i, lay in enumerate(cfg.layers):
                master = getattr(lay, f"{attr}_master")
                mean, sigma = getattr(lay, attr), getattr(lay, f"sigma_{attr}")
                dflt[name][i] = mean
                fixed[(name, i)] = mean
                if sigma > 0 and master == i:
                    users = [k for k, l2 in enumerate(cfg.layers)
                             if getattr(l2, f"{attr}_master") == i]
                    free((name, i), name, users)
                    factor([((name, i), 1.0)], mean, sigma, fixed)

        def rock(i: int, trend, sfx: str, third: str) -> None:
            """Prior for one end-member rock in layer i (sfx 'm' or 's')."""
            lay = cfg.layers[i]
            dk, lk = ("d", lay.depth_master), ("LFIV", lay.lfiv_master)
            vp_bar = (trend.vp_intercept + trend.vp_depth_coeff * lay.depth
                      + trend.vp_lfiv_coeff * lay.lfiv)
            vs_bar = trend.vs_intercept + trend.vs_slope * vp_bar
            if trend.is_reservoir:
                a3, b3, s3 = trend.phi_intercept, trend.phi_slope, trend.phi_sigma
            else:
                # Gardner power law rho = f * vp**e, linearised about vp_bar
                f, e = trend.rho_factor, trend.rho_exponent
                a3 = f * vp_bar**e * (1.0 - e)
                b3 = e * f * vp_bar ** (e - 1.0)
                s3 = trend.rho_sigma
            third_bar = a3 + b3 * vp_bar
            vpk, vsk, thk = (f"vp_{sfx}", i), (f"vs_{sfx}", i), (third, i)
            dflt[vpk[0]][i], dflt[vsk[0]][i], dflt[third][i] = vp_bar, vs_bar, third_bar
            fixed.update({vpk: vp_bar, vsk: vs_bar, thk: third_bar})
            if trend.vp_sigma > 0:
                free(vpk, vpk[0])
                factor([(vpk, 1.0), (dk, -trend.vp_depth_coeff), (lk, -trend.vp_lfiv_coeff)],
                       trend.vp_intercept, trend.vp_sigma, fixed)
            elif dk in index or lk in index:
                raise ValueError(f"rock {trend.name}: sigma_vp=0 needs fixed depth and LFIV")
            if trend.vs_sigma > 0:
                free(vsk, vsk[0])
                factor([(vsk, 1.0), (vpk, -trend.vs_slope)], trend.vs_intercept,
                       trend.vs_sigma, fixed)
            elif vpk in index:
                raise ValueError(f"rock {trend.name}: sigma_vs=0 needs sigma_vp=0")
            if s3 > 0:
                free(thk, third)
                factor([(thk, 1.0), (vpk, -b3)], a3, s3, fixed)
            elif vpk in index:
                raise ValueError(f"rock {trend.name}: zero density/porosity sigma needs sigma_vp=0")

        def scalar(key: Key, table: str, mean: float, sigma: float) -> None:
            dflt[table][key[1]] = mean
            if sigma > 0:
                free(key, table)
                factor([(key, 1.0)], mean, sigma, fixed)

        self.isopachs: list[tuple[int, float, float]] = []
        for i, lay in enumerate(cfg.layers):
            rock(i, cfg.rocks[lay.shale], "m", "rho_m")
            if lay.is_permeable:
                sand = cfg.rocks[lay.sand]
                ref = cfg.fluids[sand.reference_fluid]
                self.K_g[i], self.rho_g[i] = sand.K_grain, sand.rho_grain
                self.rho_b_ref[i], self.K_b_ref[i] = ref.rho, ref.rho * ref.vp**2
                scalar(("NG", i), "NG", lay.ng, lay.sigma_ng)
                rock(i, sand, "s", "phi_s")
                brine = cfg.fluids[lay.brine]
                scalar(("rho_b", i), "rho_b", brine.rho, brine.sigma_rho)
                scalar(("vp_b", i), "vp_b", brine.vp, brine.sigma_vp)
                if self.has_hc[i]:
                    kind = FLUID_NAMES[types[i]]
                    hc = lay.hydrocarbons[kind]
                    fl = cfg.fluids[hc.name]
                    # keys carry the fluid kind: oil and gas parameters in the
                    # same layer are different parameters
                    scalar((f"S_{kind}", i), "S_h", hc.saturation, hc.sigma_saturation)
                    scalar((f"vp_{kind}", i), "vp_h", fl.vp, fl.sigma_vp)
                    scalar((f"rho_{kind}", i), "rho_h", fl.rho, fl.sigma_rho)
            if use_isopachs and lay.thickness is not None and lay.sigma_thickness > 0:
                self.isopachs.append((i, lay.thickness, lay.sigma_thickness))

        self.A_default, self.B_default = cfg.A, cfg.B
        self.A_index = self.B_index = -1
        if cfg.sigma_A > 0:
            self.A_index = free(("A", -1))
            factor([(("A", -1), 1.0)], cfg.A, cfg.sigma_A, fixed)
        if cfg.sigma_B > 0:
            self.B_index = free(("B", -1))
            factor([(("B", -1), 1.0)], cfg.B, cfg.sigma_B, fixed)

        self.defaults = dflt
        self.dim = len(self.keys)
        self.index = index

        # ---- Gaussian prior: residual = Jp @ m - cp ------------------------
        self.Jp = np.zeros((len(factors), self.dim))
        self.cp = np.zeros(len(factors))
        for r, (coeffs, const, sigma) in enumerate(factors):
            for j, a in coeffs.items():
                self.Jp[r, j] = a / sigma
            self.cp[r] = const / sigma
        self.prior_precision = self.Jp.T @ self.Jp
        self.prior_mean = np.linalg.solve(self.prior_precision, self.Jp.T @ self.cp)
        self.prior_cov = np.linalg.inv(self.prior_precision)
        self.logdet_prior_cov = -np.linalg.slogdet(self.prior_precision)[1]

        self._slot_idx = {}
        for name in LAYER_PROPS:
            s = [(j, lyr) for j, tab, lyr in self.slots if tab == name]
            self._slot_idx[name] = (np.array([a for a, _ in s], dtype=int),
                                    np.array([b for _, b in s], dtype=int))

    # ------------------------------------------------------------------
    def tables(self, m: np.ndarray) -> dict[str, np.ndarray]:
        """Physical (clipped, time-ordered) property tables for model vectors
        ``m`` of shape (batch, dim).  Adds 't' (batch, nl+1), 'A', 'B'."""
        m = np.atleast_2d(m)
        batch = m.shape[0]
        p = {}
        for name in LAYER_PROPS:
            tab = np.broadcast_to(self.defaults[name], (batch, self.nl)).copy()
            j, lyr = self._slot_idx[name]
            if len(j):
                tab[:, lyr] = m[:, j]
            np.maximum(tab, 0.0, out=tab)
            if name in UNIT_RANGE:
                np.minimum(tab, 1.0, out=tab)
            p[name] = tab
        t = np.broadcast_to(self.time_default, (batch, self.nl + 1)).copy()
        free_t = self.time_index >= 0
        t[:, free_t] = m[:, self.time_index[free_t]]
        p["t"] = self.order_times(t)
        p["A"] = m[:, self.A_index] if self.A_index >= 0 else np.full(batch, self.A_default)
        p["B"] = m[:, self.B_index] if self.B_index >= 0 else np.full(batch, self.B_default)
        return p

    def order_times(self, t: np.ndarray) -> np.ndarray:
        """Enforce t[i] <= t[i+1] by truncation.

        Times are fixed one at a time, best-known (smallest prior sigma)
        first; each later one is clipped between the already-fixed times
        above and below it.  A layer whose top is clipped onto its base has
        zero thickness: it is pinched out.
        """
        out = t.copy()
        n = t.shape[1]
        lo = np.full_like(t, -np.inf)  # running max of fixed times above
        hi = np.full_like(t, np.inf)  # running min of fixed times below
        for k in self.time_order:
            v = np.clip(out[:, k], lo[:, k], hi[:, k])
            out[:, k] = v
            if k + 1 < n:
                lo[:, k + 1 :] = np.maximum(lo[:, k + 1 :], v[:, None])
            if k > 0:
                hi[:, :k] = np.minimum(hi[:, :k], v[:, None])
        return out

    def physical(self, m: np.ndarray) -> dict[str, np.ndarray]:
        """Property tables plus effective properties, reflection coefficients
        and bookkeeping for output."""
        p = self.tables(m)
        vp, vs, rho, valid = effective_properties(
            p, self.K_g, self.rho_g, self.K_b_ref, self.rho_b_ref, self.has_hc)
        p.update(vp_eff=vp, vs_eff=vs, rho_eff=rho, valid=valid)
        t = p["t"]
        p["present"] = t[:, 1:] > t[:, :-1]
        return p

    def _reflections(self, p: dict, denom: float) -> np.ndarray:
        """Reflection coefficient at the top of each layer, (batch, nl).

        A reflection exists at the top of every layer of non-zero thickness
        (except the first such layer) and is the contrast with the nearest
        non-zero-thickness layer above, so pinched-out layers are skipped."""
        present = p["present"]
        batch = present.shape[0]
        rc = np.zeros((batch, self.nl))
        have = np.zeros(batch, dtype=bool)
        a_vp = np.zeros(batch)
        a_vs = np.zeros(batch)
        a_rho = np.ones(batch)
        for k in range(self.nl):
            pk = present[:, k]
            vp, vs, rho = p["vp_eff"][:, k], p["vs_eff"][:, k], p["rho_eff"][:, k]
            if k > 0:
                both = pk & have
                if both.any():
                    with np.errstate(divide="ignore", invalid="ignore"):
                        r = reflection_coefficient(a_rho, a_vp, a_vs, rho, vp, vs,
                                                   p["A"], p["B"], denom)
                    rc[:, k] = np.where(both, r, 0.0)
            a_vp = np.where(pk, vp, a_vp)
            a_vs = np.where(pk, vs, a_vs)
            a_rho = np.where(pk, rho, a_rho)
            have |= pk
        return np.nan_to_num(rc)

    def evaluate(self, m: np.ndarray, ignore_seismic: bool = False) -> Evaluation:
        """Residuals of a batch of model vectors, shape (batch, dim)."""
        m = np.atleast_2d(m)
        p = self.physical(m)
        parts = [m @ self.Jp.T - self.cp]
        n_seis = 0
        if self.seismic is not None:
            for sm in self.seismic.stacks:
                rc = self._reflections(p, sm.theta_sq_denominator)
                R = self.seismic.reflectivity(p["t"][:, : self.nl], rc)
                res = (sm.observed - R @ sm.G.T) * sm.scale
                parts.append(np.zeros_like(res) if ignore_seismic else res)
                n_seis += res.shape[1]
        for i, dz, sigma in self.isopachs:
            thick = 0.5 * (p["t"][:, i + 1] - p["t"][:, i]) * p["vp_eff"][:, i]
            parts.append(((thick - dz) / sigma)[:, None])
        return Evaluation(np.concatenate(parts, axis=1), p["valid"], self.Jp.shape[0], n_seis)

    def reflection_tables(self, p: dict) -> list[np.ndarray]:
        if self.seismic is None:
            return []
        return [self._reflections(p, sm.theta_sq_denominator) for sm in self.seismic.stacks]

    def synthetics(self, m: np.ndarray) -> list[np.ndarray]:
        """Synthetic trace(s) on the model time grid, one array per stack."""
        p = self.physical(m)
        out = []
        for k, sm in enumerate(self.seismic.stacks):
            rc = self._reflections(p, sm.theta_sq_denominator)
            R = self.seismic.reflectivity(p["t"][:, : self.nl], rc)
            out.append(self.seismic.synthetic(k, R))
        return out
