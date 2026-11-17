"""The layered earth model for one trace and one fluid state.

A :class:`StateModel` knows, for a fixed choice of fluid in every layer (or,
in contacts mode, a fixed set of fluid contacts),

* which quantities are free parameters (the model vector ``m``), including
  variables shared between layers (``same_as_layer_above``,
  ``synchronise_rock_matrix``, synchronised fluid units, B tied to A);
* the Gaussian prior on ``m``, written as a list of linear residuals
  (trend curves, vertical-stacking difference terms, ...);
* how to map an *unconstrained* ``m`` to a physical model: values are clipped
  to their physical range and layer times are put in order, so that a layer
  whose base would lie above its top is **pinched out** (zero thickness); and
* the residuals of that physical model against the seismic, isopachs and
  vertical-stacking thickness constraints.

The negative log posterior density is ``0.5 * sum(residuals**2)`` plus
constants. Everything is vectorised over a leading batch axis.

Contacts mode: layer fluids are not fixed per state but follow from depth
contacts. Depths are hung from the reference layer using the effective
velocities; a layer crossed by a contact is split into sub-layers, each
fluid-substituted separately, and the contact itself reflects.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import HYDROCARBONS, ModelConfig
from .fluids import BRINE, FLUID_NAMES, FluidState
from .rockphysics import effective_properties, impermeable_properties, reflection_coefficient
from .seismic import SeismicModel

ROCK_PROPS = ("NG", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s",
              "NG_hard", "vp_hard", "vs_hard", "rho_hard")
HC_PROPS = tuple(f"{p}_{h}" for h in HYDROCARBONS for p in ("S", "vp", "rho"))
LAYER_PROPS = ("d", "LFIV", "float") + ROCK_PROPS + ("vp_b", "rho_b") + HC_PROPS
UNIT_RANGE = {"NG", "phi_s", "NG_hard", "S_oil", "S_lsg", "S_gas"}
MAX_SUBLAYERS = 4

Key = tuple[str, int]  # (parameter name, owning layer); -1 for globals


@dataclass
class Evaluation:
    """Result of evaluating a batch of model vectors."""

    residuals: np.ndarray  # (batch, n_res): prior, seismic, isopach, stacking
    valid: np.ndarray  # (batch,) False where the rock physics has no solution
    n_prior: int
    n_seismic: int

    @property
    def nll(self) -> np.ndarray:
        """Negative log posterior (up to constants); +inf where invalid."""
        v = 0.5 * np.sum(self.residuals**2, axis=1)
        return np.where(self.valid, v, np.inf)

    @property
    def seismic_nll(self) -> np.ndarray:
        s = self.residuals[:, self.n_prior:self.n_prior + self.n_seismic]
        return 0.5 * np.sum(s**2, axis=1)


class StateModel:
    def __init__(self, cfg: ModelConfig, state: FluidState, seismic: SeismicModel | None,
                 use_isopachs: bool = True, use_stacking: bool = True):
        self.cfg = cfg
        self.state = state
        self.seismic = seismic
        nl = self.nl = cfg.n_layers
        L = cfg.layers
        types = np.array(state.types)
        self.contacts_mode = cfg.contacts_mode
        self.has_hc = (types != BRINE) & (not self.contacts_mode)
        self.fluid_kind = [FLUID_NAMES[t] for t in types]
        self.perm = np.array([lay.is_permeable for lay in L])
        self.active = np.array([lay.active for lay in L])
        self.has_hard = np.array([lay.has_hard_rock for lay in L])

        self.keys: list[Key] = []
        index: dict[Key, int] = {}
        slots: list[tuple[int, str, int]] = []  # (index in m, table, layer)
        factors: list[tuple[dict[int, float], float, float]] = []
        value: dict[Key, float] = {}  # prior mean / fixed value of every key
        dflt = {name: np.zeros(nl) for name in LAYER_PROPS}
        self.K_g = np.ones(nl)
        self.rho_g = np.ones(nl)
        self.K_b_ref = np.ones(nl)
        self.rho_b_ref = np.ones(nl)

        def bind(table: str, layer: int, key: Key, mean: float, sigma: float,
                 prior: bool = True) -> None:
            """Layer ``layer``'s entry of ``table`` is the variable ``key``."""
            value.setdefault(key, mean)
            dflt[table][layer] = value[key]
            if sigma > 0:
                if key not in index:
                    index[key] = len(self.keys)
                    self.keys.append(key)
                    if prior:
                        factor([(key, 1.0)], mean, sigma)
                slots.append((index[key], table, layer))

        def factor(terms: list[tuple[Key, float]], const: float, sigma: float) -> None:
            """Residual (sum_i a_i*x_i - const)/sigma; fixed x move into const."""
            coeffs: dict[int, float] = {}
            for key, a in terms:
                if key in index:
                    coeffs[index[key]] = coeffs.get(index[key], 0.0) + a
                else:
                    const -= a * value[key]
            factors.append((coeffs, const, sigma))

        # ---- layer times (the dominant variables come first) -----------------
        # a forced pinch-out (inactive layer) ties its base time to its top
        self.time_index = np.full(nl + 1, -1)
        self.time_default = np.zeros(nl + 1)
        sigmas = np.array([lay.sigma_time for lay in L] + [cfg.sigma_t_base])
        means = np.array([lay.time for lay in L] + [cfg.t_base])
        owner = list(range(nl + 1))
        for k in range(1, nl + 1):
            if not L[k - 1].active:
                owner[k] = owner[k - 1]
        for k in range(nl + 1):
            o = owner[k]
            key = ("t", o)
            value.setdefault(key, means[o])
            self.time_default[k] = value[key]
            if sigmas[o] > 0:
                if key not in index:
                    index[key] = len(self.keys)
                    self.keys.append(key)
                    factor([(key, 1.0)], means[o], sigmas[o])
                self.time_index[k] = index[key]
        # times are fixed in order of increasing uncertainty, top-down on ties
        self.time_order = np.argsort(np.array([sigmas[o] for o in owner]), kind="stable")

        # ---- which layers' rock variables are needed ------------------------
        needed = self.active.copy()
        for j, lay in enumerate(L):
            if lay.active:
                needed[lay.rock_master] = True
        depth_users = {L[j].depth_master for j in range(nl) if needed[j]} | {L[cfg.master_depth_layer].depth_master}
        lfiv_users = {L[j].lfiv_master for j in range(nl) if needed[j]}
        float_users = {L[j].float_master for j in range(nl) if needed[j] and L[j].is_permeable}

        # ---- depth, rock-curve depth, LFIV, float fraction -------------------
        for i, lay in enumerate(L):
            m = lay.depth_master
            bind("d", i, ("d", m), L[m].depth, L[m].sigma_depth if m in depth_users else 0.0)
            m = lay.lfiv_master
            bind("LFIV", i, ("LFIV", m), L[m].lfiv, L[m].sigma_lfiv if m in lfiv_users else 0.0)
            value[("dc", lay.rock_curves_master)] = (
                lay.depth_rock_curves if lay.depth_rock_curves is not None else lay.depth)
            m = lay.float_master
            bind("float", i, ("float", m), L[m].float_fraction,
                 L[m].sigma_float_fraction if m in float_users else 0.0)

        def trend_depth_key(lay) -> Key:
            # trends use depth_rock_curves (fixed) when given, else the model depth
            if lay.depth_rock_curves is not None:
                return ("dc", lay.rock_curves_master)
            return ("d", lay.depth_master)

        def rock(i: int, trend, vp_t: str, vs_t: str, third_t: str, free: bool = True) -> None:
            """Trend-curve prior for one end-member rock of layer i. With
            ``free=False`` only the prior means are recorded (inactive layer)."""
            lay = L[i]
            o = lay.rock_master
            dk, lk = trend_depth_key(lay), ("LFIV", lay.lfiv_master)
            fk = ("float", lay.float_master)
            ff = value.get(fk, 0.0) if trend.is_reservoir else 0.0
            vp_bar = (trend.vp_intercept + trend.vp_depth_coeff * value[dk]
                      + trend.vp_lfiv_coeff * value[lk])
            if trend.is_reservoir:
                vp_bar += trend.vp_float_coeff * ff
                a3, b3, c3, s3 = trend.phi_intercept, trend.phi_slope, trend.phi_float_coeff, trend.phi_sigma
            else:
                # Gardner power law rho = f*vp**e, linearised about vp_bar
                f, e = trend.rho_factor, trend.rho_exponent
                a3, b3, c3, s3 = f * vp_bar**e * (1.0 - e), e * f * vp_bar ** (e - 1.0), 0.0, trend.rho_sigma
            vs_bar = trend.vs_intercept + trend.vs_slope * vp_bar
            third_bar = a3 + b3 * vp_bar + c3 * ff
            vpk, vsk, thk = (vp_t, o), (vs_t, o), (third_t, o)
            synched = o != i
            if not free:
                for tab, key, mean in ((vp_t, vpk, vp_bar), (vs_t, vsk, vs_bar), (third_t, thk, third_bar)):
                    bind(tab, i, key, mean, 0.0)
                return
            if trend.vp_sigma == 0 and not synched and (dk in index or lk in index):
                raise ValueError(f"rock {trend.name}: sigma_vp=0 needs fixed depth and LFIV")
            bind(vp_t, i, vpk, vp_bar, trend.vp_sigma, prior=False)
            bind(vs_t, i, vsk, vs_bar, trend.vs_sigma, prior=False)
            bind(third_t, i, thk, third_bar, s3, prior=False)
            if synched:
                return  # the master layer carries the prior
            if trend.vp_sigma > 0:
                terms = [(vpk, 1.0), (dk, -trend.vp_depth_coeff), (lk, -trend.vp_lfiv_coeff)]
                if trend.is_reservoir and trend.vp_float_coeff:
                    terms.append((fk, -trend.vp_float_coeff))
                factor(terms, trend.vp_intercept, trend.vp_sigma)
            if trend.vs_sigma > 0:
                factor([(vsk, 1.0), (vpk, -trend.vs_slope)], trend.vs_intercept, trend.vs_sigma)
            elif vpk in index:
                raise ValueError(f"rock {trend.name}: sigma_vs=0 needs sigma_vp=0")
            if s3 > 0:
                terms = [(thk, 1.0), (vpk, -b3)]
                if c3:
                    terms.append((fk, -c3))
                factor(terms, a3, s3)
            elif vpk in index:
                raise ValueError(f"rock {trend.name}: zero density/porosity sigma needs sigma_vp=0")

        # ---- rocks and fluids, layer by layer -------------------------------
        unit_master: dict[tuple[int, str], int] = {}  # (unit id, hc kind) -> layer
        if self.contacts_mode:
            for u_id, u in enumerate(cfg.fluid_units):
                if not u.synchronised:
                    continue
                kinds = {FLUID_NAMES[c] for ct in u.contacts for c in (ct.upper, ct.lower)} - {"brine"}
                for kind in kinds:
                    for j in range(u.top, u.bottom + 1):
                        if L[j].is_permeable and L[j].active:
                            unit_master[(u_id, kind)] = j
                            break

        self.isopachs: list[tuple[int, float, float]] = []
        for i, lay in enumerate(L):
            o = lay.rock_master
            on = bool(needed[i])
            if not lay.is_pure_sand:
                rock(i, cfg.rocks[lay.shale], "vp_m", "vs_m", "rho_m", free=on)
            if lay.has_hard_rock:
                hr = lay.hard_rock_for(self.fluid_kind[i])
                if hr is None:
                    raise ValueError(f"layer {i + 1}: no hard rock given for {self.fluid_kind[i]}")
                bind("NG_hard", i, ("NG_hard", o), L[o].ng_hard, L[o].sigma_ng_hard if on else 0.0,
                     prior=o == i)
                rock(i, cfg.rocks[hr], "vp_hard", "vs_hard", "rho_hard", free=on)
            if lay.is_permeable:
                sand = cfg.rocks[lay.sand]
                ref = cfg.fluids[sand.reference_fluid]
                self.K_g[i], self.rho_g[i] = sand.K_grain, sand.rho_grain
                self.rho_b_ref[i], self.K_b_ref[i] = ref.rho, ref.rho * ref.vp**2
                bind("NG", i, ("NG", o), L[o].ng, L[o].sigma_ng if on else 0.0, prior=o == i)
                rock(i, sand, "vp_s", "vs_s", "phi_s", free=on)
                brine = cfg.fluids[lay.brine]
                act = lay.active
                bind("rho_b", i, ("rho_b", i), _pick(lay.rho_brine, brine.rho),
                     _pick(lay.sigma_rho_brine, brine.sigma_rho) if act else 0.0)
                bind("vp_b", i, ("vp_b", i), _pick(lay.vp_brine, brine.vp),
                     _pick(lay.sigma_vp_brine, brine.sigma_vp) if act else 0.0)
                for kind in self._layer_hydrocarbons(i):
                    hc = lay.hydrocarbons[kind]
                    fl = cfg.fluids[hc.name]
                    h_owner = i
                    if self.contacts_mode:
                        u_id = cfg.fluid_units.index(cfg.fluid_unit_of(i))
                        h_owner = unit_master.get((u_id, kind), i)
                    src = L[h_owner].hydrocarbons[kind]
                    for tab, mean, sig in (
                            (f"S_{kind}", src.saturation, src.sigma_saturation),
                            (f"vp_{kind}", _pick(src.vp, fl.vp), _pick(src.sigma_vp, fl.sigma_vp)),
                            (f"rho_{kind}", _pick(src.rho, fl.rho), _pick(src.sigma_rho, fl.sigma_rho))):
                        bind(tab, i, (tab, h_owner), mean, sig if act else 0.0, prior=h_owner == i)
            if (use_isopachs and lay.active and lay.thickness is not None
                    and lay.sigma_thickness > 0):
                self.isopachs.append((i, lay.thickness, lay.sigma_thickness))

        # ---- contacts, AVO -------------------------------------------------
        self.contact_index = []
        self.contact_default = np.array([c.depth for c in cfg.contacts])
        for c_id, c in enumerate(cfg.contacts):
            key = ("contact", c_id)
            value[key] = c.depth
            if c.sigma_depth > 0:
                index[key] = len(self.keys)
                self.keys.append(key)
                factor([(key, 1.0)], c.depth, c.sigma_depth)
            self.contact_index.append(index.get(key, -1))
        self.contact_index = np.array(self.contact_index, dtype=int)
        self.A_default, self.B_default = cfg.A, cfg.B
        self.A_index = self.B_index = -1
        if cfg.sigma_A > 0:
            self.A_index = index[("A", -1)] = len(self.keys)
            self.keys.append(("A", -1))
            value[("A", -1)] = cfg.A
            factor([(("A", -1), 1.0)], cfg.A, cfg.sigma_A)
        if cfg.sync_B_from_A:
            self.B_index = self.A_index
        elif cfg.sigma_B > 0:
            self.B_index = index[("B", -1)] = len(self.keys)
            self.keys.append(("B", -1))
            value[("B", -1)] = cfg.B
            factor([(("B", -1), 1.0)], cfg.B, cfg.sigma_B)

        # ---- vertical model stacking ---------------------------------------
        self.thickness_couplings: list[tuple[int, int, float]] = []
        if use_stacking and cfg.stacking is not None:
            owner_of = {p: (lambda j: L[j].rock_master) for p in ROCK_PROPS}
            for b, targets in cfg.stacking.mappings:
                for prop, (sigma, applicable) in cfg.stacking.sigmas.items():
                    if applicable is not None and b not in applicable:
                        continue
                    for tl in targets:
                        if prop == "thickness":
                            self.thickness_couplings.append((b, tl, sigma))
                            continue
                        kb, kt = (prop, owner_of[prop](b)), (prop, owner_of[prop](tl))
                        if kb in index and kt in index and kb != kt:
                            factor([(kb, 1.0), (kt, -1.0)], 0.0, sigma)

        self.defaults = dflt
        self.dim = len(self.keys)
        self.index = index
        self.slots = slots

        # ---- Gaussian prior: residual = Jp @ m - cp --------------------------
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
            s = [(j, lyr) for j, tab, lyr in slots if tab == name]
            self._slot_idx[name] = (np.array([a for a, _ in s], dtype=int),
                                    np.array([b for _, b in s], dtype=int))
        self._sand = {i: cfg.rocks[L[i].sand] for i in range(nl) if L[i].is_permeable}

    # ------------------------------------------------------------------
    def _layer_hydrocarbons(self, i: int) -> list[str]:
        lay = self.cfg.layers[i]
        if not lay.active:
            return []
        if self.contacts_mode:
            u = self.cfg.fluid_unit_of(i)
            codes = {c for ct in u.contacts for c in (ct.upper, ct.lower)} - {BRINE}
            return [FLUID_NAMES[c] for c in sorted(codes)]
        kind = self.fluid_kind[i]
        if kind == "brine":
            return []
        return [kind]

    def tables(self, m: np.ndarray) -> dict[str, np.ndarray]:
        """Physical (clipped, time-ordered) property tables for model vectors
        ``m`` of shape (batch, dim)."""
        m = np.atleast_2d(m)
        batch = m.shape[0]
        p = {}
        for name in LAYER_PROPS:
            tab = np.broadcast_to(self.defaults[name], (batch, self.nl)).copy()
            j, lyr = self._slot_idx[name]
            if len(j):
                tab[:, lyr] = m[:, j]
            p[name] = tab
        self._remap_negative_float(p)
        for name in LAYER_PROPS:
            if name in ("float",):
                continue
            np.maximum(p[name], 0.0, out=p[name])
            if name in UNIT_RANGE:
                np.minimum(p[name], 1.0, out=p[name])
        t = np.broadcast_to(self.time_default, (batch, self.nl + 1)).copy()
        free_t = self.time_index >= 0
        t[:, free_t] = m[:, self.time_index[free_t]]
        p["t"] = self.order_times(t)
        p["A"] = m[:, self.A_index] if self.A_index >= 0 else np.full(batch, self.A_default)
        p["B"] = m[:, self.B_index] if self.B_index >= 0 else np.full(batch, self.B_default)
        c = np.broadcast_to(self.contact_default, (batch, len(self.contact_default))).copy()
        free_c = self.contact_index >= 0
        if free_c.any():
            c[:, free_c] = m[:, self.contact_index[free_c]]
        p["contact_depth"] = c
        return p

    def _remap_negative_float(self, p: dict) -> None:
        """A negative floating-grain fraction is snapped to zero by moving the
        sand properties back along the trend-curve axes (as in Java Delivery),
        so the truncation does not smear the distribution."""
        ff = p["float"]
        for i, sand in self._sand.items():
            neg = ff[:, i] < 0.0
            if not neg.any() or (sand.vp_float_coeff == 0 and sand.phi_float_coeff == 0):
                continue
            vp = p["vp_s"][:, i]
            vp_new = vp - sand.vp_float_coeff * ff[:, i]
            phi_new = p["phi_s"][:, i] + sand.phi_slope * (vp_new - vp) - sand.phi_float_coeff * ff[:, i]
            vs_new = p["vs_s"][:, i] + sand.vs_slope * (vp_new - vp)
            p["vp_s"][:, i] = np.where(neg, vp_new, vp)
            p["phi_s"][:, i] = np.where(neg, phi_new, p["phi_s"][:, i])
            p["vs_s"][:, i] = np.where(neg, vs_new, p["vs_s"][:, i])
        np.maximum(ff, 0.0, out=ff)
        np.minimum(ff, 1.0, out=ff)

    def order_times(self, t: np.ndarray) -> np.ndarray:
        """Enforce t[i] <= t[i+1] by truncation.

        Times are fixed one at a time, best-known (smallest prior sigma)
        first; each later one is clipped between the already-fixed times
        above and below it. A layer whose top is clipped onto its base has
        zero thickness: it is pinched out.
        """
        out = t.copy()
        n = t.shape[1]
        lo = np.full_like(t, -np.inf)
        hi = np.full_like(t, np.inf)
        for k in self.time_order:
            v = np.clip(out[:, k], lo[:, k], hi[:, k])
            out[:, k] = v
            if k + 1 < n:
                lo[:, k + 1:] = np.maximum(lo[:, k + 1:], v[:, None])
            if k > 0:
                hi[:, :k] = np.minimum(hi[:, :k], v[:, None])
        return out

    # ------------------------------------------------------------------
    def _fluid_props(self, p: dict, kind: str):
        """(S, vp, rho, has_hc) arrays of shape (batch, nl) for one fluid kind."""
        if kind == "brine":
            z = np.zeros_like(p["NG"])
            return z, np.ones_like(z), np.ones_like(z), np.zeros(self.nl, dtype=bool)
        return p[f"S_{kind}"], p[f"vp_{kind}"], p[f"rho_{kind}"], self.perm.copy()

    def physical(self, m: np.ndarray) -> dict[str, np.ndarray]:
        """Property tables plus effective properties, interval structure,
        depths and validity."""
        p = self.tables(m)
        batch, nl = p["NG"].shape
        vp_i, vs_i, rho_i = impermeable_properties(
            p["vp_m"], p["vs_m"], p["rho_m"], p["NG_hard"], p["vp_hard"], p["vs_hard"],
            p["rho_hard"], self.has_hard)
        p.update(vp_imp=vp_i, vs_imp=vs_i, rho_imp=rho_i)

        def eff(kind_by_layer):
            S = np.zeros((batch, nl))
            vph = np.ones((batch, nl))
            rhoh = np.ones((batch, nl))
            has = np.zeros(nl, dtype=bool)
            for i, kind in enumerate(kind_by_layer):
                if kind != "brine" and self.perm[i]:
                    S[:, i], vph[:, i], rhoh[:, i] = p[f"S_{kind}"][:, i], p[f"vp_{kind}"][:, i], p[f"rho_{kind}"][:, i]
                    has[i] = True
            vp, vs, rho, ok = effective_properties(
                p["NG"], p["phi_s"], p["vp_s"], p["vs_s"], vp_i, vs_i, rho_i,
                p["vp_b"], p["rho_b"], S, vph, rhoh, has,
                self.K_g, self.rho_g, self.K_b_ref, self.rho_b_ref)
            return vp, vs, rho, ok, np.where(has, S, 0.0)

        t = p["t"]
        if not self.contacts_mode:
            vp, vs, rho, ok, S = eff(self.fluid_kind)
            p.update(vp_eff=vp, vs_eff=vs, rho_eff=rho, S_h=S)
            p["valid"] = np.all(ok | ~self.active[None, :], axis=1)
            ft = np.broadcast_to(np.array(self.state.types), (batch, nl))
            p["fluid_type"] = ft.copy()
            for h, code in (("oil", 1), ("lsg", 2), ("gas", 3)):
                p[f"S_eff_{h}"] = np.where(ft == code, S, 0.0)
            p["z"], p["z_base"] = self._hang_depths(p["d"], t, vp)
            # one interval per layer
            p["iv_t"], p["iv_present"] = t[:, :nl], t[:, 1:] > t[:, :nl]
            p["iv_vp"], p["iv_vs"], p["iv_rho"] = vp, vs, rho
            p["iv_layer"] = np.arange(nl)
        else:
            self._contacts_physical(p, eff)
        return p

    def _hang_depths(self, d, t, vp):
        """Top depth of each layer, hung from the reference layer's depth."""
        nl = self.nl
        k = self.cfg.master_depth_layer
        dz = 0.5 * (t[:, 1:] - t[:, :-1]) * vp
        z = np.zeros((d.shape[0], nl + 1))
        z[:, k] = d[:, k]
        for i in range(k, nl):
            z[:, i + 1] = z[:, i] + dz[:, i]
        for i in range(k - 1, -1, -1):
            z[:, i] = z[:, i + 1] - dz[:, i]
        return z[:, :nl], z[:, nl]

    def _contacts_physical(self, p: dict, eff) -> None:
        cfg = self.cfg
        nl = self.nl
        t = p["t"]
        batch = t.shape[0]
        # effective properties of every layer for each possible fluid
        props = [eff([FLUID_NAMES[f]] * nl) for f in range(4)]
        vp4 = np.stack([q[0] for q in props], axis=-1)  # (batch, nl, 4)
        vs4 = np.stack([q[1] for q in props], axis=-1)
        rho4 = np.stack([q[2] for q in props], axis=-1)
        ok4 = np.stack([q[3] for q in props], axis=-1)
        S4 = np.stack([q[4] for q in props], axis=-1)
        cd = p["contact_depth"]
        units = cfg.fluid_units
        c0 = np.cumsum([0] + [len(u.contacts) for u in units])

        K = MAX_SUBLAYERS
        sub_t0 = np.zeros((batch, nl, K))
        sub_t1 = np.zeros((batch, nl, K))
        sub_f = np.zeros((batch, nl, K), dtype=int)
        sub_on = np.zeros((batch, nl, K), dtype=bool)
        valid = np.ones(batch, dtype=bool)
        z_top = np.zeros((batch, nl))
        z_bot = np.zeros((batch, nl))
        rows = np.arange(batch)
        k_ref = cfg.master_depth_layer
        order = list(range(k_ref, nl)) + list(range(k_ref - 1, -1, -1))
        for i in order:
            u_id = next(n for n, u in enumerate(units) if u.owns(i))
            u = units[u_id]
            ucd = cd[:, c0[u_id]:c0[u_id + 1]]
            UF = np.array([c.upper for c in u.contacts])
            LF = np.array([c.lower for c in u.contacts])
            down = i >= k_ref
            if down:
                z = p["d"][:, i] if i == k_ref else z_bot[:, i - 1]
                tf = t[:, i].copy()
            else:
                z = z_top[:, i + 1]
                tf = t[:, i + 1].copy()
            z_start = z.copy()
            running = np.ones(batch, dtype=bool)
            for s in range(K):
                f, lo, hi, ok = _fluid_at(z, ucd, UF, LF, down)
                valid &= ok | ~running
                f = np.where(ok, f, 0)
                vp = vp4[rows, i, f]
                if down:
                    z_end = z + 0.5 * (t[:, i + 1] - tf) * vp
                    done = (z_end < hi) | (s == K - 1)
                    tc = np.where(done, t[:, i + 1], tf + 2.0 * (hi - z) / vp)
                    seg0, seg1 = tf, tc
                    z_next = np.where(done, z_end, hi)
                else:
                    z_end = z - 0.5 * (tf - t[:, i]) * vp
                    done = (z_end > lo) | (s == K - 1)
                    tc = np.where(done, t[:, i], tf - 2.0 * (z - lo) / vp)
                    seg0, seg1 = tc, tf
                    z_next = np.where(done, z_end, lo)
                sub_t0[:, i, s] = np.where(running, seg0, 0.0)
                sub_t1[:, i, s] = np.where(running, seg1, 0.0)
                sub_f[:, i, s] = np.where(running, f, 0)
                sub_on[:, i, s] = running
                valid &= ok4[rows, i, f] | ~running | ~self.active[i]
                z = np.where(running, z_next, z)
                tf = np.where(running, tc, tf)
                running &= ~done
                if not running.any():
                    break
            if down:
                z_top[:, i], z_bot[:, i] = z_start, z
            else:
                z_bot[:, i], z_top[:, i] = z_start, z

        # sort sub-layers of each layer into time order (ascending walks run upward)
        key = np.where(sub_on, sub_t0, np.inf)
        o = np.argsort(key, axis=2, kind="stable")
        take = lambda a: np.take_along_axis(a, o, axis=2)  # noqa: E731
        sub_t0, sub_t1, sub_f, sub_on = take(sub_t0), take(sub_t1), take(sub_f), take(sub_on)
        thick = np.where(sub_on, np.maximum(sub_t1 - sub_t0, 0.0), 0.0)
        bi, li = np.meshgrid(np.arange(batch), np.arange(nl), indexing="ij")
        g = lambda a4: a4[bi[..., None], li[..., None], sub_f]  # noqa: E731
        s_vp, s_vs, s_rho, s_S = g(vp4), g(vs4), g(rho4), g(S4)
        # upscaled layer properties (paper Appendix 11.1)
        tt = thick.sum(axis=2)
        zt = (thick * s_vp).sum(axis=2)
        safe = tt > 0
        with np.errstate(divide="ignore", invalid="ignore"):
            vp_eff = np.where(safe, zt / tt, s_vp[..., 0])
            vs_eff = np.where(safe, (thick * s_vs).sum(axis=2) / tt, s_vs[..., 0])
            rho_eff = np.where(safe, (thick * s_vp * s_rho).sum(axis=2) / zt, s_rho[..., 0])
            S_h = np.where(safe, (thick * s_vp * s_S).sum(axis=2) / zt, s_S[..., 0])
            for h, code in (("oil", 1), ("lsg", 2), ("gas", 3)):
                w = np.where(sub_f == code, thick * s_vp * s_S, 0.0).sum(axis=2)
                p[f"S_eff_{h}"] = np.where(safe, w / zt, np.where(sub_f[..., 0] == code, s_S[..., 0], 0.0))
        # effective fluid: the thickest non-brine sub-layer, else brine
        zthick = np.where(sub_f > 0, thick * s_vp, -1.0)
        best = np.argmax(zthick, axis=2)
        ft = np.take_along_axis(sub_f, best[..., None], axis=2)[..., 0]
        ft = np.where(np.max(zthick, axis=2) >= 0, ft, sub_f[..., 0])
        ft = np.where(self.perm[None, :], ft, 0)
        p.update(vp_eff=vp_eff, vs_eff=vs_eff, rho_eff=rho_eff, S_h=S_h, fluid_type=ft,
                 valid=valid, z=z_top, z_base=z_bot[:, -1], n_sublayers=sub_on.sum(axis=2),
                 sub_t0=sub_t0, sub_t1=sub_t1, sub_f=sub_f)
        # flattened intervals for the reflectivity
        p["iv_t"] = sub_t0.reshape(batch, -1)
        p["iv_present"] = (thick > 0).reshape(batch, -1)
        p["iv_vp"], p["iv_vs"], p["iv_rho"] = (s_vp.reshape(batch, -1), s_vs.reshape(batch, -1),
                                               s_rho.reshape(batch, -1))
        p["iv_layer"] = np.repeat(np.arange(nl), K)

    # ------------------------------------------------------------------
    def interval_reflections(self, p: dict, stack_index: int) -> np.ndarray:
        """Reflection coefficient at the top of every interval (layer or
        contact sub-layer) of non-zero thickness: the contrast with the
        nearest non-zero-thickness interval above. Pinched-out intervals are
        skipped. Shape (batch, n_intervals)."""
        sm = self.seismic.stacks[stack_index]
        present = p["iv_present"]
        batch, n = present.shape
        rc = np.zeros((batch, n))
        have = np.zeros(batch, dtype=bool)
        a_vp = np.ones(batch)
        a_vs = np.zeros(batch)
        a_rho = np.ones(batch)
        for k in range(n):
            pk = present[:, k]
            vp, vs, rho = p["iv_vp"][:, k], p["iv_vs"][:, k], p["iv_rho"][:, k]
            both = pk & have
            if both.any():
                r = reflection_coefficient(a_rho, a_vp, a_vs, rho, vp, vs, p["A"], sm.theta_sq(a_vp, p["B"]))
                rc[:, k] = np.where(both, r, 0.0)
            a_vp = np.where(pk, vp, a_vp)
            a_vs = np.where(pk, vs, a_vs)
            a_rho = np.where(pk, rho, a_rho)
            have |= pk
        return rc

    def layer_reflections(self, p: dict, stack_index: int) -> tuple[np.ndarray, np.ndarray]:
        """(reflection at each layer top, (batch, nl)) and (reflections at
        contacts inside layers, (batch, nl, MAX_SUBLAYERS - 1))."""
        rc = self.interval_reflections(p, stack_index)
        if not self.contacts_mode:
            return rc, np.zeros((rc.shape[0], self.nl, 0))
        rc = rc.reshape(rc.shape[0], self.nl, MAX_SUBLAYERS)
        present = p["iv_present"].reshape(rc.shape)
        # the first present sub-layer of a layer carries the layer-top reflection
        first = np.argmax(present, axis=2)
        top = np.take_along_axis(rc, first[..., None], axis=2)[..., 0]
        internal = rc.copy()
        np.put_along_axis(internal, first[..., None], 0.0, axis=2)
        return top, internal[..., 1:]

    def evaluate(self, m: np.ndarray, ignore_seismic: bool = False) -> Evaluation:
        """Residuals of a batch of model vectors, shape (batch, dim)."""
        m = np.atleast_2d(m)
        p = self.physical(m)
        parts = [m @ self.Jp.T - self.cp]
        n_seis = 0
        if self.seismic is not None:
            for k, sm in enumerate(self.seismic.stacks):
                rc = self.interval_reflections(p, k)
                R = self.seismic.reflectivity(p["iv_t"], rc)
                res = (sm.observed - R @ sm.G.T) * sm.scale * sm.use
                parts.append(np.zeros_like(res) if ignore_seismic else res)
                n_seis += res.shape[1]
        t = p["t"]
        thick = 0.5 * (t[:, 1:] - t[:, :-1]) * p["vp_eff"]
        for i, dz, sigma in self.isopachs:
            parts.append(((thick[:, i] - dz) / sigma)[:, None])
        for b, tl, sigma in self.thickness_couplings:
            parts.append(((thick[:, b] - thick[:, tl]) / sigma)[:, None])
        res = np.nan_to_num(np.concatenate(parts, axis=1), nan=1e6, posinf=1e6, neginf=-1e6)
        return Evaluation(res, p["valid"], self.Jp.shape[0], n_seis)

    def time_traces(self, m: np.ndarray, times: np.ndarray) -> dict[str, np.ndarray]:
        """Impedance and synthetic seismic sampled at ``times`` (seconds) for a
        batch of model vectors: the traces used in a standard inversion QC.

        Returns ``Zp`` and ``Zs`` (batch, n_times): P and S impedance of the
        layer (or contact sub-layer) at each time, with the top layer extended
        upwards and the bottom layer downwards; and ``synthetic_<k>`` per stack
        (zero outside the modelled window).
        """
        m = np.atleast_2d(m)
        times = np.asarray(times, dtype=float)
        p = self.physical(m)
        Zp_iv = p["iv_vp"] * p["iv_rho"]
        Zs_iv = p["iv_vs"] * p["iv_rho"]
        present = p["iv_present"]
        first = np.argmax(present, axis=1)
        rows = np.arange(len(m))
        Zp = np.repeat(Zp_iv[rows, first][:, None], len(times), axis=1)
        Zs = np.repeat(Zs_iv[rows, first][:, None], len(times), axis=1)
        for k in range(present.shape[1]):
            below = present[:, k, None] & (times[None, :] >= p["iv_t"][:, k, None])
            Zp = np.where(below, Zp_iv[:, k, None], Zp)
            Zs = np.where(below, Zs_iv[:, k, None], Zs)
        out = {"Zp": Zp, "Zs": Zs}
        if self.seismic is not None:
            grid = self.seismic.t_start + self.seismic.dt * np.arange(self.seismic.ns)
            for k, syn in enumerate(self.synthetics(m)):
                out[f"synthetic_{k}"] = np.stack(
                    [np.interp(times, grid, row, left=0.0, right=0.0) for row in syn])
        return out

    def synthetics(self, m: np.ndarray) -> list[np.ndarray]:
        """Synthetic trace(s) on the model time grid, one array per stack."""
        p = self.physical(m)
        out = []
        for k in range(len(self.seismic.stacks)):
            R = self.seismic.reflectivity(p["iv_t"], self.interval_reflections(p, k))
            out.append(self.seismic.synthetic(k, R))
        return out


def _pick(override, default):
    return default if override is None else override


def _fluid_at(z, ucd, UF, LF, descending):
    """Fluid at depth ``z`` within one fluid unit, and the depth bracket of
    contacts around it. Returns (fluid code, lo, hi, ok)."""
    left = (ucd <= z[:, None]) if descending else (ucd < z[:, None])
    right = (z[:, None] < ucd) if descending else (z[:, None] <= ucd)
    n = ucd.shape[1]
    has_l = left.any(axis=1)
    has_r = right.any(axis=1)
    kl = n - 1 - np.argmax(left[:, ::-1], axis=1)
    kr = np.argmax(right, axis=1)
    f_l = LF[kl]
    f_r = UF[kr]
    f = np.where(has_l, f_l, f_r)
    ok = ~(has_l & has_r & (f_l != f_r))
    rows = np.arange(len(z))
    lo = np.where(has_l, ucd[rows, kl], -np.inf)
    hi = np.where(has_r, ucd[rows, kr], np.inf)
    return f, lo, hi, ok
