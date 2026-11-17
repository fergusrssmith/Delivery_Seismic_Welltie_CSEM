"""Readable tables of everything a Delivery model file specifies.

``describe(cfg)`` returns a dict of pandas DataFrames, one per group, so the
parametrisation of an inversion can be read at a glance instead of from the
XML:

* ``seismic``      - stacks, angles, wavelets, noise, AVO scaling
* ``framework``    - layer times, depths, isopachs (the geometry prior)
* ``mixture``      - which rocks make up each layer, net-to-gross
* ``trends``       - rock-physics trend curves, written as equations
* ``fluids``       - fluid properties
* ``layer_fluids`` - which fluids each layer may hold, probabilities, saturations
* ``fluid_states`` - the discrete fluid combinations and their prior probabilities
* ``contacts``     - fluid units and contacts (contacts mode only)
* ``coupling``     - shared variables and vertical-stacking couplings
* ``per_trace``    - values that vary trace by trace (given a prior table)
* ``implied``      - what the priors imply: mean velocity, density, impedance
                     of each layer for each fluid, and the reflection
                     coefficient at its top

``show(cfg)`` displays them with headings in a notebook.
"""
from __future__ import annotations

import math

import numpy as np

from .config import ModelConfig
from .fluids import FLUID_NAMES, enumerate_fluid_states, layer_fluid_probabilities


def _pm(mean, sigma, fmt="{:g}"):
    if mean is None:
        return "-"
    if not sigma:
        return fmt.format(mean) + " (fixed)"
    return f"{fmt.format(mean)} ± {fmt.format(sigma)}"


def _term(coef, name):
    """' + 0.5·depth' / ' - 0.5·depth' with the sign written properly."""
    return f" {'-' if coef < 0 else '+'} {abs(coef):g}·{name}"


def _sig(v, n=4):
    """Round to n significant figures, as an int when that is exact."""
    if v == 0 or not np.isfinite(v):
        return v
    r = float(f"{v:.{n}g}")
    return int(r) if r == int(r) and abs(r) >= 10 ** (n - 1) else r


def _ms(v):
    return None if v is None else 1000.0 * v


def describe(cfg: ModelConfig, priors: dict | None = None) -> dict:
    """Group the model's parameters into tables (pandas DataFrames)."""
    import pandas as pd

    L = cfg.layers
    nl = cfg.n_layers
    out = {}

    # ---- seismic ------------------------------------------------------------
    rows = []
    for k, s in enumerate(cfg.stacks):
        if s.angle is not None:
            angle = f"{math.degrees(s.angle):.1f}° (fixed)"
        else:
            rms = math.degrees(math.sqrt(s.offset_mean_square) / (s.reflector_time * s.stack_velocity))
            angle = f"≈{rms:.1f}° (from offsets, Dix)"
        rows.append({
            "stack": s.name, "seismic file": s.filename, "wavelet file": s.wavelet_filename,
            "offsets (m)": f"{s.min_offset:g} - {s.max_offset:g}",
            "reflector time (ms)": f"{1000 * s.reflector_time:g}",
            "stack velocity (m/s)": f"{s.stack_velocity:g}",
            "mean incidence angle": angle,
            "noise RMS": f"{s.noise_rms:g}",
            "window above/below (ms)": (f"{_ms(s.time_above) or 'wavelet'} / {_ms(s.time_below) or 'wavelet'}"),
        })
    out["seismic"] = pd.DataFrame(rows).set_index("stack")
    out["seismic"].attrs["note"] = (
        f"AVO intercept scale A = {_pm(cfg.A, cfg.sigma_A)}, gradient scale B = "
        f"{_pm(cfg.B, cfg.sigma_B)}{' (B tied to A)' if cfg.sync_B_from_A else ''}. "
        "Noise RMS sets how closely the synthetic must fit the data: smaller "
        "noise makes the seismic dominate the prior.")

    # ---- framework ------------------------------------------------------------
    rows = []
    for i, lay in enumerate(L):
        row = {
            "layer": f"{i + 1}. {lay.name}",
            "top time (ms)": _pm(_ms(lay.time), _ms(lay.sigma_time)),
            "model depth (m)": _pm(lay.depth, lay.sigma_depth),
            "depth for trends (m)": ("model depth" if lay.depth_rock_curves is None
                                     else f"{lay.depth_rock_curves:g} (fixed)"),
            "LFIV (m/s)": _pm(lay.lfiv, lay.sigma_lfiv),
            "isopach (m)": ("-" if lay.thickness is None else
                            ("forced pinch-out" if not lay.active else _pm(lay.thickness, lay.sigma_thickness))),
        }
        rows.append(row)
    rows.append({"layer": "base of model", "top time (ms)": _pm(_ms(cfg.t_base), _ms(cfg.sigma_t_base)),
                 "model depth (m)": "", "depth for trends (m)": "", "LFIV (m/s)": "", "isopach (m)": ""})
    out["framework"] = pd.DataFrame(rows).set_index("layer")
    out["framework"].attrs["note"] = (
        f"Depths are hung from layer {cfg.master_depth_layer + 1} using the layer velocities. "
        "Layer times are put in order by fixing the best-picked (smallest sigma) first; "
        "a layer whose base would rise above its top gets zero thickness (pinch-out).")

    # ---- rock mixture -------------------------------------------------------------
    rows = []
    for i, lay in enumerate(L):
        hard = "-"
        if lay.hard_rock or lay.hard_rock_by_fluid:
            name = lay.hard_rock or ", ".join(f"{k}: {v}" for k, v in lay.hard_rock_by_fluid.items())
            hard = f"{name}; fraction {_pm(lay.ng_hard, lay.sigma_ng_hard)}"
        rows.append({
            "layer": f"{i + 1}. {lay.name}",
            "permeable (reservoir) rock": lay.sand if lay.is_permeable else "-",
            "net-to-gross": _pm(lay.ng, lay.sigma_ng) if lay.is_permeable else "0 (no reservoir rock)",
            "impermeable rock": lay.shale,
            "second impermeable rock": hard,
            "floating grains": (_pm(lay.float_fraction, lay.sigma_float_fraction)
                                if (lay.float_fraction or lay.sigma_float_fraction) else "-"),
        })
    out["mixture"] = pd.DataFrame(rows).set_index("layer")
    out["mixture"].attrs["note"] = (
        "Each layer is a fine lamination of reservoir and non-reservoir rock; "
        "net-to-gross is the reservoir fraction. Properties are Backus-averaged.")

    # ---- trends -------------------------------------------------------------------
    used = {lay.shale for lay in L} | {lay.sand for lay in L if lay.is_permeable}
    used |= {r for lay in L for r in [lay.hard_rock, *lay.hard_rock_by_fluid.values()] if r}
    rows = []
    for name in sorted(used):
        r = cfg.rocks[name]
        vp = f"{r.vp_intercept:g}"
        if r.vp_depth_coeff:
            vp += _term(r.vp_depth_coeff, "depth")
        if r.vp_lfiv_coeff:
            vp += _term(r.vp_lfiv_coeff, "LFIV")
        if r.vp_float_coeff:
            vp += _term(r.vp_float_coeff, "float")
        vs = f"{r.vs_intercept:g}" + _term(r.vs_slope, "vp")
        if r.is_reservoir:
            third = f"porosity = {r.phi_intercept:g}" + _term(r.phi_slope, "vp")
            if r.phi_float_coeff:
                third += _term(r.phi_float_coeff, "float")
            third += f"  ± {r.phi_sigma:g}"
            extra = (f"grains vp {r.vp_grain:g}, vs {r.vs_grain:g}, rho {r.rho_grain:g}; "
                     f"reference fluid {r.reference_fluid}")
        else:
            third = f"density = {r.rho_factor:g}·vp^{r.rho_exponent:g}  ± {r.rho_sigma:g}"
            extra = "-"
        rows.append({"rock": name, "type": "reservoir" if r.is_reservoir else "non-reservoir",
                     "vp (m/s)": f"{vp}  ± {r.vp_sigma:g}", "vs (m/s)": f"{vs}  ± {r.vs_sigma:g}",
                     "density / porosity": third, "grain & reference fluid": extra})
    out["trends"] = pd.DataFrame(rows).set_index("rock")
    out["trends"].attrs["note"] = (
        "Trend curves are regressions from well logs; the ± values are the "
        "regression scatter and become the prior uncertainty. Reservoir-rock "
        "trends describe the rock filled with its reference fluid (brine); "
        "Gassmann substitution then puts in the actual pore fluid.")

    # ---- fluids --------------------------------------------------------------------
    rows = []
    for f in cfg.fluids.values():
        rows.append({"fluid": f.name, "vp (m/s)": _pm(f.vp, f.sigma_vp), "density": _pm(f.rho, f.sigma_rho),
                     "bulk modulus": f"{f.rho * f.vp ** 2:.3g}"})
    out["fluids"] = pd.DataFrame(rows).set_index("fluid")

    P = layer_fluid_probabilities(cfg)
    rows = []
    for i, lay in enumerate(L):
        if not lay.is_permeable:
            continue
        row = {"layer": f"{i + 1}. {lay.name}", "brine": lay.brine}
        for kind in ("oil", "lsg", "gas"):
            hc = lay.hydrocarbons.get(kind)
            if hc is None:
                row[kind] = "-"
                continue
            prob = "set by contacts" if cfg.contacts_mode else f"P = {P[i, FLUID_NAMES.index(kind)]:.2f}"
            row[kind] = f"{hc.name}: {prob}, saturation {_pm(hc.saturation, hc.sigma_saturation)}"
        if not cfg.contacts_mode:
            row["brine"] = f"{lay.brine}: P = {P[i, 0]:.2f}"
        rows.append(row)
    out["layer_fluids"] = pd.DataFrame(rows).set_index("layer") if rows else pd.DataFrame()
    out["layer_fluids"].attrs["note"] = (
        f"Density ordering: '{cfg.density_ordering}'"
        + (" (probabilities are post-ordering targets)" if cfg.post_ordering_probabilities else "")
        + ". 'lsg' is low-saturation gas.")

    states = enumerate_fluid_states(cfg)
    names = [f"{i + 1}. {lay.name}" for i, lay in enumerate(L)]
    rows = [dict(state=s.label, **{n: FLUID_NAMES[t] for n, t in zip(names, s.types)},
                 **{"prior probability": round(s.probability, 4)}) for s in states]
    out["fluid_states"] = pd.DataFrame(rows).set_index("state")
    out["fluid_states"].attrs["note"] = (
        "In contacts mode there is a single placeholder state; fluids follow from depth."
        if cfg.contacts_mode else
        "Each state is a separate model; the inversion jumps between them and "
        "reports how probable each is after seeing the seismic.")

    # ---- contacts ------------------------------------------------------------------
    if cfg.contacts_mode:
        rows = []
        for u_id, u in enumerate(cfg.fluid_units):
            for c in u.contacts:
                rows.append({"fluid unit": u_id + 1, "layers": f"{u.top + 1}-{u.bottom + 1}",
                             "contact": f"{FLUID_NAMES[c.upper]} over {FLUID_NAMES[c.lower]}",
                             "depth (m)": _pm(c.depth, c.sigma_depth),
                             "fluid properties shared in unit": u.synchronised})
        out["contacts"] = pd.DataFrame(rows)

    # ---- coupling ------------------------------------------------------------------
    rows = []
    for i, lay in enumerate(L):
        for what, master in (("depth", lay.depth_master), ("LFIV", lay.lfiv_master),
                             ("rock matrix", lay.rock_master)):
            if master != i:
                rows.append({"layer": f"{i + 1}. {lay.name}", "shares": what,
                             "with layer": f"{master + 1}. {L[master].name}"})
    if cfg.stacking is not None:
        for b, targets in cfg.stacking.mappings:
            for prop, (sigma, app) in cfg.stacking.sigmas.items():
                if app is None or b in app:
                    for t in targets:
                        rows.append({"layer": f"{t + 1}. {L[t].name}",
                                     "shares": f"{prop} (difference ± {sigma:g})",
                                     "with layer": f"{b + 1}. {L[b].name}"})
    for a, b in cfg.exclusions:
        rows.append({"layer": f"{a}..{b}", "shares": "seismic misfit excluded", "with layer": "-"})
    out["coupling"] = pd.DataFrame(rows)

    # ---- per trace -------------------------------------------------------------------
    if priors:
        rows = []
        for tag, arr in priors.items():
            a = np.asarray(arr, dtype=float)
            a2 = a.reshape(len(a), -1)
            if tag in ("time", "sigma_time"):
                a2 = np.where(a2 < 50, 1000.0 * a2, a2)  # seconds -> ms, as in the XML
            for j in range(a2.shape[1]):
                col = a2[:, j]
                if not np.isfinite(col).any():
                    continue
                if tag == "time":
                    where = "base of model" if j == nl else f"{j + 1}. {L[j].name}"
                elif tag.startswith(("contact", "sigma_contact")):
                    where = f"contact {j + 1}"
                elif a.ndim == 1:
                    where = "all"
                else:
                    where = f"{j + 1}. {L[j].name}"
                rows.append({"prior table entry": tag, "applies to": where,
                             "min": float(f"{np.nanmin(col):.5g}"), "mean": float(f"{np.nanmean(col):.5g}"),
                             "max": float(f"{np.nanmax(col):.5g}"),
                             "traces set": int(np.isfinite(col).sum())})
        out["per_trace"] = pd.DataFrame(rows)
        out["per_trace"].attrs["note"] = "These override the XML values trace by trace (times in ms)."

    out["implied"] = implied_properties(cfg)
    return out


def implied_properties(cfg: ModelConfig):
    """Mean effective properties of each layer under each possible fluid,
    with the reflection coefficient at its top (near stack, against the
    layer above with that layer's most likely fluid). Uses prior means."""
    import pandas as pd

    from .model import StateModel
    from .rockphysics import reflection_coefficient

    import copy

    s = cfg.stacks[0]
    if cfg.contacts_mode:
        # show each fluid a contact can put in each layer, one at a time
        cfg = copy.deepcopy(cfg)
        for u in cfg.fluid_units:
            codes = {c for ct in u.contacts for c in (ct.upper, ct.lower)} - {0}
            for i in range(u.top, u.bottom + 1):
                for kind in (FLUID_NAMES[c] for c in codes):
                    if kind in cfg.layers[i].hydrocarbons:
                        cfg.layers[i].hydrocarbons[kind].probability = 1.0 / (len(codes) + 1)
        cfg.fluid_units = None
        cfg.density_ordering = "none"
    states = enumerate_fluid_states(cfg)
    rows = []
    best = {}
    for st in states:
        model = StateModel(cfg, st, None)
        p = model.physical(model.prior_mean[None])
        for i, lay in enumerate(cfg.layers):
            fluid = FLUID_NAMES[st.types[i]]
            key = (i, fluid)
            if key in best:
                continue
            best[key] = (p["vp_eff"][0, i], p["vs_eff"][0, i], p["rho_eff"][0, i], st.probability)
    for (i, fluid), (vp, vs, rho, _) in sorted(best.items()):
        rows.append({"layer": f"{i + 1}. {cfg.layers[i].name}", "fluid": fluid,
                     "vp (m/s)": round(vp), "vs (m/s)": round(vs), "density": _sig(rho),
                     "P impedance": _sig(vp * rho)})
    df = pd.DataFrame(rows)
    # reflection at each layer top, against the layer above with its most likely fluid
    rcs = []
    for (i, fluid) in sorted(best):
        vp2, vs2, rho2, _ = best[(i, fluid)]
        if i == 0:
            rcs.append(np.nan)
            continue
        above = [v for (j, f), v in best.items() if j == i - 1]
        vp1, vs1, rho1, _ = max(above, key=lambda v: v[3])
        if s.angle is not None:
            th = cfg.B * s.angle ** 2
        else:
            th = cfg.B * vp1 ** 2 / (s.stack_velocity ** 4 * s.reflector_time ** 2 / s.offset_mean_square)
        rcs.append(float(reflection_coefficient(rho1, vp1, vs1, rho2, vp2, vs2, cfg.A, th)))
    df["reflection at top"] = ["-" if not np.isfinite(r) else f"{r:+.3f}" for r in rcs]
    df = df.set_index(["layer", "fluid"])
    df.attrs["note"] = (
        "Prior-mean properties after mixing and fluid substitution. The "
        "reflection column shows how strongly each layer/fluid combination "
        "would show up in the near stack: this is what lets the seismic tell "
        "the cases apart.")
    return df


TITLES = {
    "seismic": "Seismic data and noise",
    "framework": "Layer framework (geometry prior)",
    "mixture": "Rock mixture in each layer",
    "trends": "Rock-physics trend curves",
    "fluids": "Fluid properties",
    "layer_fluids": "Possible fluids per layer",
    "fluid_states": "Fluid states (discrete models)",
    "contacts": "Fluid contacts",
    "coupling": "Shared variables and couplings",
    "per_trace": "Values varying trace by trace",
    "implied": "What the priors imply (prior-mean properties)",
}


def show(cfg: ModelConfig, priors: dict | None = None, groups=None) -> dict:
    """Display the tables with headings (in Jupyter) and return them."""
    from IPython.display import HTML, display

    import pandas as pd

    tables = describe(cfg, priors)
    for key in groups or tables:
        df = tables.get(key)
        if df is None or df.empty:
            continue
        display(HTML(f"<h4>{TITLES.get(key, key)}</h4>"))
        with pd.option_context("display.max_colwidth", None):
            display(HTML(df.to_html()))
        note = df.attrs.get("note")
        if note:
            display(HTML(f"<p style='font-size:90%;color:#555'>{note}</p>"))
    return tables
