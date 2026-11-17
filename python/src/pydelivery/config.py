"""Parse a Delivery ``ModelDescription.xml`` file into plain dataclasses.

The XML layout is the one defined by ``schemas/Delivery.xsd`` in the original
Java distribution.  Times in the XML are in milliseconds and are converted to
seconds here; everything else is kept in the user's units (the examples use
SI).

Features of the Java code that this port does not implement yet raise
``NotImplementedError`` when the XML asks for them, rather than being
silently ignored.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

BAD_VALUE = -999.0
TIME_UNITS_AUTO_DETECTION_MAX = 50.0  # prior-trace times above this are ms

HYDROCARBONS = ("oil", "lsg", "gas")


@dataclass
class Fluid:
    name: str
    vp: float
    sigma_vp: float
    rho: float
    sigma_rho: float


@dataclass
class RockTrend:
    """Regression ("trend") curves for one end-member rock.

    vp  = vp_intercept + vp_depth_coeff*depth + vp_lfiv_coeff*LFIV  +/- vp_sigma
    vs  = vs_intercept + vs_slope*vp                                +/- vs_sigma
    reservoir rocks:     phi = phi_intercept + phi_slope*vp         +/- phi_sigma
    non-reservoir rocks: rho = rho_factor * vp**rho_exponent        +/- rho_sigma
    """

    name: str
    is_reservoir: bool
    vp_intercept: float
    vp_depth_coeff: float
    vp_lfiv_coeff: float
    vp_sigma: float
    vs_intercept: float
    vs_slope: float
    vs_sigma: float
    # reservoir end-member
    phi_intercept: float = 0.0
    phi_slope: float = 0.0
    phi_sigma: float = 0.0
    vp_grain: float = 0.0
    vs_grain: float = 0.0
    rho_grain: float = 0.0
    reference_fluid: str | None = None
    # non-reservoir end-member (Gardner type power law)
    rho_factor: float = 0.0
    rho_exponent: float = 0.0
    rho_sigma: float = 0.0

    @property
    def K_grain(self) -> float:
        return self.rho_grain * (self.vp_grain**2 - (4.0 / 3.0) * self.vs_grain**2)


@dataclass
class Hydrocarbon:
    name: str
    probability: float = 0.0
    saturation: float = 0.0
    sigma_saturation: float = 0.0


@dataclass
class Layer:
    name: str
    shale: str
    sand: str | None = None
    brine: str | None = None
    hydrocarbons: dict[str, Hydrocarbon] = field(default_factory=dict)
    time: float = 0.0  # seconds
    sigma_time: float = 0.0
    ng: float = 0.0
    sigma_ng: float = 0.0
    depth: float = 0.0
    sigma_depth: float = 0.0
    lfiv: float = 0.0
    sigma_lfiv: float = 0.0
    thickness: float | None = None  # isopach constraint (depth units)
    sigma_thickness: float = 0.0
    depth_master: int = -1  # layer index owning the depth variable
    lfiv_master: int = -1

    @property
    def is_permeable(self) -> bool:
        return self.sand is not None and (self.ng > 0 or self.sigma_ng > 0)


@dataclass
class Stack:
    name: str
    filename: str
    wavelet_filename: str
    min_offset: float
    max_offset: float
    reflector_time: float  # seconds
    stack_velocity: float
    noise_rms: float
    time_above: float | None = None  # seconds
    time_below: float | None = None

    @property
    def offset_mean_square(self) -> float:
        a, b = self.min_offset, self.max_offset
        if a == b:
            return a * a
        return (b**3 - a**3) / (3.0 * (b - a))


@dataclass
class ModelConfig:
    layers: list[Layer]
    t_base: float
    sigma_t_base: float
    stacks: list[Stack]
    fluids: dict[str, Fluid]
    rocks: dict[str, RockTrend]
    A: float = 1.0
    sigma_A: float = 0.0
    B: float = 1.0
    sigma_B: float = 0.0
    density_ordering: str = "full"  # none | partial | full
    links: dict[str, int] = field(default_factory=dict)  # tag -> 1-based block
    varying: list[tuple[int, str]] = field(default_factory=list)  # (layer, tag)
    base_dir: Path = Path(".")

    @property
    def n_layers(self) -> int:
        return len(self.layers)

    def path(self, filename: str) -> Path:
        return self.base_dir / filename


# --------------------------------------------------------------------------
# XML helpers


def _num(elem: ET.Element | None, tag: str, default: float | None = None) -> float:
    e = elem.find(tag) if elem is not None else None
    if e is None or e.text is None or not e.text.strip():
        if default is None:
            raise ValueError(f"missing <{tag}> in <{elem.tag if elem is not None else '?'}>")
        return default
    return float(e.text)


def _text(elem: ET.Element, tag: str, default: str | None = None) -> str | None:
    e = elem.find(tag)
    return e.text.strip() if e is not None and e.text else default


def _flag(elem: ET.Element | None, attr: str) -> bool:
    return elem is not None and elem.get(attr, "false").strip().lower() == "true"


# layer-level tags that may vary areally: xml tag -> (attribute, scale)
_LAYER_TAGS = {
    "time": ("time", 1e-3),
    "sigma_time": ("sigma_time", 1e-3),
    "net_to_gross": ("ng", 1.0),
    "sigma_net_to_gross": ("sigma_ng", 1.0),
    "depth": ("depth", 1.0),
    "sigma_depth": ("sigma_depth", 1.0),
    "LFIV": ("lfiv", 1.0),
    "sigma_LFIV": ("sigma_lfiv", 1.0),
    "thickness": ("thickness", 1.0),
    "sigma_thickness": ("sigma_thickness", 1.0),
}
_HC_TAGS = {}
for _h in HYDROCARBONS:
    _HC_TAGS[f"probability_of_{_h}"] = (_h, "probability")
    _HC_TAGS[f"saturation_of_{_h}"] = (_h, "saturation")
    _HC_TAGS[f"sigma_saturation_of_{_h}"] = (_h, "sigma_saturation")

_UNSUPPORTED = (
    "fluid_units_and_contacts",
    "stacked_models_info",
    "nonreservoir_hard_endmember",
    "fluid_dependent_nonreservoir_hard_endmembers",
    "seismic_likelihood_exclusion_range",
    "layer_mappings",
)


def load_config(xml_path: str | Path) -> ModelConfig:
    xml_path = Path(xml_path)
    root = ET.parse(xml_path).getroot()

    for tag in _UNSUPPORTED:
        if root.find(f".//{tag}") is not None:
            raise NotImplementedError(f"<{tag}> is not supported by pydelivery yet")
    for e in root.iter("float_fraction"):
        if e.text and float(e.text) != 0.0:
            raise NotImplementedError("floating-grain models are not supported yet")
    if (_text(root, ".//take_fluid_probabilities_as_post_ordering") or "false").lower() == "true":
        raise NotImplementedError("post-ordering fluid probabilities are not supported yet")

    # ---- seismic stacks -------------------------------------------------
    sd = root.find("seismic_data")
    stacks = []
    for s in sd.findall("stack"):
        w = s.find("wavelet")
        if w.find("noise_rms") is None:
            raise NotImplementedError("legacy <signal_to_noise> wavelet entries are not supported")
        ta = s.find("time_above") if s.find("time_above") is not None else w.find("time_above")
        tb = s.find("time_below") if s.find("time_below") is not None else w.find("time_below")
        stacks.append(
            Stack(
                name=_text(s, "name", "stack"),
                filename=_text(s, "filename"),
                wavelet_filename=_text(w, "filename"),
                min_offset=_num(s, "min_offset"),
                max_offset=_num(s, "max_offset"),
                reflector_time=1e-3 * _num(s, "reflector_time"),
                stack_velocity=_num(s, "stack_velocity"),
                noise_rms=_num(w, "noise_rms"),
                time_above=None if ta is None else 1e-3 * float(ta.text),
                time_below=None if tb is None else 1e-3 * float(tb.text),
            )
        )
    if not 1 <= len(stacks) <= 2:
        raise ValueError("expected one or two <stack> entries")
    avo = sd.find("AVO_terms")
    A = _num(avo, "A", 1.0)
    B = _num(avo, "B", 1.0)
    sigma_A = _num(avo, "sigma_A", 0.0)
    sigma_B = _num(avo, "sigma_B", 0.0)

    # ---- fluids and rocks -----------------------------------------------
    rfp = root.find("rock_fluid_properties")
    fluids = {}
    for f in rfp.find("fluid_properties").findall("fluid"):
        name = _text(f, "name")
        fluids[name] = Fluid(
            name, _num(f, "vp"), _num(f, "sigma_vp", 0.0), _num(f, "rho"), _num(f, "sigma_rho", 0.0)
        )
    rocks = {}
    rp = rfp.find("rock_properties")
    for r in rp.findall("reservoir_endmember") + rp.findall("nonreservoir_endmember"):
        name = _text(r, "name")
        vp, vs = r.find("vp_curve"), r.find("vs_curve")
        if _num(vp, "phi_float_coefficient", 0.0) != 0.0:
            raise NotImplementedError("floating-grain models are not supported yet")
        common = dict(
            name=name,
            vp_intercept=_num(vp, "intercept"),
            vp_depth_coeff=_num(vp, "depth_coefficient", 0.0),
            vp_lfiv_coeff=_num(vp, "LFIV_coefficient", 0.0),
            vp_sigma=_num(vp, "sigma"),
            vs_intercept=_num(vs, "intercept"),
            vs_slope=_num(vs, "slope"),
            vs_sigma=_num(vs, "sigma"),
        )
        if r.tag == "reservoir_endmember":
            ph = r.find("porosity_curve")
            rocks[name] = RockTrend(
                is_reservoir=True,
                phi_intercept=_num(ph, "intercept"),
                phi_slope=_num(ph, "slope"),
                phi_sigma=_num(ph, "sigma"),
                vp_grain=_num(r, "vp_grain"),
                vs_grain=_num(r, "vs_grain"),
                rho_grain=_num(r, "rho_grain"),
                reference_fluid=_text(r, "reference_fluid"),
                **common,
            )
        else:
            dc = r.find("density_curve")
            rocks[name] = RockTrend(
                is_reservoir=False,
                rho_factor=_num(dc, "factor"),
                rho_exponent=_num(dc, "exponent"),
                rho_sigma=_num(dc, "sigma"),
                **common,
            )

    # ---- layers ---------------------------------------------------------
    md = root.find("model_description")
    layers: list[Layer] = []
    varying: list[tuple[int, str]] = []
    t_base = sigma_t_base = None
    for i, le in enumerate(list(md)):
        if le.tag not in ("top_layer", "middle_layer", "bottom_layer"):
            continue
        if _flag(le, "synchronise_rock_matrix"):
            raise NotImplementedError("synchronise_rock_matrix is not supported yet")
        res = le.find("reservoir_endmember")
        nonres = le.find("nonreservoir_endmember")
        lay = Layer(name=_text(le, "name", f"layer{i + 1}"), shale=_text(nonres, "name"))
        if res is not None:
            lay.sand = _text(res, "name")
            b = res.find("brine")
            lay.brine = _text(b, "name") if b is not None else None
            for h in HYDROCARBONS:
                he = res.find(h)
                if he is not None:
                    lay.hydrocarbons[h] = Hydrocarbon(
                        name=_text(he, "name"),
                        probability=_num(he, f"probability_of_{h}", 0.0),
                        saturation=_num(he, f"saturation_of_{h}", 0.0),
                        sigma_saturation=_num(he, f"sigma_saturation_of_{h}", 0.0),
                    )
        idx = len(layers)
        for tag, (attr, scale) in _LAYER_TAGS.items():
            e = le.find(tag)
            if e is not None and e.text and e.text.strip():
                setattr(lay, attr, scale * float(e.text))
                if _flag(e, "varies_areally"):
                    varying.append((idx, tag))
        for tag in _HC_TAGS:
            e = le.find(f".//{tag}")
            if e is not None and _flag(e, "varies_areally"):
                varying.append((idx, tag))
        # shared ("same as layer above") depth and LFIV variables
        lay.depth_master = idx
        lay.lfiv_master = idx
        if idx > 0:
            if _flag(le.find("depth"), "same_as_layer_above"):
                lay.depth_master = layers[-1].depth_master
                lay.depth, lay.sigma_depth = layers[-1].depth, layers[-1].sigma_depth
            if _flag(le.find("LFIV"), "same_as_layer_above"):
                lay.lfiv_master = layers[-1].lfiv_master
                lay.lfiv, lay.sigma_lfiv = layers[-1].lfiv, layers[-1].sigma_lfiv
        for tag in ("time", "net_to_gross"):
            if _flag(le.find(tag), "same_as_layer_above"):
                raise NotImplementedError(f'same_as_layer_above="true" on <{tag}> is not supported yet')
        if lay.thickness is not None and lay.thickness == 0 and lay.sigma_thickness == 0:
            raise NotImplementedError("forced pinch-out via thickness N(0,0) is not supported yet")
        if lay.sigma_time <= 0:
            raise ValueError(f"layer {idx + 1}: sigma_time must be > 0")
        if le.tag == "bottom_layer":
            t_base = 1e-3 * _num(le, "t_base")
            sigma_t_base = 1e-3 * _num(le, "sigma_t_base", 0.0)
        layers.append(lay)
    if t_base is None:
        raise ValueError("no <bottom_layer> with <t_base> found")

    out = root.find("output")
    ordering = (_text(out, "density_ordering", "full") if out is not None else "full").lower()
    if ordering not in ("none", "partial", "full"):
        raise ValueError(f"unknown density_ordering '{ordering}'")
    links = {}
    pi = root.find("property_indices")
    if pi is not None:
        for e in pi:
            if e.tag.endswith("_link") and e.text and e.text.strip():
                links[e.tag[: -len("_link")]] = int(e.text)

    cfg = ModelConfig(
        layers=layers,
        t_base=t_base,
        sigma_t_base=sigma_t_base,
        stacks=stacks,
        fluids=fluids,
        rocks=rocks,
        A=A,
        sigma_A=sigma_A,
        B=B,
        sigma_B=sigma_B,
        density_ordering=ordering,
        links=links,
        varying=varying,
        base_dir=xml_path.parent,
    )
    _validate(cfg)
    return cfg


def _validate(cfg: ModelConfig) -> None:
    for i, lay in enumerate(cfg.layers):
        if lay.shale not in cfg.rocks:
            raise ValueError(f"layer {i + 1}: unknown non-reservoir rock '{lay.shale}'")
        if lay.sand is not None and lay.sand not in cfg.rocks:
            raise ValueError(f"layer {i + 1}: unknown reservoir rock '{lay.sand}'")
        if lay.is_permeable:
            ref = cfg.rocks[lay.sand].reference_fluid
            if ref not in cfg.fluids:
                raise ValueError(f"rock '{lay.sand}': unknown reference fluid '{ref}'")
            if lay.brine not in cfg.fluids:
                raise ValueError(f"layer {i + 1}: unknown brine '{lay.brine}'")
            for h in lay.hydrocarbons.values():
                if h.probability > 0 and h.name not in cfg.fluids:
                    raise ValueError(f"layer {i + 1}: unknown fluid '{h.name}'")


def apply_prior_trace(cfg: ModelConfig, data: np.ndarray) -> ModelConfig:
    """Return a copy of ``cfg`` updated with one trace of areally varying priors.

    Property ``k`` (1-based block from ``<property_indices>``) of layer ``l``
    (0-based) is ``data[l + (k-1)*n_layers]``.  The time block is last and has
    one extra value, the base time of the bottom layer.  ``-999`` means "keep
    the XML default".
    """
    out = copy.deepcopy(cfg)
    nl = cfg.n_layers
    data = np.asarray(data, dtype=float)

    def value(tag: str, layer: int) -> float | None:
        k = cfg.links.get(tag)
        if k is None:
            return None
        j = layer + (k - 1) * nl
        if j >= len(data) or data[j] == BAD_VALUE:
            return None
        return abs(float(data[j]))

    def seconds(v: float) -> float:
        return v * 1e-3 if v > TIME_UNITS_AUTO_DETECTION_MAX else v

    for layer, tag in cfg.varying:
        v = value(tag, layer)
        if v is None:
            continue
        lay = out.layers[layer]
        if tag in _LAYER_TAGS:
            attr = _LAYER_TAGS[tag][0]
            setattr(lay, attr, seconds(v) if tag in ("time", "sigma_time") else v)
        else:
            h, attr = _HC_TAGS[tag]
            if h in lay.hydrocarbons:
                setattr(lay.hydrocarbons[h], attr, v)
    if "time" in cfg.links:
        v = value("time", nl)
        if v is not None:
            out.t_base = seconds(v)
    # keep shared depth / LFIV values consistent with their master layer
    for lay in out.layers:
        lay.depth = out.layers[lay.depth_master].depth
        lay.sigma_depth = out.layers[lay.depth_master].sigma_depth
        lay.lfiv = out.layers[lay.lfiv_master].lfiv
        lay.sigma_lfiv = out.layers[lay.lfiv_master].sigma_lfiv
    return out
