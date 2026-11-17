"""Parse a Delivery ``ModelDescription.xml`` into plain dataclasses.

The XML layout is the one defined by ``schemas/Delivery.xsd`` of the Java
distribution. Times in the XML are milliseconds and are converted to seconds
here; everything else stays in the user's units.

Supported controls (Java Delivery equivalents in brackets where named
differently): layers with reservoir / non-reservoir / second "hard"
non-reservoir end-members (fixed or fluid-dependent), fluids and fluid
probabilities with density ordering, post-ordering fluid probabilities,
per-layer brine and hydrocarbon property overrides, floating-grain trend
terms, ``depth_rock_curves``, ``same_as_layer_above`` on depth / rock-curve
depth / LFIV / float fraction, ``synchronise_rock_matrix``, isopachs and forced
pinch-outs (thickness N(0, 0)), AVO terms with ``synchronise_B_from_A``,
constant-angle stacks, internal fluid contacts (``fluid_units_and_contacts``),
vertical model stacking (``stacked_models_info``) and seismic-likelihood
exclusion ranges. PS (converted-wave) stacks are not supported.
"""
from __future__ import annotations

import copy
import math
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np

TIME_UNITS_AUTO_DETECTION_MAX = 50.0  # prior-trace times above this are in ms
HYDROCARBONS = ("oil", "lsg", "gas")
FLUID_CODE = {"B": 0, "O": 1, "L": 2, "G": 3}
DENSITY_PRIORITY = (0, 2, 1, 3)  # brine heaviest, then lsg, oil, gas


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

    vp  = vp_intercept + vp_depth_coeff*depth + vp_lfiv_coeff*LFIV
          + vp_float_coeff*float_fraction                              +/- vp_sigma
    vs  = vs_intercept + vs_slope*vp                                   +/- vs_sigma
    reservoir rocks:     phi = phi_intercept + phi_slope*vp + phi_float_coeff*float_fraction
    non-reservoir rocks: rho = rho_factor * vp**rho_exponent  (linearised about the mean vp)
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
    vp_float_coeff: float = 0.0
    phi_intercept: float = 0.0
    phi_slope: float = 0.0
    phi_sigma: float = 0.0
    phi_float_coeff: float = 0.0
    vp_grain: float = 0.0
    vs_grain: float = 0.0
    rho_grain: float = 0.0
    reference_fluid: str | None = None
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
    # per-layer overrides of the fluid table (None: use the table)
    vp: float | None = None
    sigma_vp: float | None = None
    rho: float | None = None
    sigma_rho: float | None = None


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
    depth_rock_curves: float | None = None  # None: trends use the model depth
    lfiv: float = 0.0
    sigma_lfiv: float = 0.0
    thickness: float | None = None  # isopach (depth units)
    sigma_thickness: float = 0.0
    # second ("hard") impermeable end-member
    hard_rock: str | None = None
    hard_rock_by_fluid: dict[str, str] = field(default_factory=dict)  # 'brine','oil',...
    ng_hard: float = 0.0
    sigma_ng_hard: float = 0.0
    # floating grains
    float_fraction: float = 0.0
    sigma_float_fraction: float = 0.0
    # per-layer brine overrides (None: use the fluid table)
    rho_brine: float | None = None
    sigma_rho_brine: float | None = None
    vp_brine: float | None = None
    sigma_vp_brine: float | None = None
    synchronise: bool = False
    # "same as layer above" flags as given in the XML
    sync_above: dict[str, bool] = field(default_factory=dict)
    # resolved owners of shared variables (0-based layer indices)
    depth_master: int = -1
    rock_curves_master: int = -1
    lfiv_master: int = -1
    float_master: int = -1
    rock_master: int = -1  # synchronise_rock_matrix
    active: bool = True  # False: forced pinch-out

    @property
    def is_permeable(self) -> bool:
        return self.sand is not None and (self.ng > 0 or self.sigma_ng > 0)

    @property
    def is_pure_sand(self) -> bool:
        return self.is_permeable and self.ng == 1.0 and self.sigma_ng == 0.0

    @property
    def has_hard_rock(self) -> bool:
        return ((self.hard_rock is not None or bool(self.hard_rock_by_fluid))
                and not self.is_pure_sand
                and not (self.ng_hard == 0.0 and self.sigma_ng_hard == 0.0))

    def hard_rock_for(self, fluid: str) -> str | None:
        if self.hard_rock is not None:
            return self.hard_rock
        return self.hard_rock_by_fluid.get(fluid)


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
    angle: float | None = None  # constant-angle stack, radians

    @property
    def offset_mean_square(self) -> float:
        a, b = self.min_offset, self.max_offset
        if a == b:
            return a * a
        return (b**3 - a**3) / (3.0 * (b - a))


@dataclass
class Contact:
    upper: int  # fluid codes 0 brine, 1 oil, 2 lsg, 3 gas
    lower: int
    depth: float
    sigma_depth: float = 0.0

    @property
    def label(self) -> str:
        return "bolg"[self.upper] + "bolg"[self.lower]


@dataclass
class FluidUnit:
    top: int  # 0-based layer indices, inclusive
    bottom: int
    contacts: list[Contact]
    synchronised: bool = False

    def owns(self, layer: int) -> bool:
        return self.top <= layer <= self.bottom


@dataclass
class StackingInfo:
    """Vertical model stacking: couple properties of base and target layers."""

    mappings: list[tuple[int, list[int]]]  # (base layer, target layers), 0-based
    # property -> (difference sigma, applicable base layers or None for all)
    sigmas: dict[str, tuple[float, set[int] | None]]


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
    sync_B_from_A: bool = False
    density_ordering: str = "full"
    post_ordering_probabilities: bool = False
    master_depth_layer: int = 0  # 0-based
    fluid_units: list[FluidUnit] | None = None
    stacking: StackingInfo | None = None
    exclusions: list[tuple[int, int]] = field(default_factory=list)  # 1-based, as XML
    links: dict[str, int] = field(default_factory=dict)  # tag -> 1-based block
    varying: set[tuple[str, int]] = field(default_factory=set)  # (tag, layer)
    base_dir: Path = Path(".")
    # experimental lateral coupling: model key -> (target value, sigma), set
    # per trace by pydelivery.lateral
    lateral: dict = field(default_factory=dict)

    @property
    def n_layers(self) -> int:
        return len(self.layers)

    @property
    def contacts_mode(self) -> bool:
        return self.fluid_units is not None

    @property
    def contacts(self) -> list[Contact]:
        return [c for u in (self.fluid_units or []) for c in u.contacts]

    def path(self, filename: str) -> Path:
        p = Path(filename)
        return p if p.is_absolute() else self.base_dir / p

    def fluid_unit_of(self, layer: int) -> FluidUnit | None:
        for u in self.fluid_units or []:
            if u.owns(layer):
                return u
        return None


# --------------------------------------------------------------------------
# XML helpers


def _num(elem, tag: str, default: float | None = None) -> float:
    e = elem.find(tag) if elem is not None else None
    if e is None or e.text is None or not e.text.strip():
        if default is None:
            where = elem.tag if elem is not None else "?"
            raise ValueError(f"missing <{tag}> in <{where}>")
        return default
    return float(e.text.split()[0])


def _opt(elem, tag: str) -> float | None:
    e = elem.find(tag) if elem is not None else None
    if e is None or e.text is None or not e.text.strip():
        return None
    return float(e.text.split()[0])


def _text(elem, tag: str, default: str | None = None) -> str | None:
    e = elem.find(tag) if elem is not None else None
    return e.text.strip() if e is not None and e.text and e.text.strip() else default


def _flag(elem, attr: str) -> bool:
    return elem is not None and elem.get(attr, "false").strip().lower() == "true"


def _first(elem, *tags):
    for t in tags:
        e = elem.find(t)
        if e is not None:
            return e
    return None


# layer-level tags: xml tag -> (attribute, scale)
LAYER_TAGS = {
    "time": ("time", 1e-3),
    "sigma_time": ("sigma_time", 1e-3),
    "net_to_gross": ("ng", 1.0),
    "sigma_net_to_gross": ("sigma_ng", 1.0),
    "net_to_gross_hard": ("ng_hard", 1.0),
    "sigma_net_to_gross_hard": ("sigma_ng_hard", 1.0),
    "float_fraction": ("float_fraction", 1.0),
    "sigma_float_fraction": ("sigma_float_fraction", 1.0),
    "depth": ("depth", 1.0),
    "sigma_depth": ("sigma_depth", 1.0),
    "depth_rock_curves": ("depth_rock_curves", 1.0),
    "LFIV": ("lfiv", 1.0),
    "sigma_LFIV": ("sigma_lfiv", 1.0),
    "thickness": ("thickness", 1.0),
    "sigma_thickness": ("sigma_thickness", 1.0),
    "rho_brine": ("rho_brine", 1.0),
    "sigma_rho_brine": ("sigma_rho_brine", 1.0),
    "vp_brine": ("vp_brine", 1.0),
    "sigma_vp_brine": ("sigma_vp_brine", 1.0),
}
# hydrocarbon tags: xml tag -> (kind, attribute)
HC_TAGS = {}
for _h in HYDROCARBONS:
    HC_TAGS[f"probability_of_{_h}"] = (_h, "probability")
    HC_TAGS[f"saturation_of_{_h}"] = (_h, "saturation")
    HC_TAGS[f"sigma_saturation_of_{_h}"] = (_h, "sigma_saturation")
    HC_TAGS[f"vp_{_h}"] = (_h, "vp")
    HC_TAGS[f"sigma_vp_{_h}"] = (_h, "sigma_vp")
    HC_TAGS[f"rho_{_h}"] = (_h, "rho")
    HC_TAGS[f"sigma_rho_{_h}"] = (_h, "sigma_rho")
SYNC_ABOVE_TAGS = ("depth", "depth_rock_curves", "LFIV", "float_fraction")


def load_config(xml_path) -> ModelConfig:
    xml_path = Path(xml_path)
    root = ET.parse(xml_path).getroot()

    # ---- seismic -----------------------------------------------------------
    sd = root.find("seismic_data")
    stacks = []
    for s in sd.findall("stack"):
        if _flag(s, "is-PS-data"):
            raise NotImplementedError("PS (converted-wave) stacks are not supported")
        w = s.find("wavelet")
        noise = _first(w, "noise_rms", "signal_to_noise",
                       "near_stack_noise_rms", "far_stack_noise_rms",
                       "near_stack_signal_to_noise", "far_stack_signal_to_noise")
        if noise is None:
            raise ValueError("stack wavelet needs <noise_rms>")
        ta = _first(s, "time_above") if s.find("time_above") is not None else w.find("time_above")
        tb = _first(s, "time_below") if s.find("time_below") is not None else w.find("time_below")
        ang = s.get("is-constant-angle-stack-with-angle-degrees")
        angle = None
        if ang is not None and float(ang) >= 0:
            if float(ang) >= 90:
                raise ValueError("constant stack angle must be in [0, 90) degrees")
            angle = math.radians(float(ang))
        stacks.append(Stack(
            name=_text(s, "name", "stack"),
            filename=_text(s, "filename"),
            wavelet_filename=_text(w, "filename"),
            min_offset=_num(s, "min_offset", 0.0),
            max_offset=_num(s, "max_offset", 0.0),
            reflector_time=1e-3 * _num(s, "reflector_time"),
            stack_velocity=_num(s, "stack_velocity"),
            noise_rms=float(noise.text.split()[0]),
            time_above=None if ta is None else 1e-3 * float(ta.text),
            time_below=None if tb is None else 1e-3 * float(tb.text),
            angle=angle,
        ))
    if not 1 <= len(stacks) <= 2:
        raise ValueError("expected one or two <stack> entries")
    if len(stacks) == 2 and stacks[0].offset_mean_square > stacks[1].offset_mean_square:
        stacks.reverse()  # near stack first, as in the Java code
    varying: set[tuple[str, int]] = set()
    for k, s in zip(("near", "far"), sd.findall("stack")):
        w = s.find("wavelet")
        e = _first(w, "noise_rms", "signal_to_noise")
        if _flag(e, "varies_areally"):
            varying.add((f"{k}_stack_noise_rms", 0))
    avo = sd.find("AVO_terms")
    sync_B = _flag(avo, "synchronise_B_from_A")
    A, sigma_A = _num(avo, "A", 1.0), _num(avo, "sigma_A", 0.0)
    B, sigma_B = _num(avo, "B", 1.0), _num(avo, "sigma_B", 0.0)
    if sync_B:
        B, sigma_B = A, sigma_A
    for tag in ("A", "sigma_A", "B", "sigma_B"):
        if _flag(avo.find(tag), "varies_areally"):
            varying.add((tag, 0))

    # ---- fluids and rocks ----------------------------------------------------
    rfp = root.find("rock_fluid_properties")
    fluids = {}
    for f in rfp.find("fluid_properties").findall("fluid"):
        name = _text(f, "name")
        fluids[name] = Fluid(name, _num(f, "vp"), _num(f, "sigma_vp", 0.0),
                             _num(f, "rho"), _num(f, "sigma_rho", 0.0))
    rocks = {}
    rp = rfp.find("rock_properties")
    for r in rp.findall("reservoir_endmember") + rp.findall("nonreservoir_endmember"):
        name = _text(r, "name")
        vp, vs = r.find("vp_curve"), r.find("vs_curve")
        common = dict(
            name=name,
            vp_intercept=_num(vp, "intercept"),
            vp_depth_coeff=_num(vp, "depth_coefficient", 0.0),
            vp_lfiv_coeff=_num(vp, "LFIV_coefficient", 0.0),
            vp_sigma=_num(vp, "sigma"),
            vp_float_coeff=_num(vp, "phi_float_coefficient", 0.0),
            vs_intercept=_num(vs, "intercept"),
            vs_slope=_num(vs, "slope"),
            vs_sigma=_num(vs, "sigma"),
        )
        if r.tag == "reservoir_endmember":
            ph = r.find("porosity_curve")
            rocks[name] = RockTrend(
                is_reservoir=True,
                phi_intercept=_num(ph, "intercept"), phi_slope=_num(ph, "slope"),
                phi_sigma=_num(ph, "sigma"), phi_float_coeff=_num(ph, "phi_float_coefficient", 0.0),
                vp_grain=_num(r, "vp_grain"), vs_grain=_num(r, "vs_grain"),
                rho_grain=_num(r, "rho_grain"), reference_fluid=_text(r, "reference_fluid"),
                **common)
        else:
            dc = r.find("density_curve")
            rocks[name] = RockTrend(
                is_reservoir=False, rho_factor=_num(dc, "factor"),
                rho_exponent=_num(dc, "exponent"), rho_sigma=_num(dc, "sigma"), **common)

    # ---- layers --------------------------------------------------------------
    md = root.find("model_description")
    layers: list[Layer] = []
    t_base = sigma_t_base = None
    for le in md:
        if le.tag not in ("top_layer", "middle_layer", "bottom_layer"):
            continue
        idx = len(layers)
        res = le.find("reservoir_endmember")
        lay = Layer(name=_text(le, "name", f"layer {idx + 1}"),
                    shale=_text(le.find("nonreservoir_endmember"), "name"),
                    synchronise=_flag(le, "synchronise_rock_matrix"))
        if res is not None:
            lay.sand = _text(res, "name")
            b = res.find("brine")
            lay.brine = _text(b, "name") if b is not None else None
            for h in HYDROCARBONS:
                he = res.find(h)
                if he is None:
                    continue
                hc = Hydrocarbon(name=_text(he, "name"))
                for tag, (kind, attr) in HC_TAGS.items():
                    if kind != h:
                        continue
                    e = he.find(tag)
                    if e is not None and e.text and e.text.strip():
                        setattr(hc, attr, float(e.text.split()[0]))
                        if _flag(e, "varies_areally"):
                            varying.add((tag, idx))
                lay.hydrocarbons[h] = hc
        hard = le.find("nonreservoir_hard_endmember")
        if hard is not None:
            fd = hard.find("fluid_dependent_nonreservoir_hard_endmembers")
            if fd is not None:
                for kind in ("brine", "oil", "lsg", "gas"):
                    r = _text(fd, f"{kind}_case")
                    if r is not None:
                        lay.hard_rock_by_fluid[kind] = r
            else:
                lay.hard_rock = _text(hard, "name")
        for tag, (attr, scale) in LAYER_TAGS.items():
            e = le.find(f".//{tag}")
            if e is None or not (e.text and e.text.strip()):
                continue
            setattr(lay, attr, scale * float(e.text.split()[0]))
            if _flag(e, "varies_areally"):
                varying.add((tag, idx))
            if _flag(e, "same_as_layer_above"):
                if tag in SYNC_ABOVE_TAGS:
                    lay.sync_above[tag] = idx > 0
                else:
                    warnings.warn(f'layer {idx + 1}: same_as_layer_above on <{tag}> is '
                                  "ignored (as in Java Delivery)")
        if lay.sigma_time <= 0:
            raise ValueError(f"layer {idx + 1}: sigma_time must be > 0")
        if le.tag == "bottom_layer":
            t_base = 1e-3 * _num(le, "t_base")
            sigma_t_base = 1e-3 * _num(le, "sigma_t_base", 0.0)
            if _flag(le.find("t_base"), "varies_areally"):
                varying.add(("time", idx + 1))
        layers.append(lay)
    if t_base is None:
        raise ValueError("no <bottom_layer> with <t_base> found")

    # ---- output / options -----------------------------------------------------
    out = root.find("output")
    ordering = (_text(out, "density_ordering", "full") or "full").lower()
    if ordering not in ("none", "partial", "full"):
        raise ValueError(f"unknown density_ordering '{ordering}'")
    post = (_text(root, ".//take_fluid_probabilities_as_post_ordering") or "false").lower() == "true"
    master_depth = int(_num(out, "master_depth_layer_number", 1.0)) - 1 if out is not None else 0
    if not 0 <= master_depth < len(layers):
        raise ValueError("master_depth_layer_number out of range")
    links = {}
    pi = root.find("property_indices")
    if pi is not None:
        for e in pi:
            if e.tag.endswith("_link") and e.text and e.text.strip() and int(e.text) > 0:
                links[e.tag[:-5]] = int(e.text)

    cfg = ModelConfig(
        layers=layers, t_base=t_base, sigma_t_base=sigma_t_base, stacks=stacks,
        fluids=fluids, rocks=rocks, A=A, sigma_A=sigma_A, B=B, sigma_B=sigma_B,
        sync_B_from_A=sync_B, density_ordering=ordering,
        post_ordering_probabilities=post, master_depth_layer=master_depth,
        links=links, varying=varying, base_dir=xml_path.parent)
    cfg.fluid_units = _parse_fluid_units(root, cfg)
    cfg.stacking, cfg.exclusions = _parse_stacking(root, cfg)
    _resolve_sharing(cfg)
    _validate(cfg)
    return cfg


def _parse_fluid_units(root, cfg: ModelConfig) -> list[FluidUnit] | None:
    fuc = root.find("fluid_units_and_contacts")
    if fuc is None:
        return None
    units = []
    n = 0
    for ue in fuc.findall("fluid_unit"):
        contacts = []
        for ce in ue.findall("contact"):
            key = (_text(ce, "type") or "").upper()
            if len(key) != 2 or any(ch not in FLUID_CODE for ch in key):
                raise ValueError(f"illegal contact type '{key}'")
            up, lo = FLUID_CODE[key[0]], FLUID_CODE[key[1]]
            if up == lo or DENSITY_PRIORITY[up] < DENSITY_PRIORITY[lo]:
                raise ValueError(f"contact '{key}' must be a lighter-over-heavier pair from G,O,L,B")
            contacts.append(Contact(up, lo, _num(ce, "depth"), _num(ce, "sigma_depth", 0.0)))
            for tag, link in (("depth", "contact_depth"), ("sigma_depth", "sigma_contact_depth")):
                if _flag(ce.find(tag), "varies_areally"):
                    cfg.varying.add((link, n))
            n += 1
        if not contacts:
            raise ValueError("each fluid_unit needs at least one <contact>")
        units.append(FluidUnit(int(_num(ue, "top_layer")) - 1, int(_num(ue, "bottom_layer")) - 1,
                               contacts, _flag(ue, "fluid_properties_synchronised")))
    units.sort(key=lambda u: (u.top, u.bottom))
    nl = cfg.n_layers
    if units[0].top != 0 or units[-1].bottom != nl - 1:
        raise ValueError("fluid units must span all layers")
    for a, b in zip(units, units[1:]):
        if a.bottom + 1 != b.top:
            raise ValueError("fluid units must be contiguous and not overlap")
        if cfg.layers[a.bottom].is_permeable and cfg.layers[b.top].is_permeable:
            raise ValueError("adjacent fluid units must be separated by an impermeable layer")
    if n > nl:
        raise ValueError("there can be no more fluid contacts than layers")
    return units


def _parse_stacking(root, cfg: ModelConfig):
    exclusions = []
    for ex in root.iter("seismic_likelihood_exclusion_range"):
        s, e = int(_num(ex, "start_layer")), int(_num(ex, "end_layer"))
        if s > cfg.n_layers or e < 1:
            raise ValueError("seismic exclusion layers out of range")
        exclusions.append((s, e))
    smi = root.find("stacked_models_info")
    if smi is None:
        return None, exclusions
    maps = []
    for m in smi.find("layer_mappings").findall("mapping"):
        base, targets = m.text.split(":")
        b = int(base) - 1
        t = [int(x) - 1 for x in targets.split(",")]
        if not (0 <= b < cfg.n_layers and all(0 <= x < cfg.n_layers for x in t)):
            raise ValueError(f"layer mapping '{m.text}' out of range")
        maps.append((b, t))
    if not maps:
        raise ValueError("stacked_models_info needs at least one <mapping>")
    sigmas = {}
    pds = smi.find("property_difference_std_deviations")
    for prop in ("vp_m", "vp_s", "vs_m", "vs_s", "rho_m", "phi_s", "NG", "thickness"):
        e = pds.find(prop) if pds is not None else None
        if e is None or not (e.text and e.text.strip()):
            continue
        sigma = float(e.text.split()[0])
        if sigma <= 0:
            continue
        app = e.get("applicable_layers")
        layers = None
        if app:
            layers = {int(x) - 1 for x in app.replace(":", ",").split(",") if x.strip()}
        sigmas[prop] = (sigma, layers)
    return StackingInfo(maps, sigmas), exclusions


def _resolve_sharing(cfg: ModelConfig) -> None:
    """Work out which layer owns each shared variable."""
    L = cfg.layers
    nl = len(L)
    # synchronise_rock_matrix: later layers point at the earliest matching layer
    rock_master = list(range(nl))
    for i in range(nl - 1, -1, -1):
        if not L[i].synchronise:
            continue
        for j in range(i + 1, nl):
            if not L[j].synchronise or L[i].shale != L[j].shale:
                continue
            if L[i].is_permeable and L[j].is_permeable and L[i].sand == L[j].sand:
                rock_master[j] = i
            elif not L[i].is_permeable and not L[j].is_permeable:
                rock_master[j] = i
    for i, lay in enumerate(L):
        lay.rock_master = rock_master[i]

    def lowest(tag: str, i: int) -> int:
        j = i
        if L[i].sync_above.get(tag):
            j = i - 1
        if rock_master[i] != i:
            j = min(j, rock_master[i])
        return i if j == i else lowest(tag, j)

    for i, lay in enumerate(L):
        lay.depth_master = lowest("depth", i)
        lay.rock_curves_master = lowest("depth_rock_curves", i)
        lay.lfiv_master = lowest("LFIV", i)
        lay.float_master = lowest("float_fraction", i)
        # forced pinch-out: isopach N(0, 0) (the top layer is always active)
        lay.active = i == 0 or not (lay.thickness == 0 and lay.sigma_thickness == 0)
    _propagate_shared(cfg)


def _propagate_shared(cfg: ModelConfig) -> None:
    for lay in cfg.layers:
        m = cfg.layers[lay.depth_master]
        lay.depth, lay.sigma_depth = m.depth, m.sigma_depth
        m = cfg.layers[lay.lfiv_master]
        lay.lfiv, lay.sigma_lfiv = m.lfiv, m.sigma_lfiv
        m = cfg.layers[lay.rock_curves_master]
        if m.depth_rock_curves is not None:
            lay.depth_rock_curves = m.depth_rock_curves
        m = cfg.layers[lay.float_master]
        lay.float_fraction, lay.sigma_float_fraction = m.float_fraction, m.sigma_float_fraction


def _validate(cfg: ModelConfig) -> None:
    for i, lay in enumerate(cfg.layers):
        if lay.shale not in cfg.rocks:
            raise ValueError(f"layer {i + 1}: unknown non-reservoir rock '{lay.shale}'")
        if lay.sand is not None and lay.sand not in cfg.rocks:
            raise ValueError(f"layer {i + 1}: unknown reservoir rock '{lay.sand}'")
        for r in [lay.hard_rock, *lay.hard_rock_by_fluid.values()]:
            if r is not None and r not in cfg.rocks:
                raise ValueError(f"layer {i + 1}: unknown hard rock '{r}'")
        if lay.is_permeable:
            ref = cfg.rocks[lay.sand].reference_fluid
            if ref not in cfg.fluids:
                raise ValueError(f"rock '{lay.sand}': unknown reference fluid '{ref}'")
            if lay.brine not in cfg.fluids:
                raise ValueError(f"layer {i + 1}: unknown brine '{lay.brine}'")
            for h in lay.hydrocarbons.values():
                if h.name not in cfg.fluids and (h.probability > 0 or cfg.contacts_mode):
                    raise ValueError(f"layer {i + 1}: unknown fluid '{h.name}'")
    if cfg.contacts_mode:
        for u in cfg.fluid_units:
            names: dict[str, str] = {}
            for i in range(u.top, u.bottom + 1):
                lay = cfg.layers[i]
                for kind, hc in lay.hydrocarbons.items():
                    if names.setdefault(kind, hc.name) != hc.name:
                        raise ValueError(f"{kind} names differ within the fluid unit of layer {i + 1}")
                if lay.is_permeable:
                    for c in u.contacts:
                        for code in (c.upper, c.lower):
                            kind = ("brine", "oil", "lsg", "gas")[code]
                            if code and kind not in lay.hydrocarbons:
                                raise ValueError(f"layer {i + 1}: contacts need <{kind}> in the "
                                                 "reservoir end-member")


# --------------------------------------------------------------------------
# per-trace priors


def apply_prior_trace(cfg: ModelConfig, table: dict[str, np.ndarray], i: int) -> ModelConfig:
    """Return a copy of ``cfg`` with the priors of trace ``i`` from a prior
    table (see :mod:`pydelivery.io`). NaN entries keep the XML value."""
    out = copy.deepcopy(cfg)
    nl = cfg.n_layers

    def seconds(v: float) -> float:
        return v * 1e-3 if v > TIME_UNITS_AUTO_DETECTION_MAX else v

    for tag, arr in table.items():
        row = np.atleast_1d(np.asarray(arr[i], dtype=float))
        row = np.where(row == -999.0, np.nan, row)
        if tag in LAYER_TAGS:
            attr, _ = LAYER_TAGS[tag]
            for lyr in range(min(nl, len(row))):
                if np.isfinite(row[lyr]):
                    v = abs(row[lyr])
                    setattr(out.layers[lyr], attr, seconds(v) if tag in ("time", "sigma_time") else v)
            if tag == "time" and len(row) > nl and np.isfinite(row[nl]):
                out.t_base = seconds(abs(row[nl]))
        elif tag in HC_TAGS:
            kind, attr = HC_TAGS[tag]
            for lyr in range(min(nl, len(row))):
                hc = out.layers[lyr].hydrocarbons.get(kind)
                if hc is not None and np.isfinite(row[lyr]):
                    setattr(hc, attr, abs(row[lyr]))
        elif tag in ("contact_depth", "sigma_contact_depth"):
            for c, v in zip(out.contacts, row):
                if np.isfinite(v):
                    setattr(c, "depth" if tag == "contact_depth" else "sigma_depth", abs(v))
        elif tag in ("near_stack_noise_rms", "near_stack_signal_to_noise",
                     "far_stack_noise_rms", "far_stack_signal_to_noise"):
            k = 0 if tag.startswith("near") else 1
            if k < len(out.stacks) and np.isfinite(row[0]):
                out.stacks[k].noise_rms = abs(row[0])
        elif tag in ("A", "sigma_A", "B", "sigma_B"):
            if np.isfinite(row[0]):
                setattr(out, tag, abs(row[0]))
        elif tag == "t_base":
            if np.isfinite(row[0]):
                out.t_base = seconds(abs(row[0]))
        elif tag == "sigma_t_base":
            if np.isfinite(row[0]):
                out.sigma_t_base = seconds(abs(row[0]))
        else:
            raise ValueError(f"unknown prior table entry '{tag}'")
    if out.sync_B_from_A:
        out.B, out.sigma_B = out.A, out.sigma_A
    for lyr, lay in enumerate(out.layers):
        lay.active = lyr == 0 or not (lay.thickness == 0 and lay.sigma_thickness == 0)
    _propagate_shared(out)
    return out
