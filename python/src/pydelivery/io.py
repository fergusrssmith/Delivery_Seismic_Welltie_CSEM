"""File input and output.

Three formats are supported for every kind of data:

* **HDF5** (``.h5``, ``.hdf5``): the main format. Needs ``h5py``.
* **NumPy** (``.npz``): no dependency beyond NumPy.
* **SEG-Y** (``.sgy``, ``.segy``): for seismic and wavelets exchanged with
  other software. Read: IBM or IEEE float and integer sample formats, either
  byte order. Written: IEEE float, big-endian, SEG-Y revision 1.

Seismic Unix (``.su``) files from the Java Delivery distribution can still be
read, so existing projects can be imported (see :func:`import_java_project`).

Data kinds
----------
``TraceSet``
    Seismic stacks and wavelets: an array of traces with a common sample
    interval and a start time per trace. In HDF5/NPZ: ``data`` (n_traces,
    n_samples), ``dt`` (seconds), ``t0`` (seconds, per trace) and optional
    per-trace header arrays under ``headers/`` (``headers_`` prefix in NPZ).
    A wavelet is a one-trace set whose ``t0`` is negative (the precursor).

Prior tables
    Per-trace prior values that vary areally, as named arrays: XML tag name
    -> array of shape (n_traces,) or (n_traces, n_layers) or, for ``time``,
    (n_traces, n_layers + 1) with the base time last. Contact depths are
    (n_traces, n_contacts). NaN or -999 means "use the XML value". Times are
    in milliseconds, as in the XML.

Results
    Posterior samples: one array per property with shape (n_traces,
    n_samples, n_layers), written by :func:`write_results`.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

BAD_VALUE = -999.0

# ----------------------------------------------------------------------------
# trace header layout shared by SU and SEG-Y (bytes 1-180 are standard; the
# SU and SEG-Y rev 1 definitions differ after that)

_COMMON_SPEC = (
    "tracl,i,0 tracr,i,4 fldr,i,8 tracf,i,12 ep,i,16 cdp,i,20 cdpt,i,24 "
    "trid,h,28 nvs,h,30 nhs,h,32 duse,h,34 offset,i,36 gelev,i,40 selev,i,44 "
    "sdepth,i,48 gdel,i,52 sdel,i,56 swdep,i,60 gwdep,i,64 scalel,h,68 "
    "scalco,h,70 sx,i,72 sy,i,76 gx,i,80 gy,i,84 counit,h,88 wevel,h,90 "
    "swevel,h,92 sut,h,94 gut,h,96 sstat,h,98 gstat,h,100 tstat,h,102 "
    "laga,h,104 lagb,h,106 delrt,h,108 muts,h,110 mute,h,112 ns,H,114 "
    "dt,H,116 gain,h,118 igc,h,120 igi,h,122 corr,h,124 sfs,h,126 sfe,h,128 "
    "slen,h,130 styp,h,132 stas,h,134 stae,h,136 tatyp,h,138 afilf,h,140 "
    "afils,h,142 nofilf,h,144 nofils,h,146 lcf,h,148 hcf,h,150 lcs,h,152 "
    "hcs,h,154 year,h,156 day,h,158 hour,h,160 minute,h,162 sec,h,164 "
    "timbas,h,166 trwf,h,168 grnors,h,170 grnofr,h,172 grnlof,h,174 "
    "gaps,h,176 otrav,h,178"
)
_SU_EXTRA = "d1,f,180 f1,f,184 d2,f,188 f2,f,192 ungpow,f,196 unscale,f,200 ntr,i,204 mark,h,208"
_SEGY_EXTRA = "cdpx,i,180 cdpy,i,184 iline,i,188 xline,i,192 shotpoint,i,196 scalsp,h,200"
_NP = {"i": "i4", "h": "i2", "H": "u2", "f": "f4"}


def _fields(spec: str) -> dict[str, tuple[str, int]]:
    out = {}
    for item in spec.split():
        n, t, o = item.split(",")
        out[n] = (_NP[t], int(o))
    return out


SU_FIELDS = _fields(_COMMON_SPEC + " " + _SU_EXTRA)
SEGY_FIELDS = _fields(_COMMON_SPEC + " " + _SEGY_EXTRA)
# header words worth carrying along with traces
KEEP_HEADERS = ("tracl", "fldr", "ep", "cdp", "offset", "sx", "sy", "gx", "gy",
                "scalco", "cdpx", "cdpy", "iline", "xline")


def _header_dtype(fields: dict, byteorder: str) -> np.dtype:
    names = list(fields)
    return np.dtype({"names": names,
                     "formats": [byteorder + fields[n][0] for n in names],
                     "offsets": [fields[n][1] for n in names],
                     "itemsize": 240})


# ----------------------------------------------------------------------------
# traces


@dataclass
class Trace:
    """A single trace: samples, sample interval (s), start time (s)."""

    data: np.ndarray
    dt: float
    t0: float
    header: dict = field(default_factory=dict)

    @property
    def ns(self) -> int:
        return len(self.data)


@dataclass
class TraceSet:
    data: np.ndarray  # (n_traces, n_samples)
    dt: float  # seconds
    t0: np.ndarray  # (n_traces,) seconds
    headers: dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self):
        self.data = np.atleast_2d(np.asarray(self.data, dtype=float))
        self.t0 = np.broadcast_to(np.asarray(self.t0, dtype=float), (len(self.data),)).copy()

    def __len__(self) -> int:
        return len(self.data)

    def __getitem__(self, i: int) -> Trace:
        h = {k: v[i].item() for k, v in self.headers.items()}
        return Trace(self.data[i], self.dt, float(self.t0[i]), h)

    @property
    def n_samples(self) -> int:
        return self.data.shape[1]


def _kind(path) -> str:
    s = Path(path).suffix.lower()
    if s in (".h5", ".hdf5", ".hdf"):
        return "h5"
    if s == ".npz":
        return "npz"
    if s in (".sgy", ".segy"):
        return "segy"
    if s == ".su":
        return "su"
    raise ValueError(f"unknown file type '{s}' for {path} (use .h5, .npz, .sgy or .su)")


def _h5py():
    try:
        import h5py
    except ImportError as e:  # pragma: no cover
        raise ImportError("HDF5 files need h5py: pip install h5py") from e
    return h5py


def read_traces(path) -> TraceSet:
    """Read seismic or wavelet traces from .h5, .npz, .sgy/.segy or .su."""
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "r") as f:
            hdr = {k: f["headers"][k][()] for k in f["headers"]} if "headers" in f else {}
            t0 = f["t0"][()] if "t0" in f else f.attrs.get("t0", 0.0)
            return TraceSet(f["data"][()], float(f.attrs["dt"]), t0, hdr)
    if kind == "npz":
        z = np.load(path)
        hdr = {k[8:]: z[k] for k in z.files if k.startswith("headers_")}
        return TraceSet(z["data"], float(z["dt"]), z["t0"] if "t0" in z.files else 0.0, hdr)
    if kind == "segy":
        return read_segy(path)
    return read_su(path)


def write_traces(path, ts: TraceSet) -> None:
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "w") as f:
            f.create_dataset("data", data=ts.data.astype(np.float32))
            f.attrs["dt"] = ts.dt
            f.create_dataset("t0", data=ts.t0)
            if ts.headers:
                g = f.create_group("headers")
                for k, v in ts.headers.items():
                    g.create_dataset(k, data=v)
    elif kind == "npz":
        extra = {f"headers_{k}": v for k, v in ts.headers.items()}
        np.savez_compressed(path, data=ts.data.astype(np.float32), dt=ts.dt, t0=ts.t0, **extra)
    elif kind == "segy":
        write_segy(path, ts)
    else:
        raise ValueError("writing .su files is not supported; use .h5, .npz or .sgy")


def read_wavelet(path) -> TraceSet:
    ts = read_traces(path)
    return TraceSet(ts.data[:1], ts.dt, ts.t0[:1], {})


# ---- Seismic Unix (legacy import) ---------------------------------------


def read_su(path, byteorder: str = ">") -> TraceSet:
    """Read a Seismic Unix file (big-endian by default, as in Delivery)."""
    raw = Path(path).read_bytes()
    hdt = _header_dtype(SU_FIELDS, byteorder)
    rows, hdrs = [], []
    pos = 0
    while pos + 240 <= len(raw):
        h = np.frombuffer(raw, dtype=hdt, count=1, offset=pos)[0]
        ns = int(h["ns"])
        rows.append(np.frombuffer(raw, dtype=byteorder + "f4", count=ns, offset=pos + 240))
        hdrs.append(h)
        pos += 240 + 4 * ns
    if not rows:
        raise ValueError(f"no traces in {path}")
    if len({len(r) for r in rows}) != 1:
        raise ValueError(f"{path}: traces have differing lengths")
    h0 = hdrs[0]
    dt = int(h0["dt"]) * 1e-6
    t0 = np.array([h["delrt"] * 1e-3 if h["delrt"] != 0 else float(h["f1"]) for h in hdrs])
    headers = {k: np.array([h[k] for h in hdrs]) for k in KEEP_HEADERS if k in SU_FIELDS}
    return TraceSet(np.stack(rows), dt, t0, headers)


def read_su_blocks(path, byteorder: str = ">") -> np.ndarray:
    """Raw sample matrix of an SU file (used for Java-format prior traces)."""
    return read_su(path, byteorder).data


# ---- SEG-Y ----------------------------------------------------------------


def ibm_to_float(words: np.ndarray) -> np.ndarray:
    """Convert IBM System/360 single-precision floats (as uint32) to float64."""
    w = words.astype(np.uint32)
    sign = np.where(w >> 31, -1.0, 1.0)
    exponent = ((w >> 24) & 0x7F).astype(np.int64) - 64
    mantissa = (w & 0x00FFFFFF).astype(np.float64) / float(1 << 24)
    return sign * mantissa * np.power(16.0, exponent)


_SEGY_FORMATS = {1: ("u4", 4), 2: ("i4", 4), 3: ("i2", 2), 5: ("f4", 4), 8: ("i1", 1)}


def read_segy(path) -> TraceSet:
    raw = Path(path).read_bytes()
    if len(raw) < 3600:
        raise ValueError(f"{path} is too short to be SEG-Y")
    bo = ">"
    marker = np.frombuffer(raw, ">u4", 1, 3296)[0]
    if marker == 0x04030201:
        bo = "<"
    binhdr = lambda off, t: np.frombuffer(raw, bo + t, 1, 3200 + off)[0].item()  # noqa: E731
    dt_us = binhdr(16, "u2")
    ns = binhdr(20, "u2")
    fmt = binhdr(24, "u2")
    if fmt not in _SEGY_FORMATS:
        if bo == ">" and np.frombuffer(raw, "<u2", 1, 3224)[0] in _SEGY_FORMATS:
            bo = "<"
            dt_us, ns, fmt = binhdr(16, "u2"), binhdr(20, "u2"), binhdr(24, "u2")
        else:
            raise ValueError(f"{path}: unsupported SEG-Y sample format {fmt}")
    n_ext = max(binhdr(304, "i2"), 0)
    stype, size = _SEGY_FORMATS[fmt]
    hdt = _header_dtype(SEGY_FIELDS, bo)
    pos = 3600 + 3200 * n_ext
    rows, hdrs = [], []
    while pos + 240 <= len(raw):
        h = np.frombuffer(raw, dtype=hdt, count=1, offset=pos)[0]
        n = int(h["ns"]) or ns
        x = np.frombuffer(raw, bo + stype, n, pos + 240)
        rows.append(ibm_to_float(x) if fmt == 1 else x.astype(np.float64))
        hdrs.append(h)
        pos += 240 + size * n
        if not dt_us:
            dt_us = int(h["dt"])
    if not rows:
        raise ValueError(f"no traces in {path}")
    t0 = np.array([h["delrt"] * 1e-3 for h in hdrs])
    headers = {k: np.array([h[k] for h in hdrs]) for k in KEEP_HEADERS if k in SEGY_FIELDS}
    return TraceSet(np.stack(rows), dt_us * 1e-6, t0, headers)


def write_segy(path, ts: TraceSet, text: str = "Written by pydelivery") -> None:
    """Write IEEE-float big-endian SEG-Y rev 1. Start times are stored in
    ``delrt`` (whole milliseconds), so they must be multiples of 1 ms."""
    ms = ts.t0 * 1e3
    if np.any(np.abs(ms - np.round(ms)) > 1e-6):
        raise ValueError("SEG-Y stores start times in whole ms; use .h5 or .npz instead")
    dt_us = int(round(ts.dt * 1e6))
    th = (text.ljust(3200)[:3200]).encode("ascii")
    bh = np.zeros(400, dtype=np.uint8)
    put = lambda off, v, t: bh.__setitem__(slice(off, off + np.dtype(t).itemsize),  # noqa: E731
                                           np.frombuffer(np.array(v, dtype=">" + t).tobytes(), np.uint8))
    put(16, dt_us, "u2")
    put(20, ts.n_samples, "u2")
    put(24, 5, "u2")  # IEEE float
    put(54, 1, "u2")  # measurement system: metres
    put(300, 0x0100, "u2")  # SEG-Y revision 1
    put(302, 1, "u2")  # fixed-length traces
    hdt = _header_dtype(SEGY_FIELDS, ">")
    with open(path, "wb") as f:
        f.write(th)
        f.write(bh.tobytes())
        for i in range(len(ts)):
            h = np.zeros(1, dtype=hdt)
            for k, v in ts.headers.items():
                if k in SEGY_FIELDS:
                    h[k] = v[i]
            h["tracl"] = h["tracl"] or i + 1
            h["ns"] = ts.n_samples
            h["dt"] = dt_us
            h["delrt"] = int(round(ms[i]))
            f.write(h.tobytes())
            f.write(ts.data[i].astype(">f4").tobytes())


# ----------------------------------------------------------------------------
# prior tables


def read_prior_table(path, cfg=None) -> dict[str, np.ndarray]:
    """Read per-trace priors as named arrays (see module docstring).

    Java-format block files (.su or .sgy) are converted using the XML's
    ``<property_indices>`` links, which needs ``cfg``.
    """
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "r") as f:
            return {k: np.asarray(f[k][()], dtype=float) for k in f if not isinstance(f[k], _h5py().Group)}
    if kind == "npz":
        z = np.load(path)
        return {k: np.asarray(z[k], dtype=float) for k in z.files}
    if cfg is None:
        raise ValueError("block-format prior traces need the model config to decode")
    blocks = read_su(path).data if kind == "su" else read_segy(path).data
    return priors_from_blocks(blocks, cfg)


def write_prior_table(path, table: dict[str, np.ndarray]) -> None:
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "w") as f:
            for k, v in table.items():
                f.create_dataset(k, data=v)
    elif kind == "npz":
        np.savez_compressed(path, **table)
    else:
        raise ValueError("prior tables are written as .h5 or .npz")


def priors_from_blocks(blocks: np.ndarray, cfg) -> dict[str, np.ndarray]:
    """Decode Java Delivery prior traces.

    Property with link block ``k`` (1-based) for layer ``l`` is at
    ``blocks[:, l + (k-1)*n_layers]``; the time block holds n_layers+1
    values. Stack-level values sit in the layer-1 slot and contact values in
    slot ``c`` for contact number ``c``.
    """
    nl = cfg.n_layers
    out = {}
    n_contacts = len(cfg.contacts)
    for tag, k in cfg.links.items():
        name = tag.replace("signal_to_noise", "noise_rms")
        start = (k - 1) * nl
        if tag == "time":
            width = nl + 1
        elif tag in ("contact_depth", "sigma_contact_depth"):
            width = n_contacts
        elif tag in STACK_LEVEL_TAGS:
            width = 1
        else:
            width = nl
        v = blocks[:, start:start + width].astype(float)
        v[v == BAD_VALUE] = np.nan
        # as in Java Delivery, only entries flagged varies_areally are used
        used = [j for j in range(width) if (name, j) in cfg.varying]
        if not used:
            continue
        mask = np.zeros(width, dtype=bool)
        mask[used] = True
        v[:, ~mask] = np.nan
        out[name] = v[:, 0] if tag in STACK_LEVEL_TAGS else v
    return out


STACK_LEVEL_TAGS = (
    "near_stack_noise_rms", "far_stack_noise_rms", "near_stack_signal_to_noise",
    "far_stack_signal_to_noise", "A", "sigma_A", "B", "sigma_B",
)


# ----------------------------------------------------------------------------
# results


def write_results(path, results, cfg=None) -> None:
    """Write posterior samples of all traces to .h5 or .npz.

    Each property is stored as (n_traces, n_samples, n_layers) (or without
    the layer axis for per-sample scalars). Times are in seconds. Trace
    header words of the seismic are stored under ``headers``.
    """
    if not results:
        raise ValueError("no results to write")
    n = min(r.n_samples for r in results)
    keys = [k for k in results[0].samples if k not in ("mode",)]
    arrays = {k: np.stack([r.samples[k][:n] for r in results]) for k in keys}
    hdr_keys = sorted(set().union(*[r.header for r in results]) & set(KEEP_HEADERS))
    headers = {k: np.array([r.header.get(k, 0) for r in results]) for k in hdr_keys}
    states = [";".join(s.label for s in r.states) for r in results]
    layer_names = [lay.name for lay in cfg.layers] if cfg is not None else []
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "w") as f:
            g = f.create_group("samples")
            for k, v in arrays.items():
                g.create_dataset(k, data=v, compression="gzip")
            if headers:
                h = f.create_group("headers")
                for k, v in headers.items():
                    h.create_dataset(k, data=v)
            f.create_dataset("fluid_states", data=np.array(states, dtype="S"))
            f.attrs["layer_names"] = np.array(layer_names, dtype="S")
            f.attrs["format"] = "pydelivery-results-1"
    elif kind == "npz":
        np.savez_compressed(path, **{f"samples_{k}": v for k, v in arrays.items()},
                            **{f"headers_{k}": v for k, v in headers.items()},
                            fluid_states=np.array(states), layer_names=np.array(layer_names))
    else:
        raise ValueError("results are written as .h5 or .npz")


def read_results(path) -> dict:
    """Read a results file into ``{'samples': {...}, 'headers': {...},
    'fluid_states': [...], 'layer_names': [...]}``."""
    kind = _kind(path)
    if kind == "h5":
        with _h5py().File(path, "r") as f:
            samples = {k: f["samples"][k][()] for k in f["samples"]}
            headers = {k: f["headers"][k][()] for k in f["headers"]} if "headers" in f else {}
            states = [s.decode() for s in f["fluid_states"][()]]
            names = [s.decode() for s in f.attrs.get("layer_names", [])]
    else:
        z = np.load(path)
        samples = {k[8:]: z[k] for k in z.files if k.startswith("samples_")}
        headers = {k[8:]: z[k] for k in z.files if k.startswith("headers_")}
        states = list(z["fluid_states"])
        names = list(z["layer_names"])
    return {"samples": samples, "headers": headers,
            "fluid_states": [s.split(";") for s in states], "layer_names": names}


# ----------------------------------------------------------------------------
# importing Java Delivery projects


def import_java_project(xml_path, out_dir, prior_traces=None, fmt: str = "h5") -> Path:
    """Convert a Java Delivery project to Python formats.

    Copies the model XML into ``out_dir`` with its seismic and wavelet file
    names pointing at converted files (``fmt`` is 'h5', 'npz' or 'sgy'), and
    converts block-format prior traces into a named prior table
    ``priors.<h5|npz>``. Returns the path of the new XML.
    """
    from xml.etree import ElementTree as ET

    from .config import load_config

    xml_path, out_dir = Path(xml_path), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = load_config(xml_path)
    tree = ET.parse(xml_path)
    root = tree.getroot()
    done: dict[str, str] = {}

    def convert(elem) -> None:
        src = elem.text.strip()
        if src not in done:
            dst = Path(src).stem + "." + fmt
            write_traces(out_dir / dst, read_traces(cfg.path(src)))
            done[src] = dst
        elem.text = done[src]

    for s in root.find("seismic_data").findall("stack"):
        convert(s.find("filename"))
        convert(s.find("wavelet").find("filename"))
    new_xml = out_dir / xml_path.name
    tree.write(new_xml, encoding="UTF-8", xml_declaration=True)
    if prior_traces is not None:
        table = read_prior_table(prior_traces, cfg)
        write_prior_table(out_dir / ("priors." + ("npz" if fmt == "npz" else "h5")), table)
    return new_xml
