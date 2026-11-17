"""Minimal Seismic Unix (SU) trace reader/writer.

SU files are a sequence of traces, each a 240-byte SEG-Y style trace header
followed by ``ns`` 4-byte IEEE floats.  Delivery works with big-endian
("XDR") SU files, which is the default here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np

HEADER_BYTES = 240

# name -> (type, byte offset).  Same table as the Java sulib.
_FIELDS_SPEC = (
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
    "gaps,h,176 otrav,h,178 d1,f,180 f1,f,184 d2,f,188 f2,f,192 "
    "ungpow,f,196 unscale,f,200 ntr,i,204 mark,h,208 shortpad,h,210"
)
_NP = {"i": "i4", "h": "i2", "H": "u2", "f": "f4"}
FIELDS: dict[str, tuple[str, int]] = {}
for _item in _FIELDS_SPEC.split():
    _n, _t, _o = _item.split(",")
    FIELDS[_n] = (_NP[_t], int(_o))


def header_dtype(byteorder: str = ">") -> np.dtype:
    names = list(FIELDS)
    return np.dtype(
        {
            "names": names,
            "formats": [byteorder + FIELDS[n][0] for n in names],
            "offsets": [FIELDS[n][1] for n in names],
            "itemsize": HEADER_BYTES,
        }
    )


@dataclass
class Trace:
    """One SU trace: a header (dict of the standard fields) and float data."""

    header: dict = field(default_factory=dict)
    data: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.float32))

    @property
    def ns(self) -> int:
        return len(self.data)

    @property
    def dt(self) -> float:
        """Sample interval in seconds."""
        return self.header.get("dt", 0) * 1e-6

    @property
    def t0(self) -> float:
        """Time of first sample in seconds (delrt is in ms)."""
        return self.header.get("delrt", 0) * 1e-3


def read_su(path: str | Path, byteorder: str = ">") -> list[Trace]:
    """Read all traces from an SU file (traces may have differing ns)."""
    raw = Path(path).read_bytes()
    hdt = header_dtype(byteorder)
    traces: list[Trace] = []
    pos = 0
    while pos + HEADER_BYTES <= len(raw):
        h = np.frombuffer(raw, dtype=hdt, count=1, offset=pos)[0]
        ns = int(h["ns"])
        pos += HEADER_BYTES
        data = np.frombuffer(raw, dtype=byteorder + "f4", count=ns, offset=pos)
        pos += 4 * ns
        header = {n: h[n].item() for n in FIELDS}
        traces.append(Trace(header, data.astype(np.float64)))
    return traces


def write_su(path: str | Path, traces: Iterable[Trace], byteorder: str = ">") -> None:
    hdt = header_dtype(byteorder)
    with open(path, "wb") as f:
        for tr in traces:
            h = np.zeros(1, dtype=hdt)
            for k, v in tr.header.items():
                if k in FIELDS:
                    h[k] = v
            h["ns"] = len(tr.data)
            f.write(h.tobytes())
            f.write(np.asarray(tr.data, dtype=byteorder + "f4").tobytes())
