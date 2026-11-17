"""Command-line entry point: ``pydelivery invert`` and ``pydelivery summary``."""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from .fluids import FLUID_NAMES
from .inversion import Inversion, write_realisations


def _invert(a: argparse.Namespace) -> int:
    inv = Inversion.from_files(a.model, a.m)
    traces = range(inv.n_traces)
    if a.traces:
        traces = [int(t) - 1 for t in a.traces.split(",")]
    rng = np.random.default_rng(a.seed)
    results = []
    t0 = time.time()
    for i in traces:
        cfg, seis = inv.trace_setup(i)
        from .inversion import invert_trace

        r = invert_trace(cfg, seis, a.N, rng, n_chains=a.chains, burn_in=a.burn_in,
                         thin=a.thin, n_random_starts=a.random_starts,
                         ignore_seismic=a.ignore_seismic, use_isopachs=not a.ignore_isopachs,
                         header=inv.stacks[0][i].header)
        results.append(r)
        if a.verbose:
            fs = ", ".join(f"{k}:{v:.2f}" for k, v in r.fluid_state_probabilities().items() if v > 0)
            pin = " ".join(f"{p:.2f}" for p in r.pinchout_probability())
            print(f"trace {i + 1}: modes={len(r.modes)} fluid states [{fs}] "
                  f"P(pinch-out) [{pin}] accept rw={r.rw_acceptance:.2f} "
                  f"jump={r.jump_acceptance:.2f}", file=sys.stderr)
    write_realisations(a.o, results)
    if a.verbose:
        print(f"wrote {sum(r.n_samples for r in results)} realisations to {a.o} "
              f"in {time.time() - t0:.1f}s", file=sys.stderr)
    return 0


def _summary(a: argparse.Namespace) -> int:
    """Per-trace fluid and pinch-out probabilities from a realisations file."""
    from .inversion import SU_FIELDS
    from .su import read_su

    tr = read_su(a.realisations)
    nb = len(SU_FIELDS)
    nl = (tr[0].ns - 1) // nb
    data = np.stack([t.data for t in tr])
    # a new trace location starts whenever the realisation number restarts at 1
    marks = np.array([t.header["mark"] for t in tr])
    starts = list(np.flatnonzero(marks == 1)) + [len(tr)]
    ft_col = SU_FIELDS.index("fluid_type") * nl
    t_col = SU_FIELDS.index("t") * nl
    print("trace  layer  " + "  ".join(f"P({n})" for n in FLUID_NAMES) + "  P(pinch-out)")
    for k in range(len(starts) - 1):
        d = data[starts[k]:starts[k + 1]]
        t = d[:, t_col:t_col + nl + 1]
        pinch = (np.diff(t, axis=1) <= 0).mean(axis=0)
        for layer in range(nl):
            ft = d[:, ft_col + layer].astype(int)
            pf = [(ft == j).mean() for j in range(4)]
            print(f"{k + 1:5d}  {layer + 1:5d}  " + "  ".join(f"{p:8.3f}" for p in pf)
                  + f"  {pinch[layer]:12.3f}")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="pydelivery", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)
    q = sub.add_parser("invert", help="sample the posterior at each trace")
    q.add_argument("model", help="ModelDescription.xml")
    q.add_argument("-m", help="SU file of areally varying prior traces")
    q.add_argument("-o", default="realisations.su", help="output SU file")
    q.add_argument("-N", type=int, default=100, help="realisations per trace")
    q.add_argument("--traces", help="comma-separated 1-based trace numbers")
    q.add_argument("--seed", type=int)
    q.add_argument("--chains", type=int, default=64)
    q.add_argument("--burn-in", type=int, default=1000)
    q.add_argument("--thin", type=int, default=40)
    q.add_argument("--random-starts", type=int, default=0,
                   help="extra random starting points per fluid state for mode finding")
    q.add_argument("-IS", "--ignore-seismic", action="store_true", help="sample the prior")
    q.add_argument("-II", "--ignore-isopachs", action="store_true")
    q.add_argument("-v", "--verbose", action="store_true")
    q.set_defaults(func=_invert)
    s = sub.add_parser("summary", help="fluid and pinch-out probabilities per trace")
    s.add_argument("realisations")
    s.set_defaults(func=_summary)
    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
