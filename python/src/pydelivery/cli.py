"""Command line: ``pydelivery invert | summary | import-java | convert``."""
from __future__ import annotations

import argparse
import sys
import time

import numpy as np

from .fluids import FLUID_NAMES
from .inversion import Inversion, invert_trace
from .io import import_java_project, read_results, read_traces, write_results, write_traces


def _invert(a: argparse.Namespace) -> int:
    inv = Inversion.from_files(a.model, a.priors)
    traces = range(inv.n_traces)
    if a.traces:
        traces = [int(t) - 1 for t in a.traces.split(",")]
    rng = np.random.default_rng(a.seed)
    results = []
    t0 = time.time()
    for i in traces:
        cfg, seis = inv.trace_setup(i)
        r = invert_trace(cfg, seis, a.N, rng, n_chains=a.chains, burn_in=a.burn_in,
                         thin=a.thin, n_random_starts=a.random_starts,
                         ignore_seismic=a.ignore_seismic, use_isopachs=not a.ignore_isopachs,
                         use_stacking=not a.ignore_stacking, header=inv.stacks[0][i].header)
        results.append(r)
        if a.verbose:
            fs = ", ".join(f"{k}:{v:.2f}" for k, v in r.fluid_state_probabilities().items() if v > 0)
            pin = " ".join(f"{p:.2f}" for p in r.pinchout_probability())
            print(f"trace {i + 1}: modes={len(r.modes)} fluid states [{fs}] "
                  f"P(pinch-out) [{pin}] accept rw={r.rw_acceptance:.2f} "
                  f"jump={r.jump_acceptance:.2f}", file=sys.stderr)
    write_results(a.o, results, inv.cfg)
    if a.verbose:
        print(f"wrote {sum(r.n_samples for r in results)} samples to {a.o} "
              f"in {time.time() - t0:.1f}s", file=sys.stderr)
    return 0


def _summary(a: argparse.Namespace) -> int:
    """Per-trace, per-layer fluid and pinch-out probabilities, and P10/50/90
    of a chosen property."""
    res = read_results(a.results)
    s = res["samples"]
    t = np.concatenate([s["t"], s["t_base"][..., None]], axis=2)
    pinch = (np.diff(t, axis=2) <= 0).mean(axis=1)
    ft = s["fluid_type"].astype(int)
    names = res["layer_names"] or [str(k + 1) for k in range(ft.shape[2])]
    if a.property is not None and a.property not in s:
        raise SystemExit(f"no property '{a.property}' in {a.results}; have {sorted(s)}")
    prop = s[a.property] if a.property else None
    head = "trace  layer             " + "  ".join(f"P({n})" for n in FLUID_NAMES) + "  P(pinch)"
    if prop is not None:
        head += f"  {a.property} P10/P50/P90"
    print(head)
    labels = res["headers"].get("tracl", np.arange(1, ft.shape[0] + 1))
    for k in range(ft.shape[0]):
        for layer in range(ft.shape[2]):
            pf = [(ft[k, :, layer] == j).mean() for j in range(4)]
            line = (f"{int(labels[k]):5d}  {layer + 1:2d} {names[layer][:14]:14s} "
                    + "  ".join(f"{p:8.3f}" for p in pf) + f"  {pinch[k, layer]:8.3f}")
            if prop is not None and prop.ndim == 3:
                q = np.quantile(prop[k, :, layer], [0.1, 0.5, 0.9])
                line += "  " + " / ".join(f"{v:.4g}" for v in q)
            print(line)
    return 0


def _import(a: argparse.Namespace) -> int:
    xml = import_java_project(a.model, a.out_dir, a.priors, a.format)
    print(f"wrote {xml}")
    return 0


def _convert(a: argparse.Namespace) -> int:
    write_traces(a.output, read_traces(a.input))
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="pydelivery", description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("invert", help="sample the posterior at each trace")
    q.add_argument("model", help="ModelDescription.xml")
    q.add_argument("-m", "--priors", help="prior table (.h5/.npz) or Java prior traces (.su/.sgy)")
    q.add_argument("-o", default="results.h5", help="output file (.h5 or .npz)")
    q.add_argument("-N", type=int, default=100, help="samples per trace")
    q.add_argument("--traces", help="comma-separated 1-based trace numbers")
    q.add_argument("--seed", type=int)
    q.add_argument("--chains", type=int, default=64)
    q.add_argument("--burn-in", type=int, default=1000)
    q.add_argument("--thin", type=int, default=40)
    q.add_argument("--random-starts", type=int, default=None,
                   help="extra random starting points per fluid state for mode finding "
                        "(default 0, or 8 in contacts mode)")
    q.add_argument("-IS", "--ignore-seismic", action="store_true", help="sample the prior")
    q.add_argument("-II", "--ignore-isopachs", action="store_true")
    q.add_argument("-IVS", "--ignore-stacking", action="store_true",
                   help="switch off vertical model-stacking couplings")
    q.add_argument("-v", "--verbose", action="store_true")
    q.set_defaults(func=_invert)

    s = sub.add_parser("summary", help="fluid and pinch-out probabilities per trace")
    s.add_argument("results", help="results file from 'invert'")
    s.add_argument("--property", default=None,
                   help="also print P10/P50/P90 of this per-layer property (e.g. NG, phi_s, z)")
    s.set_defaults(func=_summary)

    c = sub.add_parser("import-java", help="convert a Java Delivery project to HDF5/NPZ/SEG-Y")
    c.add_argument("model", help="ModelDescription.xml of the Java project")
    c.add_argument("out_dir")
    c.add_argument("-m", "--priors", help="Java prior traces (.su)")
    c.add_argument("--format", choices=("h5", "npz", "sgy"), default="h5")
    c.set_defaults(func=_import)

    v = sub.add_parser("convert", help="convert a trace file between .su/.sgy/.h5/.npz")
    v.add_argument("input")
    v.add_argument("output")
    v.set_defaults(func=_convert)

    a = p.parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
