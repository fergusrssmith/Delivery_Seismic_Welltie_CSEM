"""Compare pydelivery posteriors with a Java ``delivery`` realisations file.

Usage:
    python compare_with_java.py MODEL.xml PRIOR_TRACES.su JAVA_REALISATIONS.su LAYER [N]

LAYER is the 1-based layer to report.  Both runs must use the same model and
prior traces.  Prints, per trace, fluid probabilities, NG, the layer's top
time and the pinch-out probabilities from each code.  Samples differ draw by
draw (different random numbers); the statistics should agree within Monte
Carlo error.
"""
import sys

import numpy as np

from pydelivery import Inversion
from pydelivery.inversion import SU_FIELDS
from pydelivery.su import read_su


def main(xml, prior, java_path, layer, n=2000):
    layer = int(layer) - 1
    inv = Inversion.from_files(xml, prior)
    nl, ntr = inv.cfg.n_layers, inv.n_traces
    java = np.stack([t.data for t in read_su(java_path)])
    java = java.reshape(ntr, -1, len(SU_FIELDS) * nl + 1)
    col = {k: SU_FIELDS.index(k) * nl for k in SU_FIELDS}
    res = inv.run(int(n), seed=1)
    np.set_printoptions(precision=3, suppress=True)
    print("trace | P(b,o,l,g) python / java | NG python / java | "
          "t (ms) mean+-sd python / java | P(pinch-out) python / java")
    for i, r in enumerate(res):
        j = java[i]
        jf = np.bincount(j[:, col["fluid_type"] + layer].astype(int), minlength=4) / len(j)
        jt = j[:, col["t"]:col["t"] + nl + 1]
        jp = (np.diff(jt, axis=1) <= 0).mean(axis=0)
        pt = 1e3 * r.samples["t"][:, layer]
        print(f"{i + 1:3d} | {r.fluid_probabilities()[layer]} {jf} | "
              f"{r.mean('NG')[layer]:.3f} {j[:, col['NG'] + layer].mean():.3f} | "
              f"{pt.mean():.1f}+-{pt.std():.1f} {jt[:, layer].mean():.1f}+-{jt[:, layer].std():.1f} | "
              f"{r.pinchout_probability()} {jp}")


if __name__ == "__main__":
    main(*sys.argv[1:])
