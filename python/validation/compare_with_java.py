"""Compare pydelivery posteriors with a Java ``delivery`` realisations file.

Usage:
    python compare_with_java.py MODEL.xml PRIOR_TRACES.su JAVA_REALISATIONS.su LAYER [N]

LAYER is the 1-based layer to report. Both runs must use the same model and
prior traces. Prints, per trace, fluid probabilities, NG, the layer's top
time, its hydrocarbon saturations and the pinch-out probabilities from each
code. Samples differ draw by draw (different random numbers); statistics
should agree within Monte Carlo error.

This is the only place Java output is read; it needs the Java realisations
file only, not Java itself.
"""
import sys

import numpy as np

from pydelivery import Inversion
from pydelivery.io import read_su

# block layout of Java Delivery realisations (format version 22)
JAVA_FIELDS = ("d", "LFIV", "NG", "vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s", "S_h",
               "rho_h", "fluid_type", "vp_eff", "rho_eff", "Z_eff", "R_near", "R_far",
               "seismic_chisquare", "t_fc", "d_fc", "R_near_fc", "R_far_fc", "S_oil",
               "S_lsg", "S_gas", "vs_eff", "Zs_eff", "t")


def main(xml, prior, java_path, layer, n=2000):
    layer = int(layer) - 1
    inv = Inversion.from_files(xml, prior)
    nl, ntr = inv.cfg.n_layers, inv.n_traces
    java = read_su(java_path).data.reshape(ntr, -1, len(JAVA_FIELDS) * nl + 1)
    col = {k: JAVA_FIELDS.index(k) * nl for k in JAVA_FIELDS}
    res = inv.run(int(n), seed=1)
    np.set_printoptions(precision=3, suppress=True)
    print("trace | P(b,o,l,g) python / java | NG | t (ms) mean+-sd | S_oil | S_gas | "
          "P(pinch-out) python / java")
    for i, r in enumerate(res):
        j = java[i]
        # Java packs sub-layer fluids into higher bits in contacts mode
        jf = np.bincount(j[:, col["fluid_type"] + layer].astype(int) & 3, minlength=4) / len(j)
        jt = j[:, col["t"]:col["t"] + nl + 1]
        jp = (np.diff(jt, axis=1) <= 0).mean(axis=0)
        pt = 1e3 * r.samples["t"][:, layer]
        print(f"{i + 1:3d} | {r.fluid_probabilities()[layer]} {jf} | "
              f"{r.mean('NG')[layer]:.3f} {j[:, col['NG'] + layer].mean():.3f} | "
              f"{pt.mean():.1f}+-{pt.std():.1f} {jt[:, layer].mean():.1f}+-{jt[:, layer].std():.1f} | "
              f"{r.mean('S_oil')[layer]:.3f} {j[:, col['S_oil'] + layer].mean():.3f} | "
              f"{r.mean('S_gas')[layer]:.3f} {j[:, col['S_gas'] + layer].mean():.3f} | "
              f"{r.pinchout_probability()} {jp}")


if __name__ == "__main__":
    main(*sys.argv[1:])
