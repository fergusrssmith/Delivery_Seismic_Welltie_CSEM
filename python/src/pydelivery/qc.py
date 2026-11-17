"""Standard inversion QC sections and plot.

The usual way to check an inversion of synthetic data is five sections side
by side: the true impedance, the input seismic, the estimated impedance, the
synthetic seismic of the estimate, and the residual (input minus synthetic).
A good inversion has an estimate that looks like the truth and a residual
that looks like the noise.
"""
from __future__ import annotations

import numpy as np


def posterior_sections(results, times, stack: int = 0, seed: int = 0) -> dict[str, np.ndarray]:
    """Stack per-trace posterior traces into sections of shape (n_traces, n_times).

    Keys: ``Zp_mean``, ``Zs_mean``, ``synthetic_mean`` (posterior means) and
    ``Zp_sample``, ``synthetic_sample`` (one posterior realisation per trace).
    """
    rng = np.random.default_rng(seed)
    rows = [r.time_traces(times, stack=stack, keep=1, rng=rng) for r in results]
    out = {k: np.stack([row[k] for row in rows]) for k in rows[0] if k.endswith("_mean")}
    for k in ("Zp", "synthetic"):
        out[f"{k}_sample"] = np.stack([row[f"{k}_samples"][0] for row in rows])
    return out


def plot_qc(times, seismic, Z_true, Z_est, synthetic, *, x=None, Z_sample=None,
            time_range=None, title=None, Z_label="P impedance", seismic_clip=None):
    """Plot the standard five-panel inversion QC (plus an optional single
    realisation). All arrays are (n_traces, n_times); ``times`` in seconds.

    Impedance panels share one colour scale; seismic, synthetic and residual
    share another, so the residual is directly comparable with the data.
    """
    import matplotlib.pyplot as plt

    times = np.asarray(times)
    n = seismic.shape[0]
    x = np.arange(1, n + 1) if x is None else np.asarray(x)
    t_ms = 1000 * times
    half = 0.5 * (t_ms[1] - t_ms[0]) if len(t_ms) > 1 else 0.5
    extent = [x[0] - 0.5, x[-1] + 0.5, t_ms[-1] + half, t_ms[0] - half]  # pixels centred on samples
    zlo, zhi = np.percentile(np.concatenate([Z_true.ravel(), Z_est.ravel()]), [1, 99])
    clip = seismic_clip or np.percentile(np.abs(seismic), 99)
    residual = seismic - synthetic
    panels = [
        ("1. Input (true) impedance", Z_true, "viridis", (zlo, zhi)),
        ("2. Input seismic", seismic, "seismic_r", (-clip, clip)),
        ("3. Estimated impedance (posterior mean)", Z_est, "viridis", (zlo, zhi)),
        ("4. Synthetic from the estimate", synthetic, "seismic_r", (-clip, clip)),
        ("5. Residual = input - synthetic", residual, "seismic_r", (-clip, clip)),
    ]
    if Z_sample is not None:
        panels.append(("One posterior realisation (impedance)", Z_sample, "viridis", (zlo, zhi)))
    fig, axes = plt.subplots(3, 2, figsize=(12, 11), sharex=True, sharey=True)
    order = [axes[0, 0], axes[1, 0], axes[0, 1], axes[1, 1], axes[2, 1], axes[2, 0]]
    for ax, (name, data, cmap, (lo, hi)) in zip(order, panels):
        im = ax.imshow(data.T, aspect="auto", cmap=cmap, vmin=lo, vmax=hi, extent=extent,
                       interpolation="nearest")
        ax.set_title(name, fontsize=10)
        fig.colorbar(im, ax=ax, label=Z_label if cmap == "viridis" else "amplitude", fraction=0.046)
    if Z_sample is None:
        axes[2, 0].axis("off")
    if time_range is not None:
        axes[0, 0].set_ylim(1000 * time_range[1], 1000 * time_range[0])
    for ax in axes[:, 0]:
        ax.set_ylabel("time (ms)")
    for ax in axes[2, :]:
        ax.set_xlabel("trace")
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    rms = lambda a: float(np.sqrt(np.mean(a**2)))  # noqa: E731
    return fig, {"rms_seismic": rms(seismic), "rms_residual": rms(residual)}


# ----------------------------------------------------------------------------
# facies classification


def kmeans(X: np.ndarray, k: int, n_init: int = 10, n_iter: int = 200, seed: int = 0):
    """Plain K-means (k-means++ starts, best of ``n_init``) on rows of ``X``.

    Returns (labels, centres). Features should be on comparable scales.
    """
    rng = np.random.default_rng(seed)
    best = (np.inf, None, None)
    for _ in range(n_init):
        c = [X[rng.integers(len(X))]]
        for _ in range(1, k):
            d2 = np.min([np.sum((X - ci) ** 2, axis=1) for ci in c], axis=0)
            c.append(X[rng.choice(len(X), p=d2 / d2.sum())])
        C = np.array(c)
        for _ in range(n_iter):
            lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), axis=1)
            Cn = np.array([X[lab == j].mean(0) if np.any(lab == j) else C[j] for j in range(k)])
            if np.allclose(Cn, C):
                break
            C = Cn
        inertia = ((X - C[lab]) ** 2).sum()
        if inertia < best[0]:
            best = (inertia, lab, C)
    return best[1], best[2]


def model_facies(cfg) -> list[dict]:
    """The facies the model can produce: one per distinct rock and fluid
    combination, with its prior-mean P and S impedance."""
    from .describe import implied_properties

    tab = implied_properties(cfg)
    out, seen = [], set()
    for (layer, fluid), row in tab.iterrows():
        i = int(layer.split(".")[0]) - 1
        lay = cfg.layers[i]
        name = f"{lay.sand} + {fluid}" if lay.is_permeable else lay.shale
        if name in seen:
            continue
        seen.add(name)
        rho = row["density"]
        out.append({"name": name, "Zp": row["vp (m/s)"] * rho, "Zs": row["vs (m/s)"] * rho,
                    "permeable": lay.is_permeable, "fluid": fluid, "layer": i})
    return out


def facies_present(results, facies: list[dict], layer: int, min_fraction: float = 0.05) -> list[dict]:
    """Facies the inversion says are present: impermeable facies always, and a
    sand-fluid facies if it is the single most probable outcome of ``layer``
    at no fewer than ``min_fraction`` of the traces."""
    from .fluids import FLUID_NAMES

    probs = np.array([r.fluid_probabilities()[layer] for r in results])  # absent, b, o, l, g
    best = np.argmax(probs, axis=1)
    keep = []
    for f in facies:
        if not f["permeable"]:
            keep.append(f)
            continue
        col = 1 + FLUID_NAMES.index(f["fluid"])
        if np.mean(best == col) >= min_fraction:
            keep.append(f)
    return keep


def classify_facies(Zp: np.ndarray, Zs: np.ndarray, facies: list[dict], seed: int = 0):
    """K-means on (Zp, Zs) samples with one cluster per facies.

    Each cluster is named after its nearest facies (prior-mean impedances,
    in standardised units). Names are not forced to be unique: if two
    clusters are both nearest to shale, K-means has spent a cluster splitting
    the shale, which is itself a diagnostic (e.g. of trace-to-trace jitter).

    Returns (labels shaped like Zp, as indices into ``facies``; cluster
    centres in (Zp, Zs) units; facies index of each cluster).
    """
    X = np.column_stack([Zp.ravel(), Zs.ravel()])
    ref = np.array([[f["Zp"], f["Zs"]] for f in facies])
    scale = X.std(axis=0)
    lab, C = kmeans(X / scale, len(facies), seed=seed)
    C = C * scale
    nearest = np.argmin((((C[:, None, :] - ref[None]) / scale) ** 2).sum(-1), axis=1)
    return nearest[lab].reshape(Zp.shape), C, nearest
