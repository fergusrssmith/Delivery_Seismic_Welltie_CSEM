"""Experimental lateral continuity between traces.

Delivery inverts each trace on its own. This module couples a line of traces
through a spatially correlated prior on chosen rock properties: the Bayesian
counterpart of the spatial (Laplacian) regularisation used in deterministic
inversions.

The prior
---------
At every trace the prior of each property stays exactly Delivery's: the
trend-curve mean and scatter (e.g. shale vp 2812 ± 100 m/s). What is added is
correlation between traces: the deviation e_i of trace i from its trend mean
follows a stationary first-order Markov chain along the line,

    e_(i+1) = rho * e_i + w,   w ~ N(0, (1 - rho^2) * Sigma),   rho = exp(-1 / L)

so neighbouring traces are similar, traces far apart are independent, and
every trace's own prior is still N(trend mean, Sigma). ``L`` is the
correlation length in traces, like the range of a variogram. ``L -> 0`` gives
back trace-by-trace Delivery. Fluid states and pinch-outs stay per trace.

How it is computed (one extra inversion pass)
---------------------------------------------
1. Invert every trace on its own (ordinary Delivery).
2. Summarise what each trace's *data* say about the coupled properties:
   its posterior divided by its prior (both as Gaussians).
3. A Kalman smoother along the line combines these with the correlated
   prior. For each trace it gives the information arriving from all the
   *other* traces' data, through the correlation.
4. Re-invert every trace with that information as one extra Gaussian prior
   term.

Step 3 is exact for the Gaussian summaries, so no information is counted
twice, and the line-average level of a property is not claimed to be better
known than the trends allow. The approximation is the Gaussian summary in
step 2, which is good for rock properties but not for discrete outcomes,
which is why fluids and pinch-outs are not coupled.

Limitations: lines only (a 3-D grid needs loopy belief propagation), and the
trend scatter Sigma is taken as the same along the line.
"""
from __future__ import annotations

import copy
import os
from dataclasses import dataclass

import numpy as np

from .inversion import Inversion, TraceResult, _mp_context, invert_trace

DEFAULT_VARIABLES = ("vp_m", "vs_m", "rho_m", "vp_s", "vs_s", "phi_s", "NG")


@dataclass
class LateralCoupling:
    """Which variables to correlate along the line, and over what distance.

    ``variables``: model variable names (vp_m, vs_m, rho_m, vp_s, vs_s, phi_s,
    NG, NG_hard, vp_hard, ..., d, LFIV), coupled in every layer where they are
    free at every trace. ``correlation_length``: in traces.
    """

    variables: tuple = DEFAULT_VARIABLES
    correlation_length: float = 20.0

    @property
    def rho(self) -> float:
        return float(np.exp(-1.0 / self.correlation_length)) if self.correlation_length > 0 else 0.0


@dataclass
class LateralRun:
    results: list[TraceResult]
    keys: list  # coupled variables, in order
    local_mean: np.ndarray  # (n_traces, n_keys) trace-by-trace posterior means
    local_sd: np.ndarray
    smoothed_mean: np.ndarray  # posterior means with lateral coupling (Gaussian smoother)
    smoothed_sd: np.ndarray
    prior_sd: np.ndarray  # (n_keys,) per-trace prior standard deviation
    difference_sd: np.ndarray  # (n_keys,) implied prior sd of neighbour differences


def _coupled_keys(results, coupling) -> list:
    common = None
    for r in results:
        for mo in r.modes:
            ks = {k for k in mo.model.index if k[0] in coupling.variables}
            common = ks if common is None else common & ks
    return sorted(common or [], key=str)


def _summaries(result: TraceResult, keys: list):
    """Posterior mean/cov of the coupled variables, and their prior mean/cov."""
    X = np.empty((result.n_samples, len(keys)))
    mode_idx = result.samples["mode"]
    for j, mo in enumerate(result.modes):
        idx = np.flatnonzero(mode_idx == j)
        if len(idx):
            cols = [mo.model.index[k] for k in keys]
            X[idx] = np.stack([result.vectors[i] for i in idx])[:, cols]
    model = result.modes[0].model
    cols = [model.index[k] for k in keys]
    return (X.mean(axis=0), np.atleast_2d(np.cov(X, rowvar=False)),
            model.prior_mean[cols], model.prior_cov[np.ix_(cols, cols)])


def _psd(A, floor=0.0):
    w, V = np.linalg.eigh(0.5 * (A + A.T))
    return (V * np.maximum(w, floor)) @ V.T


def correlated_chain(lik_J, lik_h, Sigma, rho):
    """Kalman smoother (Gaussian belief propagation) along a line.

    Deviations e_i from the trend means follow a stationary AR(1) chain with
    marginal covariance ``Sigma`` and correlation ``rho``. Trace i's data give
    the information-form likelihood (lik_J[i], lik_h[i]) on e_i.

    Returns (from_others_J, from_others_h, post_mean, post_cov): the
    information reaching each trace from the other traces' data, *beyond* its
    own prior, and the smoothed posterior of each e_i.
    """
    n, k = lik_h.shape
    Si = np.linalg.inv(Sigma)
    if rho <= 0:
        Z = np.zeros((n, k, k))
        post_cov = np.linalg.inv(lik_J + Si)
        return Z, np.zeros((n, k)), np.einsum("nij,nj->ni", post_cov, lik_h), post_cov
    Qi = Si / (1.0 - rho**2)  # transition noise precision
    Jf = np.zeros((n, k, k)); hf = np.zeros((n, k))
    Jf[0] = Si  # stationary prior enters at the first trace
    for i in range(n - 1):
        J, h = lik_J[i] + Jf[i], lik_h[i] + hf[i]
        A = np.linalg.inv(J + rho**2 * Qi)
        Jf[i + 1] = Qi - rho**2 * Qi @ A @ Qi
        hf[i + 1] = rho * Qi @ A @ h
    Jb = np.zeros((n, k, k)); hb = np.zeros((n, k))
    for i in range(n - 1, 0, -1):
        J, h = lik_J[i] + Jb[i], lik_h[i] + hb[i]
        A = np.linalg.inv(J + Qi)
        Jb[i - 1] = rho**2 * Qi - rho**2 * Qi @ A @ Qi
        hb[i - 1] = rho * Qi @ A @ h
    post_cov = np.linalg.inv(Jf + Jb + lik_J)
    post_mean = np.einsum("nij,nj->ni", post_cov, hf + hb + lik_h)
    # each trace's re-inversion already has its own prior Sigma: remove it
    return Jf + Jb - Si, hf + hb, post_mean, post_cov


def _run_one(job) -> TraceResult:
    inv, i, lateral, n_samples, seed, kwargs = job
    cfg, seis = inv.trace_setup(i)
    cfg = copy.copy(cfg)
    cfg.lateral = lateral
    return invert_trace(cfg, seis, n_samples, np.random.default_rng(seed),
                        header=inv.stacks[0][i].header, **kwargs)


def invert_lateral(inv: Inversion, coupling: LateralCoupling, n_samples: int = 500,
                   seed: int | None = None, n_jobs: int = 1,
                   start: list[TraceResult] | None = None, **kwargs) -> LateralRun:
    """Invert a line of traces with lateral correlation (see module docstring).

    ``start`` is a trace-by-trace result to begin from (computed if not
    given). Traces are taken to be in line order.
    """
    if n_jobs == -1:
        n_jobs = os.cpu_count() or 1
    ss = np.random.SeedSequence(seed)
    if start is None:
        start = inv.run(n_samples, seed=int(ss.generate_state(1)[0]), n_jobs=n_jobs, **kwargs)
    keys = _coupled_keys(start, coupling)
    if not keys:
        raise ValueError("none of the coupled variables is free at every trace")
    n, k = len(start), len(keys)
    post_m, post_C, pri_m, pri_C = map(np.array, zip(*[_summaries(r, keys) for r in start]))
    Sigma = pri_C.mean(axis=0)
    Si = np.linalg.inv(Sigma)
    # what each trace's data say: posterior / prior, as information (on e = x - trend mean)
    lik_J = np.array([_psd(np.linalg.inv(C) - Si) for C in post_C])
    lik_h = np.einsum("nij,nj->ni", np.array([np.linalg.inv(C) for C in post_C]), post_m - pri_m) \
        - np.zeros((n, k))
    oth_J, oth_h, sm_e, sm_C = correlated_chain(lik_J, lik_h, Sigma, coupling.rho)

    jobs = []
    seeds = ss.spawn(n)
    for i in range(n):
        J = _psd(oth_J[i], floor=1e-12 * np.trace(Si) / k)
        cov = np.linalg.inv(J)
        mean = pri_m[i] + cov @ oth_h[i]
        jobs.append((inv, i, {"keys": keys, "mean": mean, "cov": cov}, n_samples, seeds[i], kwargs))
    if n_jobs <= 1:
        results = [_run_one(j) for j in jobs]
    else:
        from concurrent.futures import ProcessPoolExecutor

        with ProcessPoolExecutor(max_workers=n_jobs, mp_context=_mp_context()) as ex:
            results = list(ex.map(_run_one, jobs, chunksize=max(1, len(jobs) // (4 * n_jobs))))
    sd = lambda c: np.sqrt(np.einsum("...ii->...i", c))  # noqa: E731
    prior_sd = sd(Sigma)
    return LateralRun(results, keys, post_m, sd(post_C), pri_m + sm_e, sd(sm_C), prior_sd,
                      prior_sd * np.sqrt(2.0 * (1.0 - coupling.rho)))
