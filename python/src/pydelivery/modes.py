"""Locate the local modes of the posterior for each fluid state.

For a fixed fluid state the posterior over the continuous parameters can be
multi-modal: a layer can be present or pinched out, and the seismic misfit
oscillates in the layer times (loop skipping).  Each mode is found by a
damped Gauss-Newton (Levenberg-Marquardt) fit from a deliberately chosen
start, and summarised by its peak, a local covariance and a Laplace estimate
of its probability mass.  The modes of all fluid states together are the
discrete "models" the sampler jumps between.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .model import StateModel

FD_REL_STEP = 1e-4  # finite-difference step, as a fraction of the prior sigma


@dataclass
class Mode:
    state_index: int
    model: StateModel
    mean: np.ndarray
    cov: np.ndarray
    chol: np.ndarray
    nll: float  # 0.5*|residuals|^2 at the peak
    seismic_nll: float
    log_mass: float  # Laplace log marginal likelihood (no fluid-state prior)
    present: tuple[bool, ...]  # which layers have non-zero thickness
    log_prior_prob: float = 0.0  # log of this mode's share of the state prior
    weight: float = 0.0  # estimated posterior probability of the mode

    @property
    def label(self) -> str:
        layers = "".join("1" if x else "0" for x in self.present)
        return f"{self.model.state.label}/{layers}"


def jacobian(model: StateModel, m: np.ndarray, ignore_seismic: bool = False) -> np.ndarray:
    """Central finite-difference Jacobian of the residuals, one batched call."""
    h = FD_REL_STEP * np.sqrt(np.diag(model.prior_cov))
    d = model.dim
    pts = np.concatenate([m + np.diag(h), m - np.diag(h)])
    r = model.evaluate(pts, ignore_seismic).residuals
    return ((r[:d] - r[d:]) / (2.0 * h[:, None])).T


def _fit(model: StateModel, start: np.ndarray, ignore_seismic: bool,
         max_iter: int = 80) -> tuple[np.ndarray, float] | None:
    """Levenberg-Marquardt minimisation of 0.5*|r(m)|^2 from ``start``."""
    m = start.copy()
    ev = model.evaluate(m, ignore_seismic)
    if not ev.valid[0]:
        return None
    f = float(ev.nll[0])
    r = ev.residuals[0]
    lam = 1e-3
    for _ in range(max_iter):
        J = jacobian(model, m, ignore_seismic)
        H = J.T @ J
        g = J.T @ r
        improved = False
        for _ in range(12):
            try:
                step = -np.linalg.solve(H + lam * np.diag(np.diag(H)), g)
            except np.linalg.LinAlgError:
                lam *= 10.0
                continue
            ev2 = model.evaluate(m + step, ignore_seismic)
            f2 = float(ev2.nll[0])
            if f2 < f:
                improved = True
                break
            lam *= 4.0
        if not improved:
            break
        gain = f - f2
        m, f, r = m + step, f2, ev2.residuals[0]
        lam = max(lam / 3.0, 1e-10)
        if gain < 1e-9 * (1.0 + f):
            break
    return m, f


def time_scan_starts(model: StateModel, rng: np.random.Generator, n_scan: int = 2000,
                     n_keep: int = 3, min_separation: float = 0.002) -> list[np.ndarray]:
    """Global search over the layer times (paper section 4.2).

    The seismic misfit is strongly multi-modal in the layer times (loop
    skipping), so a local optimiser started from the picked times can lock
    onto a wrong cycle. Draw ``n_scan`` sets of layer times from their prior,
    with every other parameter at its prior mean, score them all in one
    batched call, and return the best ``n_keep`` well-separated ones as
    starting points.
    """
    mu, C = model.prior_mean, model.prior_cov
    ti = np.unique(model.time_index[model.time_index >= 0])
    if not len(ti):
        return []
    L = np.linalg.cholesky(C[np.ix_(ti, ti)])
    x = np.repeat(mu[None, :], n_scan, axis=0)
    x[:, ti] = mu[ti] + rng.standard_normal((n_scan, len(ti))) @ L.T
    f = model.evaluate(x).nll
    starts: list[np.ndarray] = []
    for k in np.argsort(f):
        if not np.isfinite(f[k]) or len(starts) == n_keep:
            break
        if all(np.max(np.abs(x[k, ti] - s[ti])) > min_separation for s in starts):
            starts.append(x[k])
    return starts


def starting_points(model: StateModel, rng: np.random.Generator,
                    n_random: int = 0, time_scan: bool = True) -> list[np.ndarray]:
    """Prior mean, one pinched-out start per layer that could plausibly pinch
    out, the best layer-time configurations from a global time scan (when
    there is seismic), and optionally a few random draws from the prior."""
    mu, C = model.prior_mean, model.prior_cov
    starts = [mu.copy()]
    if time_scan and model.seismic is not None:
        starts.extend(time_scan_starts(model, rng))
    ti = model.time_index
    for i in range(model.nl):
        a, b = ti[i], ti[i + 1]
        if a < 0 or b < 0 or a == b:  # fixed, or tied by a forced pinch-out
            continue
        thick = mu[b] - mu[a]
        sd_thick = np.sqrt(max(C[a, a] + C[b, b] - 2.0 * C[a, b], 0.0))
        if thick <= 0 or sd_thick / thick > 0.5:  # roughly a >5% prior chance
            sa, sb = np.sqrt(C[a, a]), np.sqrt(C[b, b])
            s = mu.copy()
            s[a] = (sb * mu[a] + sa * mu[b]) / (sa + sb)
            s[b] = s[a] - 2.0 * FD_REL_STEP * sb
            starts.append(s)
    if n_random:
        L = np.linalg.cholesky(C)
        for _ in range(n_random):
            starts.append(mu + L @ rng.standard_normal(model.dim))
    return starts


def find_modes(model: StateModel, state_index: int, rng: np.random.Generator,
               n_random_starts: int = 0, ignore_seismic: bool = False) -> list[Mode]:
    modes: list[Mode] = []
    for start in starting_points(model, rng, n_random_starts, time_scan=not ignore_seismic):
        fit = _fit(model, start, ignore_seismic)
        if fit is None:
            continue
        m, f = fit
        J = jacobian(model, m, ignore_seismic)
        H = J.T @ J
        try:
            cov = np.linalg.inv(H)
            cov = 0.5 * (cov + cov.T)
            chol = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            continue
        # same peak as a mode we already have?
        if any(_same(m, cov, other) for other in modes):
            continue
        ev = model.evaluate(m, ignore_seismic)
        p = model.physical(m)
        log_mass = -f + 0.5 * np.linalg.slogdet(cov)[1] - 0.5 * model.logdet_prior_cov
        modes.append(Mode(state_index, model, m, cov, chol, f, float(ev.seismic_nll[0]),
                          float(log_mass), tuple(bool(x) for x in (p["t"][0, 1:] > p["t"][0, :-1]))))
    return modes


def _same(m: np.ndarray, cov: np.ndarray, other: Mode) -> bool:
    d = m - other.mean
    sd = np.sqrt(np.minimum(np.diag(cov), np.diag(other.cov)))
    return bool(np.max(np.abs(d) / sd) < 0.5)


def weigh_modes(modes: list[Mode], state_probs: list[float], min_weight: float = 1e-10) -> list[Mode]:
    """Assign mode probabilities and drop negligible modes.

    A fluid state's prior probability P_k is shared among its modes in
    proportion to their Laplace masses, so the number of modes found for a
    state does not change that state's prior weight.  The estimated posterior
    probability of mode j is then  w_j ~ mass_j * share_j.
    """
    if not modes:
        return []
    for _ in range(2):  # second pass re-shares priors after pruning
        logw = np.empty(len(modes))
        for k in set(mo.state_index for mo in modes):
            idx = [i for i, mo in enumerate(modes) if mo.state_index == k]
            lm = np.array([modes[i].log_mass for i in idx])
            share = lm - (lm.max() + np.log(np.sum(np.exp(lm - lm.max()))))
            for i, s in zip(idx, share):
                with np.errstate(divide="ignore"):
                    modes[i].log_prior_prob = float(np.log(state_probs[k]) + s)
                logw[i] = modes[i].log_mass + modes[i].log_prior_prob
        w = np.exp(logw - logw.max())
        w /= w.sum()
        for mo, wi in zip(modes, w):
            mo.weight = float(wi)
        kept = [mo for mo in modes if mo.weight >= min_weight]
        if len(kept) == len(modes):
            break
        modes = kept
    return modes
