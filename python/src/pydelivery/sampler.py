"""Mode-jumping Markov chain Monte Carlo over fluid states and layer models.

The chain state is ``(mode, m)``: a mode label (which fixes the fluid state)
and that state's model vector.  Two kinds of move are mixed 50/50:

* **within-mode**: a scaled random-walk Metropolis step using the mode's
  local covariance;
* **mode jump**: propose another mode, translate the parameters the two
  fluid states share by the offset between the mode peaks, and draw any
  parameters that exist only in the new state (e.g. gas properties) from
  their prior.  This is a reversible jump between models of different
  dimension.

Proposals are scored on the *physical* model (clipped, time-ordered), so
pinched-out layers and other boundary states are sampled with their proper
probability.  Many chains are advanced together; chains in the same mode are
evaluated in one batched forward-model call.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .modes import Mode

LOG_2PI = np.log(2.0 * np.pi)


@dataclass
class ChainSamples:
    mode: np.ndarray  # (n_samples,) mode index of each sample
    m: list[np.ndarray]  # model vector of each sample (length depends on mode)
    rw_acceptance: float
    jump_acceptance: float


class _Jump:
    """Precomputed mapping for a jump from mode a to mode b."""

    def __init__(self, a: Mode, b: Mode):
        ka, kb = a.model.index, b.model.index
        common = [k for k in b.model.keys if k in ka]
        self.src = np.array([ka[k] for k in common], dtype=int)
        self.dst = np.array([kb[k] for k in common], dtype=int)
        self.shift = b.mean[self.dst] - a.mean[self.src]
        self.new = np.array([kb[k] for k in b.model.keys if k not in ka], dtype=int)
        self.dropped = np.array([ka[k] for k in a.model.keys if k not in kb], dtype=int)
        self.new_mu = b.model.prior_mean[self.new]
        self.new_sd = np.sqrt(np.diag(b.model.prior_cov)[self.new])
        self.drop_mu = a.model.prior_mean[self.dropped]
        self.drop_sd = np.sqrt(np.diag(a.model.prior_cov)[self.dropped])


def _log_normal(x, mu, sd):
    if x.shape[-1] == 0:
        return np.zeros(x.shape[0])
    return np.sum(-0.5 * ((x - mu) / sd) ** 2 - np.log(sd) - 0.5 * LOG_2PI, axis=-1)


def run_chains(modes: list[Mode], n_samples: int, rng: np.random.Generator,
               n_chains: int = 64, burn_in: int = 1000, thin: int = 40,
               ignore_seismic: bool = False) -> ChainSamples:
    """Draw ``n_samples`` posterior samples (rounded up to a multiple of
    ``n_chains``)."""
    n_modes = len(modes)
    w = np.array([mo.weight for mo in modes])
    w = np.clip(w / w.sum(), 1e-12, 1.0 - 1e-12)
    # log target = log(mode prior share) + log prior normalisation - nll
    const = np.array([
        mo.log_prior_prob - 0.5 * mo.model.dim * LOG_2PI - 0.5 * mo.model.logdet_prior_cov
        for mo in modes
    ])
    jumps: dict[tuple[int, int], _Jump] = {}

    def nll(j: int, x: np.ndarray) -> np.ndarray:
        return modes[j].model.evaluate(x, ignore_seismic).nll

    # start each chain from the Laplace approximation of a mode
    mode = rng.choice(n_modes, size=n_chains, p=w / w.sum())
    m: list[np.ndarray] = [None] * n_chains
    f = np.empty(n_chains)
    for j in range(n_modes):
        idx = np.flatnonzero(mode == j)
        if not len(idx):
            continue
        mo = modes[j]
        x = mo.mean + rng.standard_normal((len(idx), mo.model.dim)) @ mo.chol.T
        fx = nll(j, x)
        bad = ~np.isfinite(fx)
        x[bad] = mo.mean
        fx[bad] = mo.nll
        for c, xi, fi in zip(idx, x, fx):
            m[c], f[c] = xi, fi

    n_keep = -(-n_samples // n_chains)
    out_mode: list[np.ndarray] = []
    out_m: list[np.ndarray] = []
    rw_try = rw_acc = jump_try = jump_acc = 0

    for it in range(burn_in + n_keep * thin):
        do_jump = (rng.random(n_chains) < 0.5) if n_modes > 1 else np.zeros(n_chains, bool)
        # ---- within-mode random walk ---------------------------------
        for j in np.unique(mode[~do_jump]):
            idx = np.flatnonzero((mode == j) & ~do_jump)
            mo = modes[j]
            d = mo.model.dim
            x = np.stack([m[c] for c in idx])
            prop = x + (2.4 / np.sqrt(d)) * rng.standard_normal((len(idx), d)) @ mo.chol.T
            fp = nll(j, prop)
            accept = np.log(rng.random(len(idx))) < f[idx] - fp
            rw_try += len(idx)
            rw_acc += int(accept.sum())
            for c, a, xi, fi in zip(idx, accept, prop, fp):
                if a:
                    m[c], f[c] = xi, fi
        # ---- jumps between modes -------------------------------------
        if do_jump.any():
            target = mode.copy()
            for c in np.flatnonzero(do_jump):
                q = w.copy()
                q[mode[c]] = 0.0
                target[c] = rng.choice(n_modes, p=q / q.sum())  # Q(j'|j) = w_j'/(1-w_j)
            pairs = {(int(a), int(b)) for a, b in zip(mode[do_jump], target[do_jump])}
            for a, b in pairs:
                idx = np.flatnonzero(do_jump & (mode == a) & (target == b))
                if (a, b) not in jumps:
                    jumps[(a, b)] = _Jump(modes[a], modes[b])
                jp = jumps[(a, b)]
                x = np.stack([m[c] for c in idx])
                prop = np.empty((len(idx), modes[b].model.dim))
                prop[:, jp.dst] = x[:, jp.src] + jp.shift
                new = jp.new_mu + jp.new_sd * rng.standard_normal((len(idx), len(jp.new)))
                prop[:, jp.new] = new
                fp = nll(b, prop)
                log_alpha = (
                    (const[b] - fp) - (const[a] - f[idx])
                    + _log_normal(x[:, jp.dropped], jp.drop_mu, jp.drop_sd)
                    - _log_normal(new, jp.new_mu, jp.new_sd)
                    + np.log(w[a] / (1.0 - w[b])) - np.log(w[b] / (1.0 - w[a]))
                )
                accept = np.log(rng.random(len(idx))) < log_alpha
                jump_try += len(idx)
                jump_acc += int(accept.sum())
                for c, acc, xi, fi in zip(idx, accept, prop, fp):
                    if acc:
                        m[c], f[c], mode[c] = xi, fi, b
        if it >= burn_in and (it - burn_in) % thin == thin - 1:
            out_mode.append(mode.copy())
            out_m.extend(x.copy() for x in m)

    return ChainSamples(
        np.concatenate(out_mode)[:n_samples],
        out_m[:n_samples],
        rw_acc / max(rw_try, 1),
        jump_acc / max(jump_try, 1),
    )
