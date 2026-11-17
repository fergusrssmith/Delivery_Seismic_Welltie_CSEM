"""Discrete fluid states: which fluid sits in each layer.

Each permeable layer may hold brine, oil, low-saturation gas (lsg) or gas,
with prior probabilities set per layer. A density-ordering rule then decides
which combinations across layers are geologically allowed (e.g. no oil above
gas in communicating sands). The surviving combinations are the discrete
"fluid states" the inversion chooses between.

In contacts mode (``fluid_units_and_contacts`` in the XML) fluids are set by
depth contacts instead, so there is a single placeholder state; see
:mod:`pydelivery.model`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import DENSITY_PRIORITY, ModelConfig

BRINE, OIL, LSG, GAS = 0, 1, 2, 3
FLUID_NAMES = ("brine", "oil", "lsg", "gas")
FLUID_CHARS = "bolg"


@dataclass(frozen=True)
class FluidState:
    types: tuple[int, ...]  # fluid type per layer (BRINE for impermeable layers)
    probability: float  # normalised prior probability

    @property
    def label(self) -> str:
        return "".join(FLUID_CHARS[t] for t in self.types)


def layer_fluid_probabilities(cfg: ModelConfig) -> np.ndarray:
    """(n_layers, 4) prior probability of each fluid type in each layer."""
    P = np.zeros((cfg.n_layers, 4))
    for i, lay in enumerate(cfg.layers):
        if lay.is_permeable and lay.active:
            for h, t in (("oil", OIL), ("lsg", LSG), ("gas", GAS)):
                hc = lay.hydrocarbons.get(h)
                if hc is not None and hc.probability > 0:
                    P[i, t] = hc.probability
        P[i, BRINE] = max(1.0 - P[i, 1:].sum(), 0.0)
    return P


def enumerate_fluid_states(cfg: ModelConfig) -> list[FluidState]:
    """All fluid combinations allowed by the density-ordering rule.

    ``none``    : any fluid in any permeable layer.
    ``partial`` : ordering applies between directly adjacent permeable layers.
    ``full``    : ordering applies across the whole model; impermeable layers
                  do not act as seals.
    """
    nl = cfg.n_layers
    if cfg.contacts_mode:
        return [FluidState((BRINE,) * nl, 1.0)]
    P = layer_fluid_probabilities(cfg)
    perm = [lay.is_permeable and lay.active for lay in cfg.layers]
    mode = cfg.density_ordering
    states: list[tuple[tuple[int, ...], float]] = []

    def recurse(layer: int, types: list[int], prob: float, lightest_allowed: int) -> None:
        if layer == nl:
            states.append((tuple(types), prob))
            return
        if not perm[layer]:
            nxt = lightest_allowed if mode == "full" else GAS
            recurse(layer + 1, types + [BRINE], prob * P[layer, BRINE], nxt)
            return
        for t in (BRINE, OIL, LSG, GAS):
            if P[layer, t] <= 0:
                continue
            if mode != "none" and DENSITY_PRIORITY[t] > DENSITY_PRIORITY[lightest_allowed]:
                continue  # lighter than the fluid above: not allowed
            recurse(layer + 1, types + [t], prob * P[layer, t], t)

    recurse(0, [], 1.0, GAS)
    if not states:
        raise ValueError("density ordering and fluid probabilities allow no legal fluid states")
    if cfg.post_ordering_probabilities:
        probs = post_ordering_state_probabilities([t for t, _ in states], P)
    else:
        probs = np.array([p for _, p in states])
        probs = probs / probs.sum()
    return [FluidState(t, float(p)) for (t, _), p in zip(states, probs) if p > 1e-12]


def project_to_simplex(v: np.ndarray) -> np.ndarray:
    """Euclidean projection onto {p : p >= 0, sum(p) = 1}."""
    u = np.sort(v)[::-1]
    css = np.cumsum(u) - 1.0
    k = np.arange(1, len(v) + 1)
    rho = np.nonzero(u - css / k > 0)[0][-1]
    return np.maximum(v - css[rho] / (rho + 1.0), 0.0)


def post_ordering_state_probabilities(types: list[tuple[int, ...]], P: np.ndarray,
                                      eps: float = 1e-8, n_iter: int = 20000) -> np.ndarray:
    """Fluid-state probabilities whose per-layer fluid abundances best match
    ``P`` *after* density ordering (paper Appendix 15).

    Solves min ||X p - P||^2 + eps*|p|^2 over the probability simplex, where
    X[(layer, fluid), state] = 1 if that state has that fluid in that layer.
    Inconsistent targets get the closest consistent answer.
    """
    nl = P.shape[0]
    X = np.zeros((nl * 4, len(types)))
    for k, t in enumerate(types):
        for i, f in enumerate(t):
            X[i * 4 + f, k] = 1.0
    target = P.reshape(-1)
    H = X.T @ X + eps * np.eye(len(types))
    g0 = X.T @ target
    step = 1.0 / np.linalg.eigvalsh(H)[-1]
    p = np.full(len(types), 1.0 / len(types))
    y, t_k = p.copy(), 1.0
    for _ in range(n_iter):  # accelerated projected gradient (FISTA)
        p_new = project_to_simplex(y - step * (H @ y - g0))
        t_new = 0.5 * (1.0 + np.sqrt(1.0 + 4.0 * t_k * t_k))
        y = p_new + ((t_k - 1.0) / t_new) * (p_new - p)
        if np.max(np.abs(p_new - p)) < 1e-14:
            p = p_new
            break
        p, t_k = p_new, t_new
    return p
