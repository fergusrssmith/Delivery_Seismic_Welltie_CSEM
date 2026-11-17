"""Discrete fluid states: which fluid sits in each layer.

Each permeable layer may hold brine, oil, low-saturation gas (lsg) or gas,
with prior probabilities set per layer.  A density-ordering rule then decides
which combinations across layers are geologically allowed (e.g. no oil above
gas in communicating sands).  The surviving combinations are the discrete
"fluid states" the inversion chooses between.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ModelConfig

BRINE, OIL, LSG, GAS = 0, 1, 2, 3
FLUID_NAMES = ("brine", "oil", "lsg", "gas")
FLUID_CHARS = "bolg"
# Density ordering from heaviest to lightest is brine, lsg, oil, gas.
DENSITY_PRIORITY = (0, 2, 1, 3)


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
        if lay.is_permeable:
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
    P = layer_fluid_probabilities(cfg)
    perm = [lay.is_permeable for lay in cfg.layers]
    nl = cfg.n_layers
    mode = cfg.density_ordering
    states: list[tuple[tuple[int, ...], float]] = []

    def recurse(layer: int, types: list[int], prob: float, lightest_allowed: int) -> None:
        if layer == nl:
            states.append((tuple(types), prob))
            return
        if not perm[layer]:
            # impermeable: always "brine"; resets ordering unless mode is full
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
        raise ValueError(
            "density ordering and fluid probabilities allow no legal fluid states"
        )
    total = sum(p for _, p in states)
    return [FluidState(t, p / total) for t, p in states]
