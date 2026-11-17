"""Effective elastic properties of a layer: Gassmann fluid substitution in the
permeable end-member, then Backus averaging with the impermeable end-member.

All functions are vectorised: property arrays have shape ``(batch, n_layers)``.
"""
from __future__ import annotations

import numpy as np

FOUR_THIRDS = 1.3333333333333  # as in the Java code


def effective_properties(
    p: dict[str, np.ndarray],
    K_g: np.ndarray,
    rho_g: np.ndarray,
    K_b_ref: np.ndarray,
    rho_b_ref: np.ndarray,
    has_hydrocarbon: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return ``vp_eff, vs_eff, rho_eff`` of shape (batch, n_layers) and a
    ``valid`` mask of shape (batch,).

    ``p`` holds the (already range-truncated) layer properties:
    NG, vp_m, vs_m, rho_m, vp_s, vs_s, phi_s, vp_b, rho_b, S_h, vp_h, rho_h.
    The trend curves describe the sand saturated with a *reference* brine
    (``K_b_ref``, ``rho_b_ref``); the in-situ pore fluid is a mix of the
    layer's brine and, where ``has_hydrocarbon`` is set, a hydrocarbon of
    saturation ``S_h``.  A sample is invalid when Gassmann's equation has no
    physical solution (negative dry or saturated bulk modulus).
    """
    NG, phi = p["NG"], p["phi_s"]
    vp_m, vs_m, rho_m = p["vp_m"], p["vs_m"], p["rho_m"]
    sand = NG > 0.0

    with np.errstate(divide="ignore", invalid="ignore"):
        # sand saturated with the reference brine
        rho_satb = phi * rho_b_ref + (1.0 - phi) * rho_g
        M_satb = rho_satb * p["vp_s"] ** 2
        mu_sat = rho_satb * p["vs_s"] ** 2
        K_satb = M_satb - FOUR_THIRDS * mu_sat

        bad = sand & (K_satb < 0.0) & (phi != 0.0)
        # dry-rock modulus must be positive
        K_sat_min = K_g / (phi * (K_g / K_b_ref - 1.0) + 1.0)
        bad |= sand & (K_satb < K_sat_min) & (phi != 0.0)

        # in-situ pore fluid: brine, or brine + hydrocarbon (Reuss mix)
        S = np.where(has_hydrocarbon, p["S_h"], 0.0)
        K_brine = p["rho_b"] * p["vp_b"] ** 2
        K_h = p["rho_h"] * p["vp_h"] ** 2
        K_fl = np.where(has_hydrocarbon, 1.0 / (S / K_h + (1.0 - S) / K_brine), K_brine)
        rho_fl = S * p["rho_h"] + (1.0 - S) * p["rho_b"]

        # Gassmann substitution, reference brine -> in-situ fluid
        K_satfl = K_g / (
            1.0
            + 1.0
            / (
                (1.0 / (K_g / K_fl - 1.0) - 1.0 / (K_g / K_b_ref - 1.0)) / phi
                + 1.0 / (K_g / K_satb - 1.0)
            )
        )
        K_satfl = np.where(phi != 0.0, K_satfl, K_g)
        bad |= sand & (K_satfl < 0.0)
        rho_satfl = (1.0 - phi) * rho_g + phi * rho_fl
        M_satfl = K_satfl + FOUR_THIRDS * mu_sat

        # Backus average of the sand / shale laminate
        M_m = rho_m * vp_m**2
        mu_m = rho_m * vs_m**2
        rho_mix = NG * rho_satfl + (1.0 - NG) * rho_m
        M_mix = 1.0 / (NG / M_satfl + (1.0 - NG) / M_m)
        mu_mix = 1.0 / (NG / mu_sat + (1.0 - NG) / mu_m)
        vp_mix = np.sqrt(M_mix / rho_mix)
        vs_mix = np.sqrt(mu_mix / rho_mix)

    vp_eff = np.where(sand, vp_mix, vp_m)
    vs_eff = np.where(sand, vs_mix, vs_m)
    rho_eff = np.where(sand, rho_mix, rho_m)
    finite = np.isfinite(vp_eff) & np.isfinite(vs_eff) & np.isfinite(rho_eff)
    valid = ~np.any(bad | ~finite, axis=-1)
    return vp_eff, vs_eff, rho_eff, valid


def reflection_coefficient(
    rho1, vp1, vs1, rho2, vp2, vs2, A, B, theta_sq_denominator
) -> np.ndarray:
    """Small-contrast P-P reflection coefficient to second order in angle.

    Layer 1 is above layer 2.  ``theta^2 = vp1^2 / theta_sq_denominator`` is
    the Dix estimate of the squared incidence angle for the stack; ``A`` and
    ``B`` are the AVO intercept and gradient scale factors.
    """
    vp_bar = 0.5 * (vp1 + vp2)
    vs_bar = 0.5 * (vs1 + vs2)
    d_rho = (rho2 - rho1) / (0.5 * (rho1 + rho2))
    d_vp = (vp2 - vp1) / vp_bar
    with np.errstate(divide="ignore", invalid="ignore"):
        d_vs = np.where(vs_bar > 0, (vs2 - vs1) / vs_bar, 0.0)
    beta = vs_bar / vp_bar
    theta_sq = B * vp1**2 / theta_sq_denominator
    return A * 0.5 * (d_rho + d_vp) + theta_sq * (
        0.5 * d_vp - 2.0 * beta**2 * (d_rho + 2.0 * d_vs)
    )
