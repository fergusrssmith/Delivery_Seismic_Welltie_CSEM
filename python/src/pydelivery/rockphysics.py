"""Effective elastic properties of a layer: Gassmann fluid substitution in the
permeable end-member, then Backus averaging with the impermeable end-member
(which may itself be a Backus mix of shale and a second "hard" rock).

All functions are vectorised over leading axes.
"""
from __future__ import annotations

import numpy as np

FOUR_THIRDS = 1.3333333333333  # as in the Java code


def impermeable_properties(vp_m, vs_m, rho_m, ng_hard, vp_h, vs_h, rho_h, has_hard):
    """Backus mix of the shale and the hard impermeable rock.

    ``ng_hard`` is the hard-rock volume fraction of the impermeable part.
    Layers without a hard rock (``has_hard`` False) return the shale.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        rho = ng_hard * rho_h + (1.0 - ng_hard) * rho_m
        M = 1.0 / (ng_hard / (rho_h * vp_h**2) + (1.0 - ng_hard) / (rho_m * vp_m**2))
        mu = 1.0 / (ng_hard / (rho_h * vs_h**2) + (1.0 - ng_hard) / (rho_m * vs_m**2))
        vp = np.sqrt(M / rho)
        vs = np.sqrt(mu / rho)
    use = has_hard & (ng_hard > 0.0)
    pure = has_hard & (ng_hard >= 1.0)
    vp = np.where(pure, vp_h, np.where(use, vp, vp_m))
    vs = np.where(pure, vs_h, np.where(use, vs, vs_m))
    rho = np.where(pure, rho_h, np.where(use, rho, rho_m))
    return vp, vs, rho


def effective_properties(NG, phi, vp_s, vs_s, vp_imp, vs_imp, rho_imp,
                         vp_b, rho_b, S_h, vp_h, rho_h, has_hc,
                         K_g, rho_g, K_b_ref, rho_b_ref):
    """Effective ``vp, vs, rho`` of a sand/impermeable laminate and a
    ``valid`` mask (same shape).

    The trend curves describe the sand saturated with a *reference* brine
    (``K_b_ref``, ``rho_b_ref``); the in-situ pore fluid is the layer's brine
    mixed, where ``has_hc``, with a hydrocarbon of saturation ``S_h``
    (Reuss average). Invalid where Gassmann has no physical solution
    (negative dry or saturated bulk modulus).
    """
    sand = NG > 0.0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        rho_satb = phi * rho_b_ref + (1.0 - phi) * rho_g
        M_satb = rho_satb * vp_s**2
        mu_sat = rho_satb * vs_s**2
        K_satb = M_satb - FOUR_THIRDS * mu_sat

        bad = sand & (K_satb < 0.0) & (phi != 0.0)
        K_sat_min = K_g / (phi * (K_g / K_b_ref - 1.0) + 1.0)
        bad |= sand & (K_satb < K_sat_min) & (phi != 0.0)

        S = np.where(has_hc, S_h, 0.0)
        K_brine = rho_b * vp_b**2
        K_hc = rho_h * vp_h**2
        K_fl = np.where(has_hc, 1.0 / (S / K_hc + (1.0 - S) / K_brine), K_brine)
        rho_fl = S * rho_h + (1.0 - S) * rho_b

        K_satfl = K_g / (1.0 + 1.0 / ((1.0 / (K_g / K_fl - 1.0) - 1.0 / (K_g / K_b_ref - 1.0)) / phi
                                      + 1.0 / (K_g / K_satb - 1.0)))
        K_satfl = np.where(phi != 0.0, K_satfl, K_g)
        bad |= sand & (K_satfl < 0.0)
        rho_satfl = (1.0 - phi) * rho_g + phi * rho_fl
        M_satfl = K_satfl + FOUR_THIRDS * mu_sat

        M_imp = rho_imp * vp_imp**2
        mu_imp = rho_imp * vs_imp**2
        rho_mix = NG * rho_satfl + (1.0 - NG) * rho_imp
        M_mix = 1.0 / (NG / M_satfl + (1.0 - NG) / M_imp)
        mu_mix = 1.0 / (NG / mu_sat + (1.0 - NG) / mu_imp)
        vp_mix = np.sqrt(M_mix / rho_mix)
        vs_mix = np.sqrt(mu_mix / rho_mix)

    vp = np.where(sand, vp_mix, vp_imp)
    vs = np.where(sand, vs_mix, vs_imp)
    rho = np.where(sand, rho_mix, rho_imp)
    valid = ~bad & np.isfinite(vp) & np.isfinite(vs) & np.isfinite(rho) & (rho > 0)
    return vp, vs, rho, valid


def reflection_coefficient(rho1, vp1, vs1, rho2, vp2, vs2, A, sin_sq, tan_sq=None):
    """Small-contrast (linearised) P-P reflection coefficient.

    Layer 1 is above layer 2::

        R = A * (drho + dvp) / 2 + tan2 * dvp / 2 - 2 beta^2 sin2 (drho + 2 dvs)

    with relative contrasts d* and beta = vs/vp (interface averages). With
    ``sin_sq = tan_sq = theta^2`` this is Delivery's small-angle form (Gunning
    & Glinsky eq. 22); with the true sin^2 and tan^2 it is the Aki-Richards
    linearisation, better for angles beyond about 25-30 degrees. Both are
    already scaled by the AVO gradient factor B. ``A`` scales the intercept.
    """
    if tan_sq is None:
        tan_sq = sin_sq
    with np.errstate(divide="ignore", invalid="ignore"):
        vp_bar = 0.5 * (vp1 + vp2)
        vs_bar = 0.5 * (vs1 + vs2)
        d_rho = (rho2 - rho1) / (0.5 * (rho1 + rho2))
        d_vp = (vp2 - vp1) / vp_bar
        d_vs = np.where(vs_bar > 0, (vs2 - vs1) / vs_bar, 0.0)
        beta = vs_bar / vp_bar
        r = A * 0.5 * (d_rho + d_vp) + 0.5 * tan_sq * d_vp - 2.0 * beta**2 * sin_sq * (d_rho + 2.0 * d_vs)
    return np.nan_to_num(r)
