"""Partial-angle stacks, any number of stacks, Aki-Richards reflectivity."""
import math

import numpy as np
import pytest

from helpers import three_layer, write_case
from pydelivery import Inversion, invert_trace, load_config
from pydelivery.io import write_prior_table
from pydelivery.rockphysics import reflection_coefficient

ANGLES = [(25, 35), (5, 15), (15, 25)]  # deliberately out of order


def test_stack_angle_terms_are_averages_over_the_range(tmp_path):
    cfg = load_config(write_case(tmp_path, three_layer(angles=[(5, 15)])))
    th2, s2, t2 = cfg.stacks[0].angle_terms()
    a = np.radians(np.linspace(5, 15, 200001))
    assert th2 == pytest.approx(np.mean(a ** 2), rel=1e-6)
    assert s2 == pytest.approx(np.mean(np.sin(a) ** 2), rel=1e-6)
    assert t2 == pytest.approx(np.mean(np.tan(a) ** 2), rel=1e-6)


def test_aki_richards_matches_textbook_form():
    rho1, vp1, vs1, rho2, vp2, vs2 = 2365.0, 2812.0, 1107.0, 2035.0, 2597.0, 1500.0
    th = math.radians(30)
    vp, vs, rho = (vp1 + vp2) / 2, (vs1 + vs2) / 2, (rho1 + rho2) / 2
    b2 = (vs / vp) ** 2
    textbook = (0.5 * (1 - 4 * b2 * math.sin(th) ** 2) * (rho2 - rho1) / rho
                + (vp2 - vp1) / (2 * vp * math.cos(th) ** 2)
                - 4 * b2 * math.sin(th) ** 2 * (vs2 - vs1) / vs)
    ours = reflection_coefficient(rho1, vp1, vs1, rho2, vp2, vs2, 1.0, math.sin(th) ** 2, math.tan(th) ** 2)
    assert float(ours) == pytest.approx(textbook, rel=1e-12)
    # Delivery's small-angle form is the theta -> 0 limit
    small = math.radians(2)
    d = reflection_coefficient(rho1, vp1, vs1, rho2, vp2, vs2, 1.0, small ** 2)
    ar = reflection_coefficient(rho1, vp1, vs1, rho2, vp2, vs2, 1.0, math.sin(small) ** 2, math.tan(small) ** 2)
    assert float(d) == pytest.approx(float(ar), abs=1e-6)


def test_three_partial_stacks_parse_sort_and_take_noise_by_name(tmp_path):
    xml = write_case(tmp_path, three_layer(angles=ANGLES, seismic_attrs=' reflectivity="aki-richards"'),
                     n_stacks=3)
    cfg = load_config(xml)
    assert [s.name for s in cfg.stacks] == ["a5-15", "a15-25", "a25-35"]  # near to far
    assert cfg.reflectivity == "aki-richards"
    write_prior_table(tmp_path / "pri.h5", {"noise_rms_a25-35": np.array([0.02])})
    inv = Inversion.from_files(xml, tmp_path / "pri.h5")
    cfg0, seis = inv.trace_setup(0)
    assert [s.noise_rms for s in cfg0.stacks] == [0.005, 0.005, 0.02]
    assert len(seis.stacks) == 3
    with pytest.raises(ValueError):
        write_prior_table(tmp_path / "bad.h5", {"noise_rms_nope": np.array([0.02])})
        Inversion.from_files(xml, tmp_path / "bad.h5").trace_setup(0)


def test_far_angles_constrain_shear_velocity(tmp_path):
    """A near stack alone barely constrains the shear contrast at the sand top;
    partial stacks out to 35 degrees do (vs only enters through contrasts)."""
    def dvs_sd(angles, sub):
        d = tmp_path / sub
        d.mkdir()
        text = three_layer(p_gas=0.0, noise=0.002, angles=angles,
                           seismic_attrs=' reflectivity="aki-richards"')
        xml = write_case(d, text, "bbb", n_stacks=len(angles))
        cfg, seis = Inversion.from_files(xml).trace_setup(0)
        r = invert_trace(cfg, seis, 600, np.random.default_rng(0), burn_in=400, thin=10)
        v = r.samples["vs_eff"]
        return np.std(v[:, 1] - v[:, 0])

    near = dvs_sd([(0, 8)], "near")
    partial = dvs_sd([(5, 15), (15, 25), (25, 35)], "partial")
    assert partial < 0.6 * near
