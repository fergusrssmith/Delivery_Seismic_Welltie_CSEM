import numpy as np
import pytest

from helpers import CONTACTS, model_xml, sand_layer, shale_layer, write_case
from pydelivery import Inversion, StateModel, enumerate_fluid_states, load_config
from pydelivery.modes import find_modes


def test_contacts_forward_model_matches_java_data():
    # Java made this seismic from its prior+isopach peak; at that peak our
    # contact-split synthetic must reproduce it
    inv = Inversion.from_files(CONTACTS / "ModelDescription.xml", CONTACTS / "prior_traces.su")
    rng = np.random.default_rng(0)
    for i, n_sub in ((3, [1, 2, 2]), (7, [2, 1, 2])):
        cfg, seis = inv.trace_setup(i)
        model = StateModel(cfg, enumerate_fluid_states(cfg)[0], seis)
        peak = find_modes(model, 0, rng, ignore_seismic=True)[0].mean
        assert model.evaluate(peak[None]).seismic_nll[0] < 1e-6
        assert list(model.physical(peak[None])["n_sublayers"][0]) == n_sub


def test_contacts_split_layer_and_upscale():
    inv = Inversion.from_files(CONTACTS / "ModelDescription.xml", CONTACTS / "prior_traces.su")
    cfg, seis = inv.trace_setup(5)
    model = StateModel(cfg, enumerate_fluid_states(cfg)[0], seis)
    p = model.physical(model.prior_mean[None])
    k = int(np.argmax(p["n_sublayers"][0]))
    t0, t1 = p["sub_t0"][0, k], p["sub_t1"][0, k]
    n = p["n_sublayers"][0, k]
    assert n >= 2
    assert t0[0] == pytest.approx(p["t"][0, k]) and t1[n - 1] == pytest.approx(p["t"][0, k + 1])
    np.testing.assert_allclose(t1[:n - 1], t0[1:n])  # sub-layers tile the layer
    # depth at the internal contact equals the contact depth
    vp = [model.physical(model.prior_mean[None])["iv_vp"][0, k * 4 + s] for s in range(n)]
    z = p["z"][0, k] + 0.5 * (t1[0] - t0[0]) * vp[0]
    assert z == pytest.approx(min(c.depth for c in cfg.contacts if c.depth > p["z"][0, k]), rel=1e-9)


def test_contacts_posterior_saturation():
    inv = Inversion.from_files(CONTACTS / "ModelDescription.xml", CONTACTS / "prior_traces.su")
    r = inv.run(600, traces=[0], seed=3)[0]
    # trace 1: the reservoir sits in the gas leg
    assert r.mean("S_gas")[1] == pytest.approx(0.8, abs=0.02)
    assert r.mean("S_oil")[1] < 0.02


def test_contacts_synchronised_fluids(tmp_path):
    layers = [shale_layer("top_layer", "seal", 1800),
              sand_layer("middle_layer", "s1", 1850, hc_extra="gas oil"),
              shale_layer("middle_layer", "m", 1900),
              sand_layer("middle_layer", "s2", 1950, hc_extra="gas oil"),
              shale_layer("bottom_layer", "c", 2000)]
    for sync, expect_shared in (("true", True), ("false", False)):
        extra = f"""<fluid_units_and_contacts><fluid_unit fluid_properties_synchronised="{sync}">
          <top_layer>1</top_layer><bottom_layer>5</bottom_layer>
          <contact><type>GO</type><depth>2080</depth><sigma_depth>5</sigma_depth></contact>
          <contact><type>OB</type><depth>2150</depth><sigma_depth>5</sigma_depth></contact>
          </fluid_unit></fluid_units_and_contacts>"""
        cfg = load_config(write_case(tmp_path, model_xml(layers, extra=extra)))
        assert cfg.contacts_mode and len(enumerate_fluid_states(cfg)) == 1
        m = StateModel(cfg, enumerate_fluid_states(cfg)[0], None)
        assert ("contact", 0) in m.keys and ("contact", 1) in m.keys
        assert (("vp_gas", 3) not in m.keys) == expect_shared
