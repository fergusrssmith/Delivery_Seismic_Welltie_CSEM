"""One test per XML control added beyond the basic model."""
import numpy as np
import pytest

from helpers import model_xml, sand_layer, shale_layer, three_layer, write_case
from pydelivery import Inversion, StateModel, enumerate_fluid_states, invert_trace, load_config
from pydelivery.rockphysics import impermeable_properties


def state_model(cfg, label=None, seismic=None):
    states = enumerate_fluid_states(cfg)
    st = states[0] if label is None else next(s for s in states if s.label == label)
    return StateModel(cfg, st, seismic)


def two_sands(attrs2="", extra2="", sand2_kw=None, **kw):
    return model_xml([
        shale_layer("top_layer", "seal", 1800),
        sand_layer("middle_layer", "sand A", 1880, attrs=' synchronise_rock_matrix="true"'),
        shale_layer("middle_layer", "mid shale", 1920),
        sand_layer("middle_layer", "sand B", 1960, attrs=attrs2, extra=extra2, **(sand2_kw or {})),
        shale_layer("bottom_layer", "base", 2000),
    ], **kw)


def test_synchronise_rock_matrix(tmp_path):
    cfg = load_config(write_case(tmp_path, two_sands(attrs2=' synchronise_rock_matrix="true"')))
    assert cfg.layers[3].rock_master == 1
    m = state_model(cfg)
    keys = set(m.keys)
    assert ("vp_s", 1) in keys and ("vp_s", 3) not in keys and ("NG", 3) not in keys
    x = m.prior_mean + np.random.default_rng(0).standard_normal((5, m.dim)) @ np.linalg.cholesky(m.prior_cov).T
    p = m.tables(x)
    np.testing.assert_allclose(p["vp_s"][:, 1], p["vp_s"][:, 3])
    np.testing.assert_allclose(p["NG"][:, 1], p["NG"][:, 3])
    free = state_model(load_config(write_case(tmp_path, two_sands())))
    assert free.dim > m.dim


def test_same_as_layer_above_depth(tmp_path):
    xml = model_xml([
        shale_layer("top_layer", "a", 1800).replace("<sigma_depth>0", "<sigma_depth>20"),
        shale_layer("middle_layer", "b", 1900).replace("<depth>", '<depth same_as_layer_above="true">'),
        shale_layer("bottom_layer", "c", 1960),
    ])
    cfg = load_config(write_case(tmp_path, xml))
    assert [lay.depth_master for lay in cfg.layers] == [0, 0, 2]
    m = state_model(cfg)
    assert ("d", 0) in m.keys and ("d", 1) not in m.keys


def test_depth_rock_curves_drive_trends(tmp_path):
    xml = three_layer(p_gas=0.0).replace(
        "<sigma_LFIV>0</sigma_LFIV>", "<sigma_LFIV>0</sigma_LFIV><depth_rock_curves>1000</depth_rock_curves>", 1)
    xml = xml.replace("<depth_coefficient>0</depth_coefficient>", "<depth_coefficient>0.5</depth_coefficient>")
    cfg = load_config(write_case(tmp_path, xml))
    m = state_model(cfg)
    vp_m = m.prior_mean[m.index[("vp_m", 0)]], m.prior_mean[m.index[("vp_m", 1)]]
    assert vp_m[0] == pytest.approx(2812 + 0.5 * 1000)  # layer 1 uses depth_rock_curves
    assert vp_m[1] == pytest.approx(2812 + 0.5 * 2000)  # layer 2 uses its model depth


def test_forced_pinchout(tmp_path):
    xml = model_xml([
        shale_layer("top_layer", "a", 1800),
        sand_layer("middle_layer", "gone", 1880, extra="<thickness>0</thickness><sigma_thickness>0</sigma_thickness>"),
        shale_layer("bottom_layer", "c", 1960),
    ])
    cfg = load_config(write_case(tmp_path, xml))
    assert not cfg.layers[1].active
    r = invert_trace(cfg, None, 300, np.random.default_rng(0), ignore_seismic=True, burn_in=200, thin=5)
    assert r.pinchout_probability()[1] == 1.0
    assert ("t", 2) not in r.modes[0].model.keys  # base of the layer tied to its top


def test_hard_rock_backus_mix():
    vp, vs, rho = impermeable_properties(np.array(3000.0), np.array(1500.0), np.array(2400.0),
                                         np.array(0.5), np.array(5000.0), np.array(2700.0),
                                         np.array(2700.0), np.array(True))
    M = 1 / (0.5 / (2700 * 5000**2) + 0.5 / (2400 * 3000**2))
    assert rho == pytest.approx(2550.0)
    assert vp == pytest.approx(np.sqrt(M / 2550.0))


def test_fixed_and_fluid_dependent_hard_rock(tmp_path):
    fixed = ("<nonreservoir_hard_endmember><name>limestone</name><net_to_gross_hard>0.3</net_to_gross_hard>"
             "<sigma_net_to_gross_hard>0.1</sigma_net_to_gross_hard></nonreservoir_hard_endmember>")
    dep = ("<nonreservoir_hard_endmember><fluid_dependent_nonreservoir_hard_endmembers>"
           "<brine_case>limestone</brine_case><gas_case>dolomite</gas_case>"
           "</fluid_dependent_nonreservoir_hard_endmembers><net_to_gross_hard>0.3</net_to_gross_hard>"
           "<sigma_net_to_gross_hard>0.1</sigma_net_to_gross_hard></nonreservoir_hard_endmember>")
    for extra, expect in ((fixed, (5000, 5000)), (dep, (5000, 5800))):
        xml = model_xml([shale_layer("top_layer", "a", 1800),
                         sand_layer("middle_layer", "s", 1900, p_gas=0.5, extra=extra),
                         shale_layer("bottom_layer", "c", 1960)])
        cfg = load_config(write_case(tmp_path, xml))
        means = [state_model(cfg, lab).prior_mean[state_model(cfg, lab).index[("vp_hard", 1)]]
                 for lab in ("bbb", "bgb")]
        np.testing.assert_allclose(means, expect)
        m = state_model(cfg, "bbb")
        p = m.physical(m.prior_mean[None])
        assert p["vp_imp"][0, 1] > p["vp_m"][0, 1]  # hard rock stiffens the impermeable part


def test_floating_grain_trends(tmp_path):
    xml = three_layer(p_gas=0.0, vp_float="<phi_float_coefficient>-3000</phi_float_coefficient>",
                      phi_float="<phi_float_coefficient>-0.9</phi_float_coefficient>")
    xml = xml.replace("<sigma_net_to_gross>0.1</sigma_net_to_gross>",
                      "<sigma_net_to_gross>0.1</sigma_net_to_gross><float_fraction>0.05</float_fraction>"
                      "<sigma_float_fraction>0.05</sigma_float_fraction>")
    cfg = load_config(write_case(tmp_path, xml))
    m = state_model(cfg)
    mu = dict(zip(m.keys, m.prior_mean))
    assert mu[("float", 1)] == pytest.approx(0.05)
    assert mu[("vp_s", 1)] == pytest.approx(3500 - 3000 * 0.05)
    assert mu[("phi_s", 1)] == pytest.approx(0.655 - 0.0001183 * mu[("vp_s", 1)] - 0.9 * 0.05)
    # a negative float fraction is snapped to 0 by moving along the trend axes
    x = m.prior_mean.copy()
    x[m.index[("float", 1)]] = -0.02
    p = m.tables(x[None])
    assert p["float"][0, 1] == 0.0
    assert p["vp_s"][0, 1] == pytest.approx(x[m.index[("vp_s", 1)]] - 3000 * 0.02)


def stacked_xml(sigma_vp=50.0, sigma_thick=5.0, exclusion=True):
    layers = [shale_layer("top_layer", "a1", 1800), sand_layer("middle_layer", "s1", 1850),
              shale_layer("middle_layer", "buffer", 1900),
              shale_layer("middle_layer", "a2", 1950), sand_layer("middle_layer", "s2", 2000),
              shale_layer("bottom_layer", "c2", 2050)]
    excl = ("<seismic_likelihood_exclusion_range><start_layer>4</start_layer><end_layer>3</end_layer>"
            "</seismic_likelihood_exclusion_range>") if exclusion else ""
    extra = f"""<stacked_models_info>
      <layer_mappings><mapping>1:4</mapping><mapping>2:5</mapping></layer_mappings>
      <property_difference_std_deviations>
        <vp_m applicable_layers="1">{sigma_vp}</vp_m><thickness applicable_layers="2">{sigma_thick}</thickness>
      </property_difference_std_deviations>{excl}</stacked_models_info>"""
    return model_xml(layers, t_base=2150, extra=extra)


def test_vertical_stacking_prior_and_thickness(tmp_path):
    cfg = load_config(write_case(tmp_path, stacked_xml()))
    m = state_model(cfg)
    i, j = m.index[("vp_m", 0)], m.index[("vp_m", 3)]
    C = m.prior_cov
    sd_diff = np.sqrt(C[i, i] + C[j, j] - 2 * C[i, j])
    assert sd_diff < 50  # coupled: much tighter than independent (sqrt(2)*100)
    assert m.thickness_couplings == [(1, 4, 5.0)]
    i2, j2 = m.index[("vp_m", 1)], m.index[("vp_m", 4)]
    assert abs(C[i2, j2]) < 1e-9  # mapping 2:5 only couples thickness
    from pydelivery.model import StateModel as SM
    off = SM(cfg, enumerate_fluid_states(cfg)[0], None, use_stacking=False)
    assert off.prior_cov[i, j] == pytest.approx(0.0)


def test_seismic_exclusion_window(tmp_path):
    xml = write_case(tmp_path, stacked_xml())
    inv = Inversion.from_files(xml)
    cfg, seis = inv.trace_setup(0)
    sm = seis.stacks[0]
    assert (~sm.use).any() and sm.use.any()
    w = seis.stacks[0].wavelet
    lo = 1.95 - w.reach_after - 0.03
    hi = 1.95 + w.reach_before + 0.03
    np.testing.assert_array_equal(~sm.use, (sm.times > lo) & (sm.times < hi))


def test_avo_sync_and_constant_angle(tmp_path):
    xml = three_layer(avo=' synchronise_B_from_A="true"', stacks=[' is-constant-angle-stack-with-angle-degrees="30"'])
    xml = xml.replace("<A>1.0</A>", "<A>1.0</A><sigma_A>0.1</sigma_A>")
    cfg = load_config(write_case(tmp_path, xml))
    assert cfg.sync_B_from_A and cfg.stacks[0].angle == pytest.approx(np.radians(30))
    inv = Inversion.from_files(tmp_path / "model.xml")
    _, seis = inv.trace_setup(0)
    m = state_model(cfg, "bbb", seis)
    assert m.A_index == m.B_index >= 0
    np.testing.assert_allclose(seis.stacks[0].theta_sq(np.array([3000.0]), np.array([1.0])),
                               np.radians(30) ** 2)


def test_fluid_property_overrides(tmp_path):
    xml = three_layer().replace(
        "<saturation_of_gas>0.8</saturation_of_gas>",
        "<saturation_of_gas>0.8</saturation_of_gas><vp_gas>700</vp_gas><sigma_vp_gas>20</sigma_vp_gas>")
    xml = xml.replace("<name>brine</name></brine>",
                      "<name>brine</name><rho_brine>1100</rho_brine><sigma_rho_brine>10</sigma_rho_brine></brine>", 1)
    cfg = load_config(write_case(tmp_path, xml))
    m = state_model(cfg, "bgb")
    mu = dict(zip(m.keys, m.prior_mean))
    assert mu[("vp_gas", 1)] == pytest.approx(700)
    assert mu[("rho_b", 1)] == pytest.approx(1100)


def test_post_ordering_from_xml(tmp_path):
    layers = [shale_layer("top_layer", "a", 1800),
              sand_layer("middle_layer", "s1", 1850, p_oil=0.5, p_gas=0.5),
              sand_layer("middle_layer", "s2", 1900, p_oil=0.5, p_gas=0.5),
              shale_layer("bottom_layer", "c", 1960)]
    out = ("<density_ordering>full</density_ordering>"
           "<take_fluid_probabilities_as_post_ordering>true</take_fluid_probabilities_as_post_ordering>")
    cfg = load_config(write_case(tmp_path, model_xml(layers, output=out)))
    probs = {s.label: s.probability for s in enumerate_fluid_states(cfg)}
    gas_top = sum(p for k, p in probs.items() if k[1] == "g")
    gas_bottom = sum(p for k, p in probs.items() if k[2] == "g")
    assert gas_top == pytest.approx(0.5, abs=1e-4) and gas_bottom == pytest.approx(0.5, abs=1e-4)


def test_master_depth_layer_hanging(tmp_path):
    xml = three_layer(output="<density_ordering>full</density_ordering>"
                             "<master_depth_layer_number>2</master_depth_layer_number>")
    cfg = load_config(write_case(tmp_path, xml))
    m = state_model(cfg)
    p = m.physical(m.prior_mean[None])
    assert p["z"][0, 1] == pytest.approx(2000)  # layer 2 top sits at its depth
    assert p["z"][0, 0] == pytest.approx(2000 - 0.5 * 0.1 * p["vp_eff"][0, 0])


def test_legacy_signal_to_noise_and_ps_rejected(tmp_path):
    xml = three_layer().replace("noise_rms>", "signal_to_noise>")
    assert load_config(write_case(tmp_path, xml)).stacks[0].noise_rms == 0.005
    with pytest.raises(NotImplementedError):
        load_config(write_case(tmp_path, three_layer(stacks=[' is-PS-data="true"'])))
