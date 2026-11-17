from pathlib import Path

import numpy as np
import pytest

from pydelivery import Inversion, StateModel, enumerate_fluid_states, invert_trace, load_config
from pydelivery.config import apply_prior_trace
from pydelivery.seismic import build_seismic_model, load_wavelet
from pydelivery.su import Trace, read_su, write_su

REPO = Path(__file__).resolve().parents[2]
NG_WEDGE = REPO / "au/csiro/JamesGunning/Tests/InverterTests/NG_Wedge"

TWO_LAYER_SANDS = """<?xml version="1.0"?>
<inversion>
  <seismic_data><stack><name>near</name><filename>stack.su</filename>
    <min_offset>0</min_offset><max_offset>1200</max_offset>
    <reflector_time>1800</reflector_time><stack_velocity>2200</stack_velocity>
    <wavelet><filename>{wavelet}</filename><noise_rms>{noise}</noise_rms></wavelet></stack>
  </seismic_data>
  <model_description>
    <top_layer><name>seal</name>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time>1800</time><sigma_time>10</sigma_time><net_to_gross>0</net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
    </top_layer>
    <middle_layer><name>sand A</name>
      <reservoir_endmember><name>sand</name><brine><name>brine</name></brine>
        <oil><name>oil</name><probability_of_oil>{p_oil}</probability_of_oil><saturation_of_oil>0.8</saturation_of_oil></oil>
        <gas><name>gas</name><probability_of_gas>{p_gas}</probability_of_gas><saturation_of_gas>0.8</saturation_of_gas></gas>
      </reservoir_endmember>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time>1900</time><sigma_time>10</sigma_time>
      <net_to_gross>0.8</net_to_gross><sigma_net_to_gross>0.1</sigma_net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
    </middle_layer>
    <bottom_layer><name>base shale</name>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time>1960</time><sigma_time>10</sigma_time><net_to_gross>0</net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
      <t_base>2100</t_base><sigma_t_base>5</sigma_t_base>
    </bottom_layer>
  </model_description>
  <rock_fluid_properties>
    <fluid_properties>
      <fluid><name>brine</name><vp>1750</vp><sigma_vp>0</sigma_vp><rho>1050</rho><sigma_rho>0</sigma_rho></fluid>
      <fluid><name>oil</name><vp>1300</vp><sigma_vp>0</sigma_vp><rho>800</rho><sigma_rho>0</sigma_rho></fluid>
      <fluid><name>gas</name><vp>550</vp><sigma_vp>50</sigma_vp><rho>150</rho><sigma_rho>30</sigma_rho></fluid>
    </fluid_properties>
    <rock_properties>
      <reservoir_endmember><name>sand</name><reference_fluid>brine</reference_fluid>
        <vp_curve><intercept>3500</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>150</sigma></vp_curve>
        <vs_curve><intercept>-733</intercept><slope>0.768</slope><sigma>100</sigma></vs_curve>
        <porosity_curve><intercept>0.655</intercept><slope>-0.0001183</slope><sigma>0.02</sigma></porosity_curve>
        <vp_grain>6110</vp_grain><vs_grain>4131</vs_grain><rho_grain>2649</rho_grain>
      </reservoir_endmember>
      <nonreservoir_endmember><name>shale</name>
        <vp_curve><intercept>2812</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>100</sigma></vp_curve>
        <vs_curve><intercept>-496</intercept><slope>0.57</slope><sigma>43</sigma></vs_curve>
        <density_curve><factor>300</factor><exponent>0.26</exponent><sigma>20</sigma></density_curve>
      </nonreservoir_endmember>
    </rock_properties>
  </rock_fluid_properties>
  <output><density_ordering>full</density_ordering></output>
</inversion>
"""


def make_case(tmp_path, truth_fluid: str, p_oil=0.0, p_gas=0.5, noise=0.005):
    """Write a 3-layer model and a noise-free synthetic made with ``truth_fluid``
    in the sand (at the prior mean of everything else)."""
    xml = tmp_path / "model.xml"
    xml.write_text(TWO_LAYER_SANDS.format(wavelet=NG_WEDGE / "near-wavelet.su",
                                          noise=noise, p_oil=p_oil, p_gas=p_gas))
    blank = Trace({"dt": 4000, "delrt": 1700, "tracl": 1}, np.zeros(127))
    write_su(tmp_path / "stack.su", [blank])
    cfg = load_config(xml)
    w = load_wavelet(cfg.path(cfg.stacks[0].wavelet_filename))
    seis = build_seismic_model(cfg, [blank], [w])
    state = next(s for s in enumerate_fluid_states(cfg) if s.label[1] == truth_fluid)
    model = StateModel(cfg, state, seis)
    syn = model.synthetics(model.prior_mean[None])[0][0]
    data = np.zeros(127)
    j0 = int(round((seis.t_start - 1.7) / seis.dt))
    data[j0:j0 + seis.ns] = syn
    write_su(tmp_path / "stack.su", [Trace(blank.header, data)])
    return xml


def test_su_roundtrip(tmp_path):
    tr = read_su(NG_WEDGE / "near-stack.su")
    write_su(tmp_path / "x.su", tr)
    back = read_su(tmp_path / "x.su")
    assert len(back) == len(tr)
    assert back[3].header["delrt"] == tr[3].header["delrt"]
    np.testing.assert_allclose(back[3].data, tr[3].data, rtol=1e-6)


def test_fluid_state_enumeration_matches_paper(tmp_path):
    # Delivery paper, footnote 4: two oil-possible layers under full ordering
    xml = make_case(tmp_path, "b", p_oil=0.6, p_gas=0.0)
    cfg = load_config(xml)
    cfg.layers[2].sand, cfg.layers[2].ng, cfg.layers[2].brine = "sand", 0.8, "brine"
    cfg.layers[2].hydrocarbons = dict(cfg.layers[1].hydrocarbons)
    states = {s.label[1:]: s.probability for s in enumerate_fluid_states(cfg)}
    assert set(states) == {"bb", "ob", "oo"}  # no brine above oil
    assert states["bb"] == pytest.approx(0.16 / 0.76)
    assert states["oo"] == pytest.approx(0.36 / 0.76)


def test_time_ordering_follows_uncertainty_order(tmp_path):
    # Delivery paper section 2.2.1: sort order {1,4,2,3}, t1 < t4 < t3 < t2
    cfg = load_config(make_case(tmp_path, "b"))
    model = StateModel(cfg, enumerate_fluid_states(cfg)[0], None)
    model.time_order = np.array([0, 3, 1, 2])
    t = np.array([[1.0, 4.0, 3.0, 2.0]])
    np.testing.assert_allclose(model.order_times(t), [[1.0, 2.0, 2.0, 2.0]])


def test_forward_model_reproduces_noise_free_example():
    # NG_Wedge's seismic was generated by Delivery from its own prior mean
    inv = Inversion.from_files(NG_WEDGE / "ModelDescription.xml", NG_WEDGE / "prior_traces.su")
    for i in (0, 7, 14):
        cfg, seis = inv.trace_setup(i)
        model = StateModel(cfg, enumerate_fluid_states(cfg)[0], seis)
        m = model.prior_mean.copy()
        # the truth case was made with each trace's NG, not NG clipped at 1
        m[model.index[("NG", 1)]] = cfg.layers[1].ng
        ev = model.evaluate(m[None])
        assert ev.seismic_nll[0] < 1e-4


def test_pinchout_has_finite_probability(tmp_path):
    cfg = load_config(make_case(tmp_path, "b"))
    cfg.layers[1].time = 1.955  # sand 5 ms thick, sigma 10 ms: often pinched out
    r = invert_trace(cfg, None, 2000, np.random.default_rng(0), ignore_seismic=True)
    p = r.pinchout_probability()[1]
    assert 0.2 < p < 0.6
    thick = r.samples["t"][:, 2] - r.samples["t"][:, 1]
    assert np.all(thick >= 0)


@pytest.mark.parametrize("truth", ["g", "b"])
def test_seismic_identifies_fluid(tmp_path, truth):
    inv = Inversion.from_files(make_case(tmp_path, truth))
    r = inv.run(2000, seed=1)[0]
    p = r.fluid_probabilities()[1]
    expected = 3 if truth == "g" else 0
    # brine in softer sand and gas in stiffer sand give similar amplitudes, so
    # about 10% ambiguity remains whichever is the truth
    assert p[expected] > 0.8
    # the chain's fluid frequencies should agree with the Laplace mode weights
    laplace = {}
    for mo in r.modes:
        laplace[mo.model.state.label] = laplace.get(mo.model.state.label, 0.0) + mo.weight
    for label, prob in r.fluid_state_probabilities().items():
        assert prob == pytest.approx(laplace.get(label, 0.0), abs=0.05)


def test_prior_fluid_odds_without_seismic(tmp_path):
    inv = Inversion.from_files(make_case(tmp_path, "g"))
    cfg, seis = inv.trace_setup(0)
    r = invert_trace(cfg, seis, 3000, np.random.default_rng(2), ignore_seismic=True)
    assert r.fluid_probabilities()[1][3] == pytest.approx(0.5, abs=0.06)


def test_realisations_file_layout(tmp_path):
    inv = Inversion.from_files(NG_WEDGE / "ModelDescription.xml", NG_WEDGE / "prior_traces.su")
    res = inv.run(10, traces=[0], seed=0, burn_in=100, thin=5, n_chains=10)
    from pydelivery import write_realisations

    write_realisations(tmp_path / "r.su", res)
    tr = read_su(tmp_path / "r.su")
    assert len(tr) == 10 and tr[0].ns == 28 * 3 + 1  # same layout as Java Delivery
    assert tr[0].header["duse"] == 22
    assert 1700 < tr[0].data[-1] < 2300  # t_base, ms
