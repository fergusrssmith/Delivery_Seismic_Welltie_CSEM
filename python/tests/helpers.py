"""Builders for small synthetic Delivery models used by the tests."""
from pathlib import Path

import numpy as np

from pydelivery import StateModel, enumerate_fluid_states, load_config
from pydelivery.io import TraceSet, write_traces
from pydelivery.seismic import build_seismic_model, load_wavelet

REPO = Path(__file__).resolve().parents[2]
EXAMPLES = REPO / "au/csiro/JamesGunning/Tests/InverterTests"
NG_WEDGE = EXAMPLES / "NG_Wedge"
CONTACTS = EXAMPLES / "Simple3LayerContactsTestProblem"

SHALE = '<nonreservoir_endmember><name>{shale}</name></nonreservoir_endmember>'
DEPTH = ('<depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV>'
         '<sigma_LFIV>0</sigma_LFIV>')


def shale_layer(tag, name, time, sigma_time=10, extra="", shale="shale", attrs=""):
    return (f'<{tag}{attrs}><name>{name}</name>{SHALE.format(shale=shale)}'
            f'<time>{time}</time><sigma_time>{sigma_time}</sigma_time>'
            f'<net_to_gross>0</net_to_gross>{DEPTH}{extra}</{tag}>')


def sand_layer(tag, name, time, sigma_time=10, ng=0.8, sigma_ng=0.1, p_oil=0.0, p_gas=0.0,
               extra="", hc_extra="", sand="sand", attrs=""):
    hcs = ""
    if p_oil or "oil" in hc_extra:
        hcs += (f'<oil><name>oil</name><probability_of_oil>{p_oil}</probability_of_oil>'
                f'<saturation_of_oil>0.8</saturation_of_oil></oil>')
    if p_gas or "gas" in hc_extra:
        hcs += (f'<gas><name>gas</name><probability_of_gas>{p_gas}</probability_of_gas>'
                f'<saturation_of_gas>0.8</saturation_of_gas></gas>')
    return (f'<{tag}{attrs}><name>{name}</name>'
            f'<reservoir_endmember><name>{sand}</name><brine><name>brine</name></brine>{hcs}'
            f'</reservoir_endmember>{SHALE.format(shale="shale")}'
            f'<time>{time}</time><sigma_time>{sigma_time}</sigma_time>'
            f'<net_to_gross>{ng}</net_to_gross><sigma_net_to_gross>{sigma_ng}</sigma_net_to_gross>'
            f'{DEPTH}{extra}</{tag}>')


ROCKS = """
<reservoir_endmember><name>sand</name><reference_fluid>brine</reference_fluid>
  <vp_curve><intercept>3500</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>150</sigma>{vp_float}</vp_curve>
  <vs_curve><intercept>-733</intercept><slope>0.768</slope><sigma>100</sigma></vs_curve>
  <porosity_curve><intercept>0.655</intercept><slope>-0.0001183</slope><sigma>0.02</sigma>{phi_float}</porosity_curve>
  <vp_grain>6110</vp_grain><vs_grain>4131</vs_grain><rho_grain>2649</rho_grain>
</reservoir_endmember>
<nonreservoir_endmember><name>shale</name>
  <vp_curve><intercept>2812</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>100</sigma></vp_curve>
  <vs_curve><intercept>-496</intercept><slope>0.57</slope><sigma>43</sigma></vs_curve>
  <density_curve><factor>300</factor><exponent>0.26</exponent><sigma>20</sigma></density_curve>
</nonreservoir_endmember>
<nonreservoir_endmember><name>limestone</name>
  <vp_curve><intercept>5000</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>200</sigma></vp_curve>
  <vs_curve><intercept>0</intercept><slope>0.55</slope><sigma>100</sigma></vs_curve>
  <density_curve><factor>1700</factor><exponent>0.05</exponent><sigma>30</sigma></density_curve>
</nonreservoir_endmember>
<nonreservoir_endmember><name>dolomite</name>
  <vp_curve><intercept>5800</intercept><depth_coefficient>0</depth_coefficient><LFIV_coefficient>0</LFIV_coefficient><sigma>200</sigma></vp_curve>
  <vs_curve><intercept>0</intercept><slope>0.55</slope><sigma>100</sigma></vs_curve>
  <density_curve><factor>1900</factor><exponent>0.05</exponent><sigma>30</sigma></density_curve>
</nonreservoir_endmember>
"""

FLUIDS = """
<fluid><name>brine</name><vp>1750</vp><sigma_vp>0</sigma_vp><rho>1050</rho><sigma_rho>0</sigma_rho></fluid>
<fluid><name>oil</name><vp>1300</vp><sigma_vp>0</sigma_vp><rho>800</rho><sigma_rho>0</sigma_rho></fluid>
<fluid><name>gas</name><vp>550</vp><sigma_vp>50</sigma_vp><rho>150</rho><sigma_rho>30</sigma_rho></fluid>
"""


def model_xml(layers, *, t_base=2100, sigma_t_base=5, noise=0.005, stacks=None, avo="",
              extra="", output="<density_ordering>full</density_ordering>", vp_float="",
              phi_float="", wavelet=NG_WEDGE / "near-wavelet.su", angles=None, seismic_attrs=""):
    """Assemble a model XML. ``layers`` are layer XML strings; the last must be
    a bottom_layer without its t_base (added here)."""
    layers = list(layers)
    last = layers[-1]
    layers[-1] = last.replace("</bottom_layer>",
                              f"<t_base>{t_base}</t_base><sigma_t_base>{sigma_t_base}</sigma_t_base></bottom_layer>")
    if stacks is None:
        stacks = [""]
    stack_xml = ""
    if angles is not None:  # partial-angle stacks, e.g. [(5, 15), (15, 25)]
        for k, (a, b) in enumerate(angles):
            stack_xml += (f'<stack><name>a{a}-{b}</name><filename>stack{k}.h5</filename>'
                          f'<min_angle>{a}</min_angle><max_angle>{b}</max_angle>'
                          f'<wavelet><filename>{wavelet}</filename><noise_rms>{noise}</noise_rms></wavelet></stack>')
        stacks = []
    for k, attrs in enumerate(stacks):
        lo, hi = (0, 1200) if k == 0 else (1500, 3000)
        stack_xml += (f'<stack{attrs}><name>s{k}</name><filename>stack{k}.h5</filename>'
                      f'<min_offset>{lo}</min_offset><max_offset>{hi}</max_offset>'
                      f'<reflector_time>1800</reflector_time><stack_velocity>2200</stack_velocity>'
                      f'<wavelet><filename>{wavelet}</filename><noise_rms>{noise}</noise_rms></wavelet></stack>')
    return f"""<?xml version="1.0"?>
<inversion>
  <seismic_data{seismic_attrs}>{stack_xml}<AVO_terms{avo}><A>1.0</A><B>1.0</B></AVO_terms></seismic_data>
  <model_description>{''.join(layers)}</model_description>
  {extra}
  <rock_fluid_properties>
    <fluid_properties>{FLUIDS}</fluid_properties>
    <rock_properties>{ROCKS.format(vp_float=vp_float, phi_float=phi_float)}</rock_properties>
  </rock_fluid_properties>
  <output>{output}</output>
</inversion>
"""


def write_case(tmp_path, xml_text, truth_label=None, n_stacks=1, truth_m=None):
    """Write a model and seismic. With ``truth_label`` (a fluid-state label
    such as 'bgb'), the seismic is the noise-free synthetic of that state at
    its prior mean; otherwise it is zero."""
    xml = tmp_path / "model.xml"
    xml.write_text(xml_text)
    blank = TraceSet(np.zeros((1, 127)), 0.004, 1.7)
    for k in range(n_stacks):
        write_traces(tmp_path / f"stack{k}.h5", blank)
    if truth_label is None:
        return xml
    cfg = load_config(xml)
    w = [load_wavelet(cfg.path(s.wavelet_filename)) for s in cfg.stacks]
    seis = build_seismic_model(cfg, [blank[0]] * n_stacks, w)
    state = next(s for s in enumerate_fluid_states(cfg) if s.label == truth_label)
    model = StateModel(cfg, state, seis)
    m = model.prior_mean if truth_m is None else truth_m(model)
    j0 = int(round((seis.t_start - 1.7) / seis.dt))
    for k, syn in enumerate(model.synthetics(m[None])):
        data = np.zeros((1, 127))
        data[0, j0:j0 + seis.ns] = syn[0]
        write_traces(tmp_path / f"stack{k}.h5", TraceSet(data, 0.004, 1.7))
    return xml


def three_layer(p_oil=0.0, p_gas=0.5, **kw):
    return model_xml([
        shale_layer("top_layer", "seal", 1800),
        sand_layer("middle_layer", "sand A", 1900, p_oil=p_oil, p_gas=p_gas),
        shale_layer("bottom_layer", "base", 1960),
    ], **kw)
