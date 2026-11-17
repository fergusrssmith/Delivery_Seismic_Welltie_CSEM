"""Build the synthetic gas-sand wedge used by the example notebooks.

``build_wedge(work)`` writes the model XML, wavelet, seismic and prior table
for the 200-trace wedge of ``synthetic_wedge.ipynb`` (same random seed, so the
data are identical) and returns the true model for comparison.
"""
from pathlib import Path

import numpy as np

from pydelivery import StateModel, enumerate_fluid_states, load_config
from pydelivery.config import apply_prior_trace
from pydelivery.io import TraceSet, write_prior_table, write_traces
from pydelivery.seismic import build_seismic_model, load_wavelet

XML = """<?xml version="1.0"?>
<inversion>
  <seismic_data>
    <stack>
      <name>near</name>
      <filename>stack.h5</filename>
      <min_offset>0</min_offset><max_offset>1000</max_offset>
      <reflector_time>1900</reflector_time><stack_velocity>2500</stack_velocity>
      <wavelet><filename>wavelet.h5</filename><noise_rms>0.01</noise_rms></wavelet>
    </stack>
    <AVO_terms><A>1.0</A><B>1.0</B></AVO_terms>
  </seismic_data>

  <model_description>
    <top_layer>
      <name>seal</name>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time>1800</time><sigma_time>5</sigma_time>
      <net_to_gross>0</net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
    </top_layer>
    <middle_layer>
      <name>reservoir sand</name>
      <reservoir_endmember>
        <name>sand</name>
        <brine><name>brine</name></brine>
        <gas><name>gas</name>
          <probability_of_gas>0.5</probability_of_gas>
          <saturation_of_gas>0.8</saturation_of_gas>
        </gas>
      </reservoir_endmember>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time varies_areally="true">1900</time><sigma_time>5</sigma_time>
      <net_to_gross>0.8</net_to_gross><sigma_net_to_gross>0.1</sigma_net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
    </middle_layer>
    <bottom_layer>
      <name>base shale</name>
      <nonreservoir_endmember><name>shale</name></nonreservoir_endmember>
      <time varies_areally="true">1940</time><sigma_time>10</sigma_time>
      <net_to_gross>0</net_to_gross>
      <depth>2000</depth><sigma_depth>0</sigma_depth><LFIV>2000</LFIV><sigma_LFIV>0</sigma_LFIV>
      <t_base>2050</t_base><sigma_t_base>5</sigma_t_base>
    </bottom_layer>
  </model_description>

  <rock_fluid_properties>
    <fluid_properties>
      <fluid><name>brine</name><vp>1750</vp><sigma_vp>0</sigma_vp><rho>1050</rho><sigma_rho>0</sigma_rho></fluid>
      <fluid><name>gas</name><vp>550</vp><sigma_vp>50</sigma_vp><rho>150</rho><sigma_rho>30</sigma_rho></fluid>
    </fluid_properties>
    <rock_properties>
      <reservoir_endmember><name>sand</name><reference_fluid>brine</reference_fluid>
        <vp_curve><intercept>3000</intercept><depth_coefficient>0</depth_coefficient>
                  <LFIV_coefficient>0</LFIV_coefficient><sigma>150</sigma></vp_curve>
        <vs_curve><intercept>-733</intercept><slope>0.768</slope><sigma>100</sigma></vs_curve>
        <porosity_curve><intercept>0.655</intercept><slope>-0.0001183</slope><sigma>0.02</sigma></porosity_curve>
        <vp_grain>6110</vp_grain><vs_grain>4131</vs_grain><rho_grain>2649</rho_grain>
      </reservoir_endmember>
      <nonreservoir_endmember><name>shale</name>
        <vp_curve><intercept>2812</intercept><depth_coefficient>0</depth_coefficient>
                  <LFIV_coefficient>0</LFIV_coefficient><sigma>100</sigma></vp_curve>
        <vs_curve><intercept>-496</intercept><slope>0.57</slope><sigma>43</sigma></vs_curve>
        <density_curve><factor>300</factor><exponent>0.26</exponent><sigma>20</sigma></density_curve>
      </nonreservoir_endmember>
    </rock_properties>
  </rock_fluid_properties>

  <output><density_ordering>full</density_ordering></output>
</inversion>
"""


def build_wedge(work, n_traces=200, seed=42, noise=0.01):
    """Write the wedge project into ``work``; return a dict with the truth."""
    work = Path(work)
    work.mkdir(exist_ok=True)
    rng = np.random.default_rng(seed)
    dt, f = 0.004, 30.0
    t_w = np.arange(-16, 17) * dt
    wavelet = (1 - 2 * (np.pi * f * t_w) ** 2) * np.exp(-(np.pi * f * t_w) ** 2)
    write_traces(work / "wavelet.h5", TraceSet(wavelet[None, :], dt, t_w[0]))
    (work / "ModelDescription.xml").write_text(XML)
    cfg = load_config(work / "ModelDescription.xml")

    x = np.arange(1, n_traces + 1)
    top_true = np.full(n_traces, 1.900)
    thick_true = np.clip(0.040 * (1 - (x - 1) / 129), 0, None)
    base_true = top_true + thick_true
    top_pick = top_true + rng.normal(0, 0.003, n_traces)
    base_pick = base_true + rng.normal(0, 0.006, n_traces)
    times_table = np.full((n_traces, 4), np.nan)
    times_table[:, 1] = 1000 * top_pick
    times_table[:, 2] = 1000 * base_pick
    write_prior_table(work / "priors.h5", {"time": times_table})

    table = {"time": times_table}
    w = load_wavelet(work / "wavelet.h5")
    t0, ns = 1.700, 150
    times = t0 + dt * np.arange(ns)
    clean = np.zeros((n_traces, ns))
    Z_true = np.zeros((n_traces, ns))
    for i in range(n_traces):
        cfg_i = apply_prior_trace(cfg, table, i)
        blank = TraceSet(np.zeros((1, ns)), dt, t0)[0]
        seis = build_seismic_model(cfg_i, [blank], [w])
        gas = next(s for s in enumerate_fluid_states(cfg_i) if s.label == "bgb")
        model = StateModel(cfg_i, gas, seis)
        m = model.prior_mean.copy()
        m[model.index[("t", 1)]] = top_true[i]
        m[model.index[("t", 2)]] = base_true[i]
        tr = model.time_traces(m[None], times)
        clean[i], Z_true[i] = tr["synthetic_0"][0], tr["Zp"][0]
    data = clean + rng.normal(0, noise, clean.shape)
    write_traces(work / "stack.h5", TraceSet(data, dt, t0))
    return dict(x=x, times=times, data=data, Z_true=Z_true, thick_true=thick_true,
                top_true=top_true, base_true=base_true, times_table=times_table)
