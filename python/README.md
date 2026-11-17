# pydelivery

A pure-Python port of the seismic-inversion part of **Delivery** (James
Gunning, CSIRO / BHP Billiton): Bayesian, trace-by-trace inversion of stacked
seismic into a layered earth model, with posterior samples of layer times,
net-to-gross, rock properties and pore fluids.

No Java and no Seismic Unix are needed. Data are read and written as HDF5,
NumPy `.npz` or SEG-Y; old Java projects (`.su` files) can be imported once.

The method is described in `../Docs/DeliveryPaper.pdf`. This port keeps the
geology of that model. It does not replace it with a generic smooth inversion.

## What makes Delivery different, and how the port keeps it

| Geological feature | How it is represented |
|---|---|
| **Discrete fluid states.** Each permeable layer holds brine, oil, low-saturation gas or gas, with per-layer prior probabilities | `fluids.enumerate_fluid_states` lists every combination allowed by the density-ordering rule (`none`, `partial`, `full`). Each combination is a separate model with its own parameter vector. Post-ordering probabilities (paper Appendix 15) are solved as a constrained least-squares fit. |
| **Fluid contacts in depth** (paper Appendix 11) | Fluids follow depth contacts in fluid units. Depths are hung from the reference layer; a layer crossed by a contact is split into sub-layers, each fluid-substituted separately, and the contact reflects. Contact depths can be uncertain; fluid properties can be shared across a unit. |
| **Pinch-outs.** Layers can thin to zero | Layer times have Gaussian priors but are put in order by truncation (best-picked horizons first). A pinched-out layer is skipped by the reflectivity, so a pinch-out has finite probability. A zero isopach N(0, 0) forces a permanent pinch-out. |
| **Truncated properties** (NG = 0 or 1, saturations at 0 or 1) | Clipped the same way, so boundary values carry probability mass. |
| **Rock physics** | Trend curves vs depth (or `depth_rock_curves`) and LFIV, Gassmann fluid substitution, Backus averaging of sand/shale laminates; optional second "hard" impermeable rock, fixed or fluid-dependent (Appendix 14); floating-grain terms (Appendix 9). |
| **Shared variables** | `same_as_layer_above` (depth, rock-curve depth, LFIV, float fraction), `synchronise_rock_matrix`, `synchronise_B_from_A`. |
| **Vertical model stacking** (Appendix 10) | Property-difference terms in the prior, a thickness-difference likelihood, and seismic-likelihood exclusion windows. |
| **Multi-modality** | `modes.py` finds local modes for every fluid state, including deliberately pinched-out starts, with Laplace masses. |
| **Jumps between fluid states** | `sampler.py`: reversible-jump, mode-jumping MCMC (paper section 4.3, Appendix 6). |

## Seismic stacks

Any number of PP stacks can be inverted together. A stack is defined either
as Delivery defines it, by an offset range plus a reflector time and stacking
velocity (the angle follows from Dix), or by an **angle range**, as partial
stacks are usually delivered:

```xml
<seismic_data reflectivity="aki-richards">
  <stack>
    <name>near</name><filename>near.sgy</filename>
    <min_angle>5</min_angle><max_angle>15</max_angle>
    <wavelet><filename>near_wavelet.h5</filename><noise_rms>0.01</noise_rms></wavelet>
  </stack>
  <stack>... mid 15-25 ...</stack>
  <stack>... far 25-35 ...</stack>
</seismic_data>
```

Each stack has its own wavelet and noise level (the noise can vary trace by
trace through `noise_rms_<stack name>` in the prior table). Stacks are put in
order near to far. `reflectivity` is `delivery` (the default: the original
small-angle form, θ² for both sin²θ and tan²θ) or `aki-richards` (sin²θ and
tan²θ averaged over the stack's angles), which is the better choice beyond
about 25°.

PS (converted-wave) stacks are not supported. Depth self-consistency is not
implemented; it is also switched off in the Java build.

## Install

```bash
cd python
pip install -e .[test]      # numpy and h5py
pytest
```

## Tutorial

`examples/synthetic_wedge.ipynb` walks through a complete run step by step:
build a 200-trace gas sand wedge that pinches out, make its seismic, invert it
from rough picks, and check it with the standard QC display (true impedance,
input seismic, estimated impedance, synthetic from the estimate, residual)
before looking at thickness, pinch-out and gas probabilities. It is committed
with its outputs; re-run it with `jupyter notebook examples/synthetic_wedge.ipynb`
(about five minutes on four cores).

`examples/fluid_contacts.ipynb` is a harder example from the original Java
distribution (`Simple3LayerContactsTestProblem`): a dipping sand crossing
gas-oil and oil-brine contacts. It rebuilds the known ground truth, inverts
with the contacts known, then starts from wrong contact depths to show the
seismic finding the true ones, and why depth control matters. It includes
the same QC display (about three minutes on four cores).

`examples/lateral_continuity.ipynb` (experimental) compares trace-by-trace
inversion of the wedge with `pydelivery.lateral`, which adds a spatially
correlated prior on rock properties along a line, the Bayesian counterpart of
spatial (Laplacian) regularisation. It removes most of the trace-to-trace
striping without changing the fit, the geometry or the fluid and pinch-out
answers. A K-means facies classification of the results (`pydelivery.qc.
classify_facies`) tests whether the striping changes the interpretation; here
it does not (about 13 minutes).

`examples/partial_stacks.ipynb` records the same wedge as three partial-angle
stacks (5-15°, 15-25°, 25-35°) and inverts them simultaneously, against the
5-15° stack alone. It shows what the far angles add: the shear contrast and
the sand's vp/vs become constrained by the data rather than by the
rock-physics prior, which is what makes an AI-SI cross-plot meaningful.

For your own QC, `pydelivery.qc.posterior_sections(results, times)` returns
posterior-mean impedance and synthetic sections, and `pydelivery.qc.plot_qc`
draws the five panels; `pydelivery.qc.plot_stack_qc` shows input, synthetic
and residual for every stack. `Inversion.run(..., n_jobs=-1)` inverts traces in
parallel.

## Data formats

| Kind | HDF5 (`.h5`) / NumPy (`.npz`) contents | SEG-Y |
|---|---|---|
| Seismic stack | `data` (n_traces, n_samples), `dt` (s), `t0` (s, per trace), optional `headers/…` | yes |
| Wavelet | same, one trace, `t0` < 0 (the precursor) | yes, if `t0` is a whole number of ms |
| Prior table (areally varying priors) | one array per XML tag: `time` (n_traces, n_layers + 1, ms, base last), `net_to_gross` (n_traces, n_layers), `contact_depth` (n_traces, n_contacts), `noise_rms_<stack>` or `near_stack_noise_rms` (n_traces,), … NaN = use the XML value | no |
| Results | `samples/<property>` (n_traces, n_samples, n_layers) or (n_traces, n_samples), `headers/…`, `fluid_states`, layer names | no |

The model itself stays in the Delivery `ModelDescription.xml` format; its
`<filename>` entries may name `.h5`, `.npz`, `.sgy` or `.su` files. SEG-Y is
read with IBM or IEEE floats and either byte order, and written as IEEE
big-endian revision 1.

## Run

```bash
# one-off: convert a Java project (XML + .su files) to HDF5
EX=../au/csiro/JamesGunning/Tests/InverterTests/Simple3LayerContactsTestProblem
pydelivery import-java $EX/ModelDescription.xml contacts_py -m $EX/prior_traces.su

pydelivery invert contacts_py/ModelDescription.xml -m contacts_py/priors.h5 -N 500 -o results.h5 -v
pydelivery summary results.h5 --property S_oil     # fluid / pinch-out probabilities, P10/50/90
pydelivery convert line.su line.sgy                # trace files between formats
```

From Python:

```python
from pydelivery import Inversion
inv = Inversion.from_files("ModelDescription.xml", "priors.h5")
results = inv.run(n_samples=500, seed=1)
r = results[0]
r.fluid_state_probabilities()   # e.g. {'bbb': 0.12, 'bgb': 0.88}
r.fluid_probabilities()         # (layers, [brine, oil, lsg, gas])
r.pinchout_probability()        # per layer
r.quantile("thickness", [0.1, 0.5, 0.9])
r.samples["z"]                  # layer-top depths hung from the reference layer
inv.save("results.h5", results)
```

Main output properties per sample: `t`, `z` (layer-top depth), `NG`, `vp_m`,
`vs_m`, `rho_m`, `vp_s`, `vs_s`, `phi_s`, `NG_hard`, `vp_hard`, `vs_hard`,
`rho_hard`, `float`, `vp_b`, `rho_b`, `vp_eff`, `vs_eff`, `rho_eff`, `Z_eff`,
`Zs_eff`, `S_h`, `S_oil`, `S_lsg`, `S_gas`, `fluid_type`, `R_<stack>` (reflection
at each layer top), `seismic_chisquare`, `t_base`, `z_base`, `A`, `B`, and in
contacts mode `contact_depth`, `n_sublayers`, `t_contact`, `R_contact_<stack>`.

## Validation against the Java code

`validation/compare_with_java.py` prints the posterior statistics of both
codes side by side; it needs only a Java realisations file. The codes use
different random numbers, so results agree statistically, not draw by draw.

* **Forward model:** reproduces Java-generated noise-free seismic to a misfit
  below 1e-4 on `NG_Wedge`, and to zero on the contacts example (traces where
  contacts split layers into two sub-layers).
* **NG_Wedge** (15 traces): NG means within 0.01 of Java; pinch-out
  probability of the wedge layer matches (e.g. 0.072 vs 0.072).
* **GasOilFlank** (brine → oil → gas up-dip, near + far stack): fluid
  probabilities match trace by trace (e.g. P(oil) 0.62 vs 0.60, P(gas) 1.00
  vs 1.00).
* **Simple3LayerContactsTestProblem** (contacts mode, gas over oil over
  brine): times, NG, gas and oil saturations and fluid probabilities match,
  e.g. trace 6 S_oil 0.493 vs 0.490, S_gas 0.307 vs 0.310.

Known differences: Java samples the far tails of layer times more widely
(1st/99th percentiles), and in a few contacts traces where a contact sits
near a layer boundary the saturation split differs by a few hundredths.
Medians and 10–90% ranges agree.

## Design notes and GPU

* **Linear forward operator.** For each trace and stack the synthetic at the
  misfit sample points is `G @ R`, with `G` built once; no FFT in the sampler.
* **Batched.** Every function takes a leading batch axis; the sampler runs 64
  chains at once. On CPU, 1000 samples take about 1–3 s per trace on the
  examples (contacts mode a few seconds more, because of the wider mode
  search). Java is faster per trace on these small examples.
* **GPU path.** The numerics are plain NumPy array operations, so they can move
  to JAX or CuPy without restructuring; the natural GPU batch is many traces
  times many chains. Not done yet: on the examples the cost is Python
  overhead, not arithmetic. The discrete fluid, contact and pinch-out
  structure stays as above in any GPU version; only continuous within-mode
  moves would gain from gradients.

## Licence

GPL v2 or later, as the original Delivery code (see `../LICENCE_GPL`).
Delivery and its method are by James Gunning; this is an independent port.
