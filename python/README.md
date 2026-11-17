# pydelivery

A Python port of the seismic-inversion part of **Delivery** (James Gunning,
CSIRO / BHP Billiton): Bayesian, trace-by-trace inversion of stacked seismic
into a layered earth model, with posterior samples of layer times,
net-to-gross, rock properties and pore fluids.

The method is described in `../Docs/DeliveryPaper.pdf`. This port keeps the
geology of that model. It does not replace it with a generic smooth inversion.

## What makes Delivery different, and how the port keeps it

| Geological feature | How it is represented |
|---|---|
| **Discrete fluid states.** Each permeable layer holds brine, oil, low-saturation gas or gas, with per-layer prior probabilities | `fluids.enumerate_fluid_states` lists every fluid combination allowed by the density-ordering rule (`none`, `partial`, `full`; no oil under brine, no gas under oil in communicating sands). Each combination is a separate model with its own parameter vector (gas properties exist only where there is gas). |
| **Pinch-outs.** Layers can thin to zero | Layer times have Gaussian priors but are put in order by truncation (`StateModel.order_times`): best-picked horizons are fixed first, and a layer whose base would cross its top gets zero thickness. The reflectivity skips pinched-out layers and takes the contrast between the layers either side. A pinch-out therefore has finite prior and posterior probability, not zero. |
| **Truncated properties.** NG = 0 (pure shale), NG = 1, saturations at 0 or 1 | Clipped the same way, so boundary values carry probability mass, as in the Java code. |
| **Rock physics.** Trend curves vs depth/LFIV, Gassmann fluid substitution, Backus averaging of sand/shale laminates | `rockphysics.py`; prior correlations between vp, vs and porosity/density come from the trend curves (`StateModel` builds them as linear residuals). |
| **Multi-modality.** Loop skips, present vs pinched-out layers, competing fluids | `modes.py` finds local modes for every fluid state, including from deliberately pinched-out starting times (paper Appendix 3), and gives each a Laplace mass. |
| **Jumps between fluid states.** | `sampler.py`: mode-jumping MCMC (paper section 4.3 and Appendix 6). Shared parameters are translated between mode peaks; parameters new to the target state are drawn from their prior; the reversible-jump acceptance ratio accounts for the change in dimension. |

## Install and run

```bash
cd python
pip install -e .[test]
pytest

# invert the NG wedge example shipped with the Java code
EX=../au/csiro/JamesGunning/Tests/InverterTests/NG_Wedge
pydelivery invert $EX/ModelDescription.xml -m $EX/prior_traces.su -N 500 -o realisations.su -v
pydelivery summary realisations.su        # fluid and pinch-out probabilities per trace
```

From Python:

```python
from pydelivery import Inversion, write_realisations
inv = Inversion.from_files("ModelDescription.xml", "prior_traces.su")
results = inv.run(n_samples=500, seed=1)
r = results[0]
r.fluid_state_probabilities()   # e.g. {'bbb': 0.12, 'bgb': 0.88}
r.fluid_probabilities()         # (layers, [brine, oil, lsg, gas])
r.pinchout_probability()        # per layer
r.quantile("thickness", [0.1, 0.5, 0.9])
write_realisations("realisations.su", results)
```

Inputs are the Java code's own inputs: the same `ModelDescription.xml`, SU
seismic stacks, SU wavelets and SU prior traces. The output
`realisations.su` uses the Java Delivery block layout (format version 22), so
existing SU-based tooling can read it.

## Validation against the Java code

`validation/compare_with_java.py` runs both codes on one example and prints
the statistics side by side. The two codes use different random numbers, so
results agree statistically, not draw by draw.

* **Forward model:** reproduces the noise-free `NG_Wedge` seismic (generated
  by Java Delivery) to a misfit below 1e-4 (`test_forward_model_reproduces_noise_free_example`).
* **NG_Wedge**, 15 traces, 2000 samples each: posterior NG mean within 0.01 of
  Java on every trace; top-time means within a few ms; pinch-out probability of
  the wedge layer matches (e.g. 0.071 vs 0.072, 0.029 vs 0.025).
* **GasOilFlank** (brine to oil to gas up-dip, near + far stack): posterior
  fluid probabilities match Java trace by trace, e.g. P(oil) 0.60 vs 0.60 in
  the oil zone, P(gas) 1.00 vs 1.00 in the gas zone.

Known difference: in some well-constrained traces the Python chain gives a
narrower spread of layer times than Java (e.g. trace 12 of NG_Wedge: 3.6 ms vs
7.0 ms standard deviation). The mean, NG and fluid results are unaffected.
Longer chains (`--thin`, `--burn-in`) or `--random-starts` reduce the gap.

## Not ported yet

The XML parser refuses these features rather than ignoring them:
internal fluid contacts (`fluid_units_and_contacts`), vertical model
stacking, hard-rock (second impermeable) end-members, floating-grain models,
`synchronise_rock_matrix`, `same_as_layer_above` on times/NG, forced
pinch-outs via a zero isopach, seismic exclusion zones, post-ordering fluid
probabilities, PS (converted-wave) data. Also not ported: the
`deliveryAnalyser` tool (beyond `pydelivery summary` and the `TraceResult`
statistics), the independence sampler, `--make-synthetic`, and the
`waveletExtractor` (planned next).

## Design notes and GPU

* **Linear forward operator.** Convolution is linear, so for each trace and
  stack the synthetic at the misfit sample points is `G @ R`, with `G` built
  once. No FFT is needed inside the sampler.
* **Batched.** Every function takes a leading batch axis. The sampler runs
  64 chains at once and evaluates all chains in the same mode in one call.
  On CPU, 2000 samples take about 2 s per trace on NG_Wedge and about 8 s
  on GasOilFlank (three fluid states, two stacks). Java takes about 0.6 s per
  trace on NG_Wedge, so the Java code is faster on these small examples.
* **GPU path.** The numerics are plain array operations (NumPy), chosen so
  they can move to JAX or CuPy without restructuring. On a GPU the natural
  batch is many traces times many chains. That has not been done yet: on the
  shipped examples the cost is dominated by Python overhead, not arithmetic,
  so a GPU would not help until traces are batched together. The discrete
  fluid/pinch-out structure stays exactly as above in any GPU version; only
  the continuous within-mode moves would gain from gradients (e.g. HMC),
  and the mode jumps would remain reversible-jump moves.

## Licence

GPL v2 or later, as the original Delivery code (see `../LICENCE_GPL`).
Delivery and its method are by James Gunning; this is an independent port.
