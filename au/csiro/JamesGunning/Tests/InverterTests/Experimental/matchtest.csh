#!/bin/tcsh
# make  synthetics etc from prior + isopach only
java au.csiro.JamesGunning.seismicinverter.SeismicInverter ModelDescription.xml RockandFluidProperties.xml -v 4 -m prior_traces.su -N 1 -p -S
# collect synthetics
cat near_stack_synth.su | suswapbytes | suwind key=mark max=1 | suswapbytes format=1 >! 1near_stack_synth.su
ln -sf 1near_stack_synth.su near-stack.su

cat far_stack_synth.su | suswapbytes | suwind key=mark max=1 | suswapbytes format=1 >! 1far_stack_synth.su
ln -sf 1far_stack_synth.su far-stack.su


# then try to match -- omit tight isopachs to make it a challenge!
java au.csiro.JamesGunning.seismicinverter.SeismicInverter ModelDescription.xml RockandFluidProperties.xml -v 4 -m prior_traces.su -N 3000 -S -PD -RWS -II

# then see...
suswapbytes < near_stack_synth.su | suxgraph linewidth=4,1

