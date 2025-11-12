#!/usr/bin/env bash
#
# requires some or all of the all the usual delivery arguments, e.g. -PD -RWS -N 100
#
# Example: make_spaghetti_plots.sh -PD -RWS -N 100


# endian independent script
TOBIG_ENDIANFILTER=cat
FROMBIG_ENDIANFILTER=cat

if uname -m | grep 86 > /dev/null; then
  TOBIG_ENDIANFILTER='suswapbytes format=1'
  FROMBIG_ENDIANFILTER='suswapbytes'
  echo 'little-endian detected'
fi

# trace to focus on
TR_NO=4

# priors for named fluid
delivery ModelDescription.xml -PD -m prior_traces.su -v 3 -o junk.su --trace-filter "tracl=$TR_NO" $*
# time range of desired spaghetti plot
TMIN=1.9
TMAX=2.3
(cat near-stack.su | $FROMBIG_ENDIANFILTER | suwind tmin=$TMIN tmax=$TMAX key=tracl min=$TR_NO max=$TR_NO ; \
deliveryAnalyser -i junk.su -s $TMIN $TMAX R_near wavelet.su | $FROMBIG_ENDIANFILTER | suwind tmin=$TMIN tmax=$TMAX ) | suxgraph reverse=1 linewidth=5,1 &
