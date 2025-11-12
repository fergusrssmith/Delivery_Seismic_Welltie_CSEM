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
TR_NO=8

# make realisations
delivery ModelDescription.xml -PD -m prior_traces.su -v 3 -o junk.su --trace-filter "tracl=$TR_NO" $*
# time range of desired spaghetti plot
TMIN=1.7
TMAX=2.2
#deliveryAnalyser -i junk.su -S $TMIN $TMAX R_near wavelet.su near-stack.su --trace-filter "tracl=$TR_NO" | $FROMBIG_ENDIANFILTER | suxgraph reverse=1 linewidth=5,1 &

# PS you can also do this with 
deliveryAnalyser -i junk.su -S $TMIN $TMAX R_near wavelet.su near-stack.su --trace-filter "tracl=$TR_NO" -G 