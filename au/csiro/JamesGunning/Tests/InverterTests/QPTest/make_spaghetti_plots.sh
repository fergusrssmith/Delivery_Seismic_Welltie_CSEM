#!/usr/bin/env bash
#
# requires some or all of the all the usual delivery arguments -II -IS -N 100 --IF etc

# endian independent script
TOBIG_ENDIANFILTER=cat
FROMBIG_ENDIANFILTER=cat

if uname -m | grep 86 > /dev/null; then
  TOBIG_ENDIANFILTER='suswapbytes format=1'
  FROMBIG_ENDIANFILTER='suswapbytes'
  echo 'little-endian detected'
fi


# priors for named fluid
delivery ModelDescription.xml -v 4 -PD -RWS -N 50 -o junk.su $*
# time range of desired spaghetti plot
TMIN=3.75
TMAX=4.25
(cat near-stack.su | $FROMBIG_ENDIANFILTER | suwind tmin=$TMIN tmax=$TMAX ; \
deliveryAnalyser -i junk.su -s $TMIN $TMAX R_near wavelet.su | $FROMBIG_ENDIANFILTER | suwind tmin=$TMIN tmax=$TMAX ) | suxgraph reverse=1 linewidth=3,1 &
