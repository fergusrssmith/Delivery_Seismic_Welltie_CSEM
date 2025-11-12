#!/usr/bin/env bash

TOBIG_ENDIANFILTER=cat
FROMBIG_ENDIANFILTER=cat

if uname -m | grep 86 > /dev/null; then
  TOBIG_ENDIANFILTER='suswapbytes format=1'
  FROMBIG_ENDIANFILTER='suswapbytes'
fi

export nlayers=$*
cat  realisations.su | $FROMBIG_ENDIANFILTER | suwind dt=10 itmin=0 itmax=`echo $nlayers-1|bc` | sustrip | b2a n1=1 \
| awk -v nd=$nlayers '{if ((NR-1)%nd==0){print "\n",$0}else{print $0}}'
