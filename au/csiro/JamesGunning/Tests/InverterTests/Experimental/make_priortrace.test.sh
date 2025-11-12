#!/usr/bin/env bash
nt=`cat prior_traces.test.txt|awk 'END{print NR}'`
echo $nt entries

TOBIG_ENDIANFILTER=cat
FROMBIG_ENDIANFILTER=cat

if uname -m | grep 86 > /dev/null; then
  TOBIG_ENDIANFILTER='suswapbytes format=1'
  FROMBIG_ENDIANFILTER='suswapbytes'
  echo 'little-endian detected'
fi

sunull nt=$nt ntr=1 | $TOBIG_ENDIANFILTER >| prior_traces.su
cat prior_traces.su | $FROMBIG_ENDIANFILTER | sustrip head=suhead.su > /dev/null
cat prior_traces.test.txt | a2b n1=1 | supaste ns=$nt head=suhead.su  | $TOBIG_ENDIANFILTER >| prior_traces.su
echo check...
cat prior_traces.su | $FROMBIG_ENDIANFILTER | surange
cat prior_traces.su | $FROMBIG_ENDIANFILTER | sustrip | b2a n1=1
echo "wrote prior_traces.su"


