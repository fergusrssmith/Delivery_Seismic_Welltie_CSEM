#!/usr/bin/env bash
# A test for traces with large nt values
nt=5000

TOBIG_ENDIANFILTER=cat
FROMBIG_ENDIANFILTER=cat

if uname -m | grep 86 > /dev/null; then
  TOBIG_ENDIANFILTER='suswapbytes format=1'
  FROMBIG_ENDIANFILTER='suswapbytes'
  echo 'little-endian detected'
fi

sunull nt=$nt ntr=100 | $TOBIG_ENDIANFILTER >| bigtraces.su
cat bigtraces.su | $FROMBIG_ENDIANFILTER | sustrip head=suhead.su > /dev/null
perl -e 'for ($i=0;$i<5000*100;$i++) {print rand(),"\n";}' | a2b n1=1 | supaste ns=$nt head=suhead.su  | $TOBIG_ENDIANFILTER >| bigtraces.su
echo checking with SuRange
cat bigtraces.su | java au.csiro.JamesGunning.sulib.SuRange
echo wrote bigtraces.su

