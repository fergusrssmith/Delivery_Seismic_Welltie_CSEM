#!/usr/bin/env bash

# Assemble profile_final.MS.txt files in each CMP directory by parsing 
# ModelStudy.txt files for optimal model, as ranked by MML
# Very hacky script: watch out if ModelStudy.txt format changes: assumes
# last columns in  ModelStudy.txt are log10rho_i

j=0
for i in CMP* ; 
do
  dir="CMP"$j
  cd $dir
  awk 'BEGIN{min=1e10} \
       /# m/{for(i=1;i<=NF;i++){if($i=="d_0"){ds=i-1;}; if($i=="log10rho_0"){rs=i-1;}; if($i=="logMMLapprox"){ms=i-1;}}} \
       NR>1 && $ms<min{min=$ms;BL=$0} \
       END {j=ds;$0=BL;for(i=rs;i<=NF;i++){printf("%g %g\n%g %g\n",$j,$i,$(j+1),$i);j++};}' ModelStudy.txt \
       >| Profile_final.MS.txt
  cd ..
  ((j++))
done



