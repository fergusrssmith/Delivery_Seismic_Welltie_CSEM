#!/usr/bin/env bash

# clean up. Preserve ModelStudy.* files

for i in CMP* ; 
do
  if [ -d $i ]
  then 
	cd $i
    rm -f CorrCoeff.* Covariance.* SummaryStats.* Mode[^l]*  resp_final.* iter_final.* \
		iter_final_* profile_final.* profile_final_*
	cd ..
  fi	
done
