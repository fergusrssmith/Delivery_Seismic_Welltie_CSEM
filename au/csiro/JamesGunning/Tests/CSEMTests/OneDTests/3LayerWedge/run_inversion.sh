#!/usr/bin/env bash
#./clean.sh

j=0
for i in CMP* ; 
do
  dir="CMP"$j
  cd ${dir}

    csem1Dinverter -v 2 -P -q --LM -v 2 -s --MS r,2 --log-grid 30

    #csem1Dinverter -v 2 -P -q --LM -v 2 -s --MS d:r,10 --auto-split |grep Gridding:

    # model studies down to bottom depth bz
    # bz=`echo 500+$j*21|bc -l`
	# csem1Dinverter -v 2 --model-study 1,800,10,1200:2,10,20,$bz --LM -q --startup startup --NB
    cd ..
  ((j++))
done


