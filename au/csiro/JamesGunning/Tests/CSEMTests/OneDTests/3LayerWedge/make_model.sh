
#!/usr/bin/env bash
./clean.sh

freqs=".25,.5,.75,1,1.25,1.5,2"
offsets="1000,500,15000"


thicknesses=`seq -s " " 10 25 450`
echo "Thicknesses: " ${thicknesses}

#main reservoir thickness
log10res=2



rm GeomCMP.txt
j=0
for i in ${thicknesses}; 
do
  # accumulate GeomCMP.txt file, make directories for inversions
  dir="CMP"$j
  echo ${dir} " " $j " 0 0" >> GeomCMP.txt
  if [ ! -d ${dir} ]; then mkdir ${dir} ; fi
  cd ${dir} 
    # clean 
    rm -f *
    # make a data set template
    java au.csiro.JamesGunning.CSEM.OneDinversion.DataSetMaker -f ${freqs} --off ${offsets}
    # make models & startup files	
	last_thick=`echo 1500-$i|bc`
    java au.csiro.JamesGunning.CSEM.OneDinversion.DataSetMaker --model 1000,0.32:1000,0:$i,${log10res}:${last_thick},0 -s
    # make data
    csem1Dinverter -q -P --make-synthetic 5 2e-16 -v 1
    mv synthetic_data.dat data.dat
	rm iter_final
  cd ..	
  ((j++))
done


