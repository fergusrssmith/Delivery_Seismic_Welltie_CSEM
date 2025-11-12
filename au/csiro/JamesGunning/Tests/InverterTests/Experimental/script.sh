#!/bin/bash
rm -f junk.txt
j=0
while  [ $j -le 12 ] 
do
	echo "run number " $j
   java au.csiro.JamesGunning.seismicinverter.SeismicInverter ModelDescription.xml RockandFluidProperties.xml SeismicModel.xml -v 4 -m prior_traces.su -I| grep 'starting times\|Best' >> junk.txt
    j=$[$j+1]
done
awk '/Fluid/{gsub(/.*times/,"");print $0}' junk.txt | extract.pl -2 1,2,3,4,5 | gr -a -h 0.7 -L starts
awk '/Best/{gsub(/.*times/,"");print $0}' junk.txt | extract.pl -2 1,2,3,4,5 | gr -a -h 0.7 -L solutions
awk '/Best/{gsub(/,/,"");print $6}' junk.txt | gr -a -h 0.7 -L chisquare
