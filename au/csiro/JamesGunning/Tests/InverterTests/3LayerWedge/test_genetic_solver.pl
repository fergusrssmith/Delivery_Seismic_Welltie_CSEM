#!/usr/bin/perl

# Script to check robustness of genetic solver
# How reproducable is the mode finding under random variation of seed points???
#
# We test trace 14, where the modes are well separated. Fool around as you please...
#
system("rm junk.txt");
for ($i=0;$i<20;$i++)
{
# loop over mode fnding, with variable seeds to Genetic solver
	system("delivery ModelDescription.xml -m prior_traces.su -v 4 -PD -o junk.su -N 1 -ID --IF 0 --seed $i --trace-filter 'tracl=1' |            awk '/Best_seismic_chisq/{print \$3}' >> junk.txt");
}
system("graph -TX --frame-color wheat --bg-color black -a junk.txt");
