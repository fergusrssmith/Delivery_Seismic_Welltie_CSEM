#!/usr/bin/perl

# Script to make set of prior traces like Mike Glinsk's standard wedge
# test model, with 4 layer times, one making a pinchout
# (must be big endian)
# number of traces
$NT=200;   # 200 traces
# make blanks
system("sunull nt=4 ntr=$NT >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");

#make data for layer times
# 400 ft drop across 200 traces, plug speed 7134 ft/s
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";
$t2=1.9112;
for ($i=0;$i<$NT;$i++)
{
	$t3 = $t2+2*0.000281756*$i;
	$t3 = ($t3 < $t2 )? $t2 : $t3;
	## to make a wedge
	print FILE "1.801\n",$t2,"\n",$t3,"\n2.101\n";
	print "trace ",($i+1),":\n1.801\n",$t2,"\n",$t3,"\n2.101\n\n";
	## to make a flat structure
    #print FILE "1.801\n",$t2,"\n",1.95,"\n2.101\n";
	#print "1.801\n",$t2,"\n",1.95,"\n2.101\n\n";
}
close FILE;
# glue these onto the headers
system("cat prior_data.txt | a2b n1=1 | supaste ns=4 head=suhead.su >| prior_traces.su");
# swap to big endian if os is little (detects intel linux box)
if (`uname -m` =~ m/86/)
{
	print "swapping bytes\n";
	system("cat prior_traces.su | suswapbytes format=1 >| tmp.su; mv -f tmp.su prior_traces.su");
}

