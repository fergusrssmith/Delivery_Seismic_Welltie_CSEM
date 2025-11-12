#!/usr/bin/perl

# Script to make set of prior traces for a 3 layer sandwich with
# fluid contact in the middle

# number of traces
$NT=15;   # best at 15
# number of layers
$NL=3;    # 3
# no. times (incl. base)
$N_times = $NL + 1;
# make blanks
system("sunull nt=$N_times ntr=$NT >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");

$dz=1000;   # thickness default std dev (wide & uninformative when applied)

#make data for layer times
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";
for ($i=0;$i<$NT;$i++)
{
	$nsamples=0;
# for testing variable contact detail
# for ($j=0;$j<$NL;$j++)
# 	{
# 		print FILE ($j==0) ? 2050+$i : (($j==1) ? 2120+$i : 0),"\n";
# 		$nsamples++;
# 	}

# for testing variable depths
	$d1=1900;
	for ($j=0;$j<$NL;$j++)
	{
#		$d_l = $d1 + 5*$i+(50*$j);
		$d_l = $d1;
		print FILE $d_l,"\n";
		$nsamples++;
	}

	# times 
	$t1=2.0;
	for ($j=0;$j<=$NL;$j++)
	{
		$t_l = $t1 + 0.01*$i+(0.05*$j);
		print FILE ($j==0) ? $t1 : (($j==$NL) ? 2.4 : $t_l),"\n";
#		print FILE ($j==0) ? ($t1+0.01*$i) : (($j==$NL) ? (2.2 + 0.01*$i) : $t_l),"\n";
		$nsamples++;
	}
}
close FILE;
# glue these onto the headers
system("cat prior_data.txt | a2b n1=1 | supaste ns=$nsamples head=suhead.su >| prior_traces.su");
# swap to big endian if os is little (detects intel linux box)
if (`uname -m` =~ m/86/)
{
	print "swapping bytes\n";
	system("cat prior_traces.su | suswapbytes format=1 >| tmp.su; mv -f tmp.su prior_traces.su");
	print "output written to prior_traces.su:\n";	
}
system("rm prior_data.txt suhead.su");
