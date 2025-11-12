#!/usr/bin/perl

# Script to make set of prior traces for a 3 layer sandwich with
# fluid contact in the middle

# number of traces
$NT=15;   # best at 15
# number of layers
$NL=3;    # 3
# no. times (incl. base)
$N_times = $NL + 1;

# depths arrays for surfaces as function of offset
@d0 = (1900,1964.5,2024.7014);
#@d_slope = (0,12.915,12.915);
@d_slope = (12.915,12.915,12.915);

#make data for layer times
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";

for ($i=0;$i<$NT;$i++)
{
	$nsamples=0;
#depths 
	for ($j=0;$j<$NL;$j++)
	{	
		$z = $d0[$j] + $i*$d_slope[$j];
		print FILE $z,"\n"; 
		$nsamples++;	
	}
# times 
	$tstart=2.0;
	$tpad=0.05;
    # reservoir thickness in secs
	$DTR=0.05;
	for ($j=0;$j<=$NL;$j++)
	{
		#$t_l = $t1 + 0.01*$i+($DTR*$j);
		$t1 = $tstart + 0.01*$i - ($tpad-$DTR);
		$t_l = $tstart + 0.01*$i+($DTR*$j);	
		#print FILE ($j==0) ? $t1 : (($j==$NL) ? (2.4+.01*$i) : $t_l),"\n";
		print FILE ($j==0) ? $t1 : (($j==$NL) ? ($t_l+$tpad-$DTR) : $t_l),"\n";	
		$nsamples++;
	}
}
close FILE;

# make blanks
system("sunull nt=$nsamples ntr=$NT >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");
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
