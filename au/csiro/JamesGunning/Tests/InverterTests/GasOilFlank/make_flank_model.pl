#!/usr/bin/perl

# Script to make set of prior traces  with 3 layer times, with fluid transition B-O-G 
# - as if up a flank
# (must be big endian)
# number of traces
$NT=15;   # best at 15
# make blanks
system("sunull nt=9 ntr=$NT | sushw key=cdp a=1 b=1 >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");

#make data for layer times
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";
for ($i=0;$i<$NT;$i++)
{
	## oil and gas probabilities
	# Oil in middle third
	$p_oil = ($i>$NT/3 && $i<2*$NT/3) ? 1.0 : 0;
	# Gas in last third
	$p_gas = ($i>=2*$NT/3) ? 1.0 : 0;
	$vp_oil = 3000+2.0*$i;
    ## to make the fluid probs. vp_oil is just a test of other functionality (ignore)
	print FILE "0\n",$p_oil,"\n0\n",$p_gas,"\n",$vp_oil/2.0,"\n",$vp_oil,"\n";

	## to make the flank
	$t3 = 1.904-0.006*$i;
	print FILE "1.801\n",$t3,"\n2.101\n";
	print "trace ",($i+1),":\n p_oil=",$p_oil,", p_gas=",$p_gas," 1.801\n",$t3,"\n2.101\n\n";
}
close FILE;
# glue these onto the headers
system("cat prior_data.txt | a2b n1=1 | supaste ns=9 head=suhead.su >| prior_traces.su");
# swap to big endian if os is little (detects intel linux box)
if (`uname -m` =~ m/86/)
{
	print "swapping bytes\n";
	system("cat prior_traces.su | suswapbytes format=1 >| tmp.su; mv -f tmp.su prior_traces.su");
}
print "wrote prior_traces.su\n";

