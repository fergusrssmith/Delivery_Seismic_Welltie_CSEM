#!/usr/bin/perl

# Script to make set of prior traces for a 3 layer sandwich with
$make_saturation_flat=0;
$S_oil_agnostic=0.6;

if ($#ARGV >= 0) 
{
	if ($ARGV[0] eq '-A')
	{
		$make_saturation_flat=1;
	} 
	else 
	{
		print " Usage: make_stress_model.pl [-A]\n where -A makes agnostic saturation model)\n";
		exit;
	}
}


# number of traces
$NT=21;   # best at 15
# number of layers
$NL=3;    # 3

#make ascii data for layer times, properties etc
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";

print "Making model, agnostic_saturation = ",$make_saturation_flat,"\n";

for ($i=0;$i<$NT;$i++)
{
	$nsamples=0;
	# Blocks are
	# (1) pore-pressure ("pseudo-normalised")
	for ($j=0;$j<$NL;$j++)
	{
		$s = 2*($i-0.5*($NT-1))/($NT-1);
		print FILE $s,"\n";
		$nsamples++;
	}
	# (2) oil saturation V1
	for ($j=0;$j<$NL;$j++)
	{
		$so = ($j!=1) ? 0 : (($i < 0.667*$NT) ? 0.7 : 0);
		if ($make_saturation_flat==1) {$so = (($j!=1) ? 0 : $S_oil_agnostic);}
		print FILE $so,"\n";
		$nsamples++;
	}
    # (3) oil saturation V2
	for ($j=0;$j<$NL;$j++)
	{
		$so = ($j!=1) ? 0 : (($i < 0.333*$NT) ? 0.7 : ((($i < 0.667*$NT) ? 0.1: 0)));
		if ($make_saturation_flat==1) {$so = (($j!=1) ? 0 : $S_oil_agnostic);}
		print FILE $so,"\n";
		$nsamples++;
	}
	$t1=2.0;
   # (4) depths
	for ($j=0;$j<$NL;$j++)
	{
		$t_l = $t1 + 0.01*$i+($DTR*$j);
		$z = 1900+0.5*(0.01*$i)*2583.0;
		print FILE $z,"\n";
		$nsamples++;
	}
	# (5) times 

	# reservoir thickness in secs
	$DTR=0.03;
	for ($j=0;$j<=$NL;$j++)
	{
		$t_l = $t1 + 0.01*$i+($DTR*$j);
		# times
		print FILE  (($j==$NL) ? $t_l+0.03: (($j==0) ? ($t_l-0.03) : $t_l)),"\n";
		$nsamples++;
	}
}
close FILE;
print "nsamples=",$nsamples,"\n";

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
}
print "removing: prior_data.txt suhead.su\n";
print "wrote prior_traces.su\n";
system("rm prior_data.txt suhead.su");
