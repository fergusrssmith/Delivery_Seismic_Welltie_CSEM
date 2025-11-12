#!/usr/bin/perl

# Script to make set of prior traces for a 3 layer sandwich with
$make_saturation_flat=0;
$S_oil_agnostic=0.7;
$S_true=0.9;

if ($#ARGV >= 0) 
{
	if ($ARGV[0] eq '-A')
	{
		$make_saturation_flat=1;
	} 
	else 
	{
		print " Usage: make_noncontact_model.pl [-A]\n where -A makes agnostic saturation model)\n";
		exit;
	}
}

# number of traces
$NT=21;   # 
# number of layers
$NL=3;    # 3
# no. times (incl. base)
$N_times = $NL + 1;


#make ascii data for layer times, properties etc
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";


$LE = 0.33*$NT;  # oil edge
$RE = 0.66*$NT;

for ($i=0;$i<$NT;$i++)
{
	$nsamples=0;
	# Blocks are
	# (1) saturations V1. Write in all slots so we can test delta_pore_pressure etc too
	for ($j=0;$j<$NL;$j++)
	{
		$s = ($j==1 && $i<$RE) ? $S_true : 0.0;
		if ($make_saturation_flat==1) {$s = (($j==1) ? $S_oil_agnostic: 0);}
		print FILE $s,"\n";
		$nsamples++;
	}
	# (2) saturations V2
	for ($j=0;$j<$NL;$j++)
	{
		$s = ($j==1 && $i<$LE) ? $S_true : 0.0;
		if ($make_saturation_flat==1) {$s = (($j==1) ? $S_oil_agnostic: 0);}
		print FILE $s,"\n";
		$nsamples++;
	}
	# (3) oil probability
	for ($j=0;$j<$NL;$j++)
	{
		$po = ($j!=1) ? 0 : (($i < $GE) ? 1 : 0);
		print FILE $po,"\n";
		$nsamples++;
	}
	# (4) times 
	$tstart=2.0;
	$tpad=0.05;
	# reservoir thickness in secs
	$DTR=0.03;
	for ($j=0;$j<=$NL;$j++)
	{
		$t1 = $tstart + 0.01*$i - ($tpad-$DTR);
		$t_l = $tstart + 0.01*$i+($DTR*$j);		
		# times
		#print FILE ($j==0) ? $t1 : (($j==$NL) ? 2.4 : $t_l),"\n";
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
}
print "removing: prior_data.txt suhead.su\n";
print "wrote prior_traces.su\n";
system("rm prior_data.txt suhead.su");
