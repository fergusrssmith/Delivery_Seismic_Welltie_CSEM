#!/usr/bin/perl

# Script to make set of prior traces  with 3 thickness std deviations + 4 layer times, one making a pinchout
# (must be big endian)
$make_flat=0;
$force_PO=0;

if ($ARGV[0] eq '-f') {
    $make_flat=1;
} elsif ($ARGV[0] eq '-w') {
    $make_flat=0;
} else {
	print " Usage: make_wedge_model [-f,-w] [-force]\n where -f is flat model, -w is a wedge, adding -force forces pinchout\n";
	exit;
}

if (($#ARGV > 0) && ($ARGV[1] eq '-force')) {
    $force_PO=1;
}

# number of traces
$NT=15;   # best at 15
# make blanks
system("sunull nt=7 ntr=$NT >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");

$dz=1000;   # thickness default std dev (wide & uninformative when applied)

#make data for layer times
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";
$t2=1.9112;
for ($i=0;$i<$NT;$i++)
{
	$t3 = 1.904+0.006*$i; 
	$t3 = ($t3 < $t2 )? $t2 : $t3;
# to introduce bias, add +/-0.005ms if desired....
	$t3 += 0.00;
	if ($make_flat==1) {
		## to make a flat structure
		print FILE $dz,"\n",($force_PO==1 && ($t3==$t2))?0:$dz,"\n",$dz,"\n1.801\n",$t2,"\n",1.95,"\n2.101\n";
		print $dz,"\n",($force_PO==1 && ($t3==$t2))?0:$dz,"\n",$dz,"\n1.801\n",$t2,"\n",1.95,"\n2.101\n\n";
	} else {
		## to make a wedge
		print FILE $dz,"\n",($force_PO==1 && ($t3==$t2))?0:$dz,"\n",$dz,"\n1.801\n",$t2,"\n",$t3,"\n2.101\n";
		print "trace ",($i+1),":\n",$dz,"\n",($force_PO==1 && ($t3==$t2))?0:$dz,"\n",$dz,"\n1.801\n",$t2,"\n",$t3,"\n2.101\n\n";
	}
}
close FILE;
# glue these onto the headers
system("cat prior_data.txt | a2b n1=1 | supaste ns=7 head=suhead.su >| prior_traces.su");
# swap to big endian if os is little (detects intel linux box)
if (`uname -m` =~ m/86/)
{
	print "swapping bytes\n";
	system("cat prior_traces.su | suswapbytes format=1 >| tmp.su; mv -f tmp.su prior_traces.su");
	print "output written to prior_traces\n";
}

