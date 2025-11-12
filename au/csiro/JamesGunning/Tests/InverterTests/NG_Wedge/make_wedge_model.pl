#!/usr/bin/perl

# Script to make set of prior traces  with 3 layers, middle
# layer making a transition from NG=0 to NG=1.
# (must be big endian)
# number of traces
$NT=15;   # best at 15
# make blanks
system("sunull nt=4 ntr=$NT >| prior_traces.su");
# extract the headers
system("cat prior_traces.su | sustrip head=suhead.su > /dev/null");

#make data for layer times
open(FILE, ">prior_data.txt") or die "Can't open prior_data.txt";
for ($i=0;$i<$NT;$i++)
{
	#NG
	$ng=($i-2)/($NT-4);
	# True varying NG model
    $ng=(($ng<0) ? 0 : ($ng>1 ? 1 : $ng));
	print FILE "0\n",$ng,"\n0\n";
    print "0\n",$ng,"\n0\n";

    # tests for varying spatial trends, eg GGG factor
	# $ppty_of_interest=(($i<$NT/2) ? 300 : 250);
	# print FILE "$ppty_of_interest\n",$ppty_of_interest,"\n$ppty_of_interest\n";
	# print "$ppty_of_interest\n",$ppty_of_interest,"\n$ppty_of_interest\n";

    #times
	print FILE "1.8\n1.9\n1.95\n2.10\n";
	print "1.8\n1.9\n2.0\n2.10\n";
}
close FILE;
# glue these onto the headers
system("cat prior_data.txt | a2b n1=1 | supaste ns=7 head=suhead.su >| prior_traces.su");
# swap to big endian if os is little (detects intel linux box)
if (`uname -m` =~ m/86/)
{
	print "swapping bytes\n";
	system("cat prior_traces.su | suswapbytes format=1 >| tmp.su; mv -f tmp.su prior_traces.su");
}

print "Wrote prior_traces.su\n";
