#!/usr/bin/perl -w

use strict;

my @lines;
my $filename = $ARGV[0];
if (open(TEXT_FILE,$filename)) {
    my $line;
    my $i = 0;
    while(defined($line = <TEXT_FILE>)) {
	#print $line;
	$_ = $line;
	s/:::\{052\}/&/g;
	s/:::\{053\}/+/g;
	s/:::\{054\}/-/g;
	$lines[$i] = $_;
	$i++;
        #print $_;
	#print "<>\n";
    }    
    close(TEXT_FILE);
}

if (open(TEXT_FILE,">$filename")) {
    my $i = 0;
    for ($i = 0; $i < @lines; $i++) {
	my $line = $lines[$i];
	print TEXT_FILE $line;
    }
    close(TEXT_FILE);
}
