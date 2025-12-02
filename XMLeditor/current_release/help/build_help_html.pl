#!/usr/bin/perl -w
use strict;

my @lines;
my $index = 0;
#my $filename = $ARGV[0];
if (open(TEXT_FILE,"generate_html.txt")) {
    my $line;
    while (defined($line = <TEXT_FILE>)) {
	my $filename = $line;
	chomp($filename) if defined($filename);
	my $title = <TEXT_FILE>;
	chomp($title) if defined($title);
	my $text = <TEXT_FILE>;
	chomp($text) if defined($text);
	my $params = "";
	my $space = <TEXT_FILE>;
	chomp($space) if defined($space);
	while (defined($space) && $space ne "") {
	    $params .= $space;
	    $params .= "</P>\n        <P>";
	    $space = <TEXT_FILE>;
	    chomp($space);
	}
	if (open(WRITE,">$filename") && open(TEMPLATE_FILE,"template.html")) {
	    while (defined($line = <TEMPLATE_FILE>)) {
		$_ = $line;
		s/[<]!--title--[>]/$title/;
		s/[<]!--text--[>]/$text/;
		s/[<]!--parameters--[>]/$params/;
		print WRITE $_;
	    }
	    close(WRITE);
	    print "wrote file $filename\n";
	}
    }
}
close(TEMPLATE_FILE);
close(TEXT_FILE);
