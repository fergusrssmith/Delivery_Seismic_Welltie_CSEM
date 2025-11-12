use strict;
#one line println's only!!!!
my @lines;
my $index = 0;
my $filename = $ARGV[0];
if (open(FILE,$filename)) {
    my $line;
    while (defined($line = <FILE>)) {
	$_ = $line;
	if (/\/\/System\.out\.print/ || /\/\/System\.err\.print/) {
	    #print "comment:   $_";
	} else {
	    $lines[$index] = $_;
	    $index++;
	    #print "no comment:$_";
	}
    }
    close(FILE);
    if (open(WRITE,">$filename")) {
	foreach my $line (@lines) {
	    print WRITE $line;    
	}
    }
    close(WRITE);
}
