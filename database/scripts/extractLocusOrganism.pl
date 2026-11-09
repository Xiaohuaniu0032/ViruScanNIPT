#!/usr/bin/env perl
use strict;
use warnings;

my ($input) = @ARGV;
die "Usage: perl $0 viral.1.genomic.gbff[.gz]\n"
    unless defined $input && @ARGV == 1;
-f $input && -r $input
    or die "[ERROR] Cannot read input file: $input\n";

my $fh;
if ($input =~ /\.gz\z/i) {
    open $fh, '-|', 'gzip', '-dc', '--', $input
        or die "[ERROR] Cannot start gzip for $input: $!\n";
} else {
    open $fh, '<', $input
        or die "[ERROR] Cannot open $input: $!\n";
}

my ($locus, $organism) = ('', '');
my $records_written = 0;

while (my $line = <$fh>) {
    $line =~ s/\r?\n\z//;

    if ($line =~ /^LOCUS\s+(\S+)/) {
        $locus = $1;
    }
    elsif ($line =~ /^  ORGANISM\s+(.+)/) {
        $organism = $1;
    }
    elsif ($line =~ m{^//\s*$}) {
        if (length $locus && length $organism) {
            print $locus, "\t", $organism, "\n";
            ++$records_written;
        }
        ($locus, $organism) = ('', '');
    }
}

close $fh
    or die "[ERROR] Input read/decompression failed for $input: $!\n";

if (length $locus && length $organism) {
    print $locus, "\t", $organism, "\n";
    ++$records_written;
}

warn "[INFO] Records written: $records_written\n";
