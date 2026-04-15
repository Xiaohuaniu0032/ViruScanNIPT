#!/usr/bin/env perl
use strict;
use warnings;

my ($in) = @ARGV;
die "Usage: perl $0 viral.1.genomic.gbff[.gz]\n" unless $in;

my $fh;
if ($in =~ /\.gz$/) {
    open($fh, "-|", "gzip -dc $in") or die "Cannot open $in: $!\n";
} else {
    open($fh, "<", $in) or die "Cannot open $in: $!\n";
}

my ($locus, $org) = ("", "");

while (my $line = <$fh>) {
    chomp $line;

    if ($line =~ /^LOCUS\s+(\S+)/) {
        $locus = $1;
    }
    elsif ($line =~ /^  ORGANISM\s+(.+)/) {
        $org = $1;
    }
    elsif ($line =~ m{^//}) {
        if ($locus ne "") {
            print $locus, "\t", $org, "\n";
        }
        ($locus, $org) = ("", "");
    }
}

close $fh;

