#!/usr/bin/env perl
use strict;
use warnings;

# --------------------------------------------------
# Usage:
# perl filter_viral_genomic_dna_host.pl input.tsv output.tsv
# --------------------------------------------------

my ($input_file, $output_file) = @ARGV;

die "Usage: perl $0 input.tsv output.tsv\n"
    unless defined $input_file && defined $output_file;

open my $IN,  '<', $input_file  or die "Cannot open $input_file: $!\n";
open my $OUT, '>', $output_file or die "Cannot write $output_file: $!\n";

# Read header
my $header = <$IN>;
die "Input file is empty!\n" unless defined $header;

chomp $header;
$header =~ s/\r$//;

my @header = split /\t/, $header, -1;

# Identify column indexes
my %col;
for my $i (0 .. $#header) {
    $col{$header[$i]} = $i;
}

die "Missing required column: mol_type\n"
    unless exists $col{'mol_type'};

my $mol_idx = $col{'mol_type'};

# Output header
print $OUT "$header\n";

# Statistics
my $total_count = 0;
my $dna_count   = 0;

while (my $line = <$IN>) {

    chomp $line;
    $line =~ s/\r$//;

    next if $line eq '';

    $total_count++;

    my @fields = split /\t/, $line, -1;

    die "Invalid number of columns at line $. \n"
        if @fields <= $mol_idx;

    # Retain genomic DNA
    next unless $fields[$mol_idx] eq 'genomic DNA';

    print $OUT "$line\n";

    $dna_count++;
}

close $IN;
close $OUT;

# Print statistics
print "Filtering completed!\n";
print "Total records:          $total_count\n";
print "Genomic DNA records:    $dna_count\n";
print "Output file: $output_file\n";
