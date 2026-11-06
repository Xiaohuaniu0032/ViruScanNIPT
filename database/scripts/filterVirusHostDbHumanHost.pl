#!/usr/bin/env perl
use strict;
use warnings;

# --------------------------------------------------
# Usage:
# perl filter_virushostdb_human_host.pl input.tsv output.tsv
# --------------------------------------------------

my ($input_file, $output_file) = @ARGV;

die "Usage: perl $0 input.tsv output.tsv\n"
    unless defined $input_file && defined $output_file;

open my $IN, '<', $input_file
    or die "Cannot open $input_file: $!\n";

open my $OUT, '>', $output_file
    or die "Cannot write $output_file: $!\n";

# Read header
my $header = <$IN>;
die "Input file is empty!\n" unless defined $header;

chomp $header;
$header =~ s/\r$//;

my @header = split /\t/, $header, -1;

# Identify column index
my %col;
for my $i (0 .. $#header) {
    $col{$header[$i]} = $i;
}

die "Missing required column: VirusHostDB_host_name\n"
    unless exists $col{'VirusHostDB_host_name'};

my $host_idx = $col{'VirusHostDB_host_name'};

# Output header
print $OUT "$header\n";

# Statistics
my $total_count  = 0;
my $matched_count = 0;

while (my $line = <$IN>) {

    chomp $line;
    $line =~ s/\r$//;

    next if $line eq '';

    $total_count++;

    my @fields = split /\t/, $line, -1;

    die "Invalid number of columns at line $.\n"
        if @fields <= $host_idx;

    # Retain Homo sapiens only
    next unless $fields[$host_idx] eq 'Homo sapiens';

    print $OUT "$line\n";
    $matched_count++;
}

close $IN;
close $OUT;

# Print statistics
print "Filtering completed!\n";
print "Total records:   $total_count\n";
print "Human records:   $matched_count\n";
print "Output file:     $output_file\n";
