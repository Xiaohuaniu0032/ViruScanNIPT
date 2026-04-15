#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path getcwd);
use List::Util qw(sum);

#==================================================
# Script Name:
#   VSN_11_group_distribution_stats.pl
#
# Description:
#   Summarize count distribution statistics for positive Groupname records.
#
#   Main functions:
#     1. read Virus_Group_Statistics_Final.xls
#     2. keep only Status = Positive rows
#     3. output sample-level detail table for each Groupname
#     4. calculate numerical statistics for each Groupname:
#          - Sample_Size
#          - Min
#          - Max
#          - Mean
#          - Median
#          - Q1_25%
#          - Q3_75%
#
# Output:
#   1. Groupname_Sample_Count_Detail.xls
#   2. Groupname_Numerical_Statistics.xls
#==================================================

my %opt;
GetOptions(
    "input=s"  => \$opt{input},
    "outdir=s" => \$opt{outdir},
    "help"     => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/input outdir/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{input}  = abs_path($opt{input});
$opt{outdir} = normalize_path($opt{outdir});

-f $opt{input} or die "[ERROR] input file not found: $opt{input}\n";
mkdir_if_not_exists($opt{outdir});

my $detail_file = "$opt{outdir}/Groupname_Sample_Count_Detail.xls";
my $stats_file  = "$opt{outdir}/Groupname_Numerical_Statistics.xls";

#--------------------------------------------------
# storage
#--------------------------------------------------
my %group_counts;   # $group_counts{Groupname} = [count1, count2, ...]
my %group_samples;  # $group_samples{Groupname} = [ [sample, count], ... ]

#--------------------------------------------------
# read input
# expected columns:
#   sample Groupname Groupname_Count Status
#--------------------------------------------------
open my $fh_in, "<", $opt{input} or die "[ERROR] Cannot open $opt{input}\n";
my $header = <$fh_in>;
defined $header or die "[ERROR] input file is empty: $opt{input}\n";

while (<$fh_in>) {
    chomp;
    next if /^\s*$/;

    my ($sample, $group, $count, $status) = split /\t/, $_;

    next unless defined $status && $status eq 'Positive';
    next unless defined $group  && $group ne 'NA';
    next unless defined $count  && $count =~ /^\d+$/;

    push @{$group_counts{$group}}, $count;
    push @{$group_samples{$group}}, [$sample, $count];
}
close $fh_in;

#--------------------------------------------------
# sort groups by sample size descending
#--------------------------------------------------
my @sorted_groups = sort {
       scalar(@{$group_counts{$b}}) <=> scalar(@{$group_counts{$a}})
    || $a cmp $b
} keys %group_counts;

#--------------------------------------------------
# output detail file
#--------------------------------------------------
open my $fh_det, ">", $detail_file or die "[ERROR] Cannot write $detail_file\n";
print $fh_det join("\t", qw(Groupname Sample Count)), "\n";

for my $group (@sorted_groups) {
    my @samples = sort {
           $b->[1] <=> $a->[1]
        || $a->[0] cmp $b->[0]
    } @{$group_samples{$group}};

    for my $pair (@samples) {
        print $fh_det join("\t", $group, $pair->[0], $pair->[1]), "\n";
    }
}
close $fh_det;

#--------------------------------------------------
# output stats file
#--------------------------------------------------
open my $fh_stat, ">", $stats_file or die "[ERROR] Cannot write $stats_file\n";
print $fh_stat join("\t", qw(Groupname Sample_Size Min Max Mean Median Q1_25% Q3_75%)), "\n";

for my $group (@sorted_groups) {
    my @vals = sort { $a <=> $b } @{$group_counts{$group}};
    my $n = scalar @vals;

    my $min_val = $vals[0];
    my $max_val = $vals[-1];
    my $mean    = sum(@vals) / $n;
    my $median  = calculate_quantile(\@vals, 0.50);
    my $q1      = calculate_quantile(\@vals, 0.25);
    my $q3      = calculate_quantile(\@vals, 0.75);

    printf $fh_stat "%s\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\n",
        $group, $n, $min_val, $max_val, $mean, $median, $q1, $q3;
}
close $fh_stat;

exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl VSN_11_group_distribution_stats.pl \
    --input Virus_Group_Statistics_Final.xls \
    --outdir 03_summary_analysis/06_group_distribution

Description:
  Summarize count distribution statistics for positive Groupname records.

Outputs:
  1. Groupname_Sample_Count_Detail.xls
  2. Groupname_Numerical_Statistics.xls

USAGE
    exit(1);
}

sub normalize_path {
    my ($path) = @_;
    if ($path =~ m{^/}) {
        return $path;
    } else {
        return abs_path(getcwd()) . "/$path";
    }
}

sub mkdir_if_not_exists {
    my ($dir) = @_;
    return if -d $dir;
    mkdir $dir or die "[ERROR] cannot create directory: $dir\n";
}

sub calculate_quantile {
    my ($data_ref, $quantile) = @_;
    my @data = @$data_ref;
    my $n = scalar @data;

    return undef if $n == 0;

    my $pos    = $quantile * ($n - 1);
    my $lower  = int($pos);
    my $upper  = $lower + 1;
    my $weight = $pos - $lower;

    if ($upper >= $n) {
        return $data[$lower];
    } else {
        return $data[$lower] * (1 - $weight) + $data[$upper] * $weight;
    }
}

