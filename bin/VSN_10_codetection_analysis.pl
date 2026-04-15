#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_10_codetection_analysis.pl
#
# Description:
#   Perform co-detection analysis to detect pairs or higher combinations of viruses
#   present in the same sample.
#
#   Main functions:
#     1. read Virus_Group_Statistics_Final.xls
#     2. detect co-detected virus groups in each sample
#     3. generate detailed co-detection file and summary file
#     4. rank the co-detection combinations by their frequency
#
# Output:
#   Sample_CoDetection_Detail.xls
#   Virus_CoDetection_Summary.xls
#==================================================

my %opt;
GetOptions(
    "input=s"    => \$opt{input},
    "outdir=s"   => \$opt{outdir},
    "help"       => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/input outdir/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{input}  = abs_path($opt{input});
$opt{outdir} = normalize_path($opt{outdir});

-f $opt{input} or die "[ERROR] input file not found: $opt{input}\n";
mkdir_if_not_exists($opt{outdir});

my $detail_file = "$opt{outdir}/Sample_CoDetection_Detail.xls";
my $summary_file = "$opt{outdir}/Virus_CoDetection_Summary.xls";

#--------------------------------------------------
# storage for co-detection combinations
#--------------------------------------------------
my %sample_pos_groups;  # $sample_pos_groups{sample} = [group1, group2, ...]
my %combo_freq;         # $combo_freq{combo} = frequency

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

    next if !defined $status || $status eq "NA" || !defined $group || $group eq "NA";

    # Track positive samples with their virus groups
    if (defined $status && $status eq 'Positive') {
        push @{$sample_pos_groups{$sample}}, $group;
    }
}
close $fh_in;

#--------------------------------------------------
# generate detailed co-detection results
#--------------------------------------------------
open my $fh_det, ">", $detail_file or die "[ERROR] Cannot write $detail_file\n";
print $fh_det "SampleID\tPositive_Group_Count\tCoDetection_Combination\n";

foreach my $s (sort keys %sample_pos_groups) {
    my @groups = sort @{$sample_pos_groups{$s}};
    my $count  = scalar @groups;

    if ($count >= 2) {
        my $combo = join(" + ", @groups);
        print $fh_det "$s\t$count\t$combo\n";
        $combo_freq{$combo}++;
    }
}
close $fh_det;

#--------------------------------------------------
# generate summary file with co-detection frequency ranking
#--------------------------------------------------
open my $fh_sum, ">", $summary_file or die "[ERROR] Cannot write $summary_file\n";
print $fh_sum "Virus_Combination\tDetected_Sample_Count\n";

# Sort combinations by frequency (descending)
my @sorted_combos = sort { $combo_freq{$b} <=> $combo_freq{$a} } keys %combo_freq;

foreach my $combo (@sorted_combos) {
    my $freq = $combo_freq{$combo};
    print $fh_sum "$combo\t$freq\n";
}
close $fh_sum;

exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl VSN_10_codetection_analysis.pl \
    --input Virus_Group_Statistics_Final.xls \
    --outdir 03_summary_analysis/05_codetection

Description:
  Perform co-detection analysis for all samples, identifying virus combinations
  and counting their occurrences.

Output:
  Sample_CoDetection_Detail.xls
  Virus_CoDetection_Summary.xls

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

