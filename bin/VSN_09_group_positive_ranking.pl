#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_09_group_positive_ranking.pl
#
# Description:
#   Summarize positive sample counts for each Groupname and
#   generate a ranked output table.
#
#   Main functions:
#     1. read Virus_Group_Statistics_Final.xls
#     2. count total sample number
#     3. count total positive sample number
#     4. count positive sample number for each Groupname
#     5. rank Groupname by PositiveSampleCount (descending)
#
# Output:
#   Groupname_Positive_Ranking.xls
#
# Output columns:
#   Groupname
#   PositiveSampleCount
#   TotalSampleCount
#   Positive_Percentage
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

my $outfile = "$opt{outdir}/Groupname_Positive_Ranking.xls";

#--------------------------------------------------
# storage
#--------------------------------------------------
my %group_positive_samples;   # $group_positive_samples{group}{sample} = 1
my %total_samples;            # all sample IDs
my %positive_samples;         # all positive sample IDs

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

    next unless defined $sample;
    $total_samples{$sample} = 1;

    if (defined $status && $status eq 'Positive') {
        $positive_samples{$sample} = 1;
        $group_positive_samples{$group}{$sample} = 1 if defined $group && $group ne 'NA';
    }
}
close $fh_in;

#--------------------------------------------------
# overall stats
#--------------------------------------------------
my $total_count    = scalar keys %total_samples;
my $positive_count = scalar keys %positive_samples;
my $overall_rate   = $total_count > 0 ? ($positive_count / $total_count * 100) : 0;

print STDERR "=" x 30, "\n";
print STDERR "Overall positive summary:\n";
print STDERR "Total samples: $total_count\n";
print STDERR "Positive samples: $positive_count\n";
printf STDERR "Overall positive rate: %.2f%%\n", $overall_rate;
print STDERR "=" x 30, "\n";

#--------------------------------------------------
# prepare ranking results
#--------------------------------------------------
my @results;
for my $group (keys %group_positive_samples) {
    my $pos_sample_count = scalar keys %{$group_positive_samples{$group}};
    my $percentage = $total_count > 0 ? ($pos_sample_count / $total_count * 100) : 0;

    push @results, {
        name     => $group,
        count    => $pos_sample_count,
        raw_perc => $percentage,
        perc     => sprintf("%.2f%%", $percentage),
    };
}

@results = sort {
       $b->{count}    <=> $a->{count}
    || $a->{name} cmp $b->{name}
} @results;

#--------------------------------------------------
# write output
#--------------------------------------------------
open my $fh_out, ">", $outfile or die "[ERROR] Cannot write $outfile\n";
print $fh_out join("\t", qw(Groupname PositiveSampleCount TotalSampleCount Positive_Percentage)), "\n";

for my $res (@results) {
    print $fh_out join("\t",
        $res->{name},
        $res->{count},
        $total_count,
        $res->{perc}
    ), "\n";
}
close $fh_out;

exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl VSN_09_group_positive_ranking.pl \
    --input Virus_Group_Statistics_Final.xls \
    --outdir 03_summary_analysis/04_group_ranking

Description:
  Summarize positive sample counts for each Groupname and generate:
    Groupname_Positive_Ranking.xls

Output columns:
  Groupname
  PositiveSampleCount
  TotalSampleCount
  Positive_Percentage

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

