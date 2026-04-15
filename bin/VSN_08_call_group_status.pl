#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_08_call_group_status.pl
#
# Description:
#   Summarize Groupname hit counts for each sample and assign
#   Positive/Negative status.
#
#   Main functions:
#     1. read Annotated_Human_Virus_Results.xls
#     2. count hits by sample and Groupname
#     3. assign status using positive-hit-cutoff
#     4. if one sample has no valid Groupname, output one NA row
#
# Output:
#   Virus_Group_Statistics_Final.xls
#
# Output columns:
#   sample
#   Groupname
#   Groupname_Count
#   Status
#==================================================

my %opt;
GetOptions(
    "input=s"                => \$opt{input},
    "outdir=s"               => \$opt{outdir},
    "positive-hit-cutoff=i"  => \$opt{positive_hit_cutoff},
    "help"                   => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/input outdir positive_hit_cutoff/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{input}  = abs_path($opt{input});
$opt{outdir} = normalize_path($opt{outdir});

-f $opt{input} or die "[ERROR] input file not found: $opt{input}\n";
mkdir_if_not_exists($opt{outdir});

my $outfile = "$opt{outdir}/Virus_Group_Statistics_Final.xls";

#--------------------------------------------------
# storage
#--------------------------------------------------
my %stats;         # $stats{sample}{group} = count
my @sample_order;  # preserve sample order
my %seen_sample;

#--------------------------------------------------
# read input
# expected columns:
#   sample query subject identity alignment_length evalue viral_name Groupname
#--------------------------------------------------
open my $fh_in, "<", $opt{input} or die "[ERROR] Cannot open $opt{input}\n";
my $header = <$fh_in>;
defined $header or die "[ERROR] input file is empty: $opt{input}\n";

while (<$fh_in>) {
    chomp;
    next if /^\s*$/;

    my @fields = split /\t/, $_;
    next unless @fields >= 8;

    my $sample = $fields[0];
    my $query  = $fields[1];
    my $group  = $fields[-1];

    if (!$seen_sample{$sample}) {
        push @sample_order, $sample;
        $seen_sample{$sample} = 1;
    }

    next if !defined $query || $query eq 'NA';
    next if !defined $group || $group eq 'NA';

    $stats{$sample}{$group}++;
}
close $fh_in;

#--------------------------------------------------
# output
#--------------------------------------------------
open my $fh_out, ">", $outfile or die "[ERROR] Cannot write $outfile\n";
print $fh_out join("\t", qw(sample Groupname Groupname_Count Status)), "\n";

for my $sample (@sample_order) {
    if (exists $stats{$sample} && %{$stats{$sample}}) {
        for my $group (sort keys %{$stats{$sample}}) {
            my $count  = $stats{$sample}{$group};
            my $status = ($count >= $opt{positive_hit_cutoff}) ? "Positive" : "Negative";

            print $fh_out join("\t", $sample, $group, $count, $status), "\n";
        }
    } else {
        print $fh_out join("\t", $sample, "NA", "NA", "NA"), "\n";
    }
}
close $fh_out;

exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl VSN_08_call_group_status.pl \
    --input Annotated_Human_Virus_Results.xls \
    --outdir 03_summary_analysis/03_group_status \
    --positive-hit-cutoff 2

Description:
  Count Groupname hits for each sample and assign Positive/Negative
  status using the specified cutoff.

Output:
  Virus_Group_Statistics_Final.xls

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

