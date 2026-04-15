#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
use Data::Dumper qw(Dumper);
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_05_parse_blast.pl
#
# Description:
#   Parse BLAST output for one sample, filter valid hits by:
#     1. percent identity
#     2. alignment length
#     3. e-value
#
#   Then map accession/locus ID to viral organism name using
#   locus_organism.tsv, and generate:
#     1. raw blast result table
#     2. parsed blast result table
#
# Output files:
#   1. <sample>_blast_raw_result.xls
#   2. <sample>_blast_out_parsed.txt
#==================================================

my %opt;
GetOptions(
    "sample=s"                  => \$opt{sample},
    "blast-out=s"               => \$opt{blast_out},
    "outdir=s"                  => \$opt{outdir},
    "locus-organism=s"          => \$opt{locus_organism},
    "identity-cutoff=f"         => \$opt{identity_cutoff},
    "alignment-length-cutoff=i" => \$opt{alignment_length_cutoff},
    "evalue-cutoff=s"           => \$opt{evalue_cutoff},
    "help"                      => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/sample blast_out outdir locus_organism identity_cutoff alignment_length_cutoff evalue_cutoff/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{blast_out}      = abs_path($opt{blast_out});
$opt{outdir}         = normalize_path($opt{outdir});
$opt{locus_organism} = abs_path($opt{locus_organism});

-f $opt{blast_out}      or die "[ERROR] blast output not found: $opt{blast_out}\n";
-f $opt{locus_organism} or die "[ERROR] locus_organism file not found: $opt{locus_organism}\n";
mkdir_if_not_exists($opt{outdir});

my $sample = $opt{sample};

my $blast_raw = "$opt{outdir}/${sample}_blast_raw_result.xls";
my $outfile   = "$opt{outdir}/${sample}_blast_out_parsed.txt";

#--------------------------------------------------
# load accession -> organism mapping
#--------------------------------------------------
my %acc2name;
open my $fh_map, "<", $opt{locus_organism} or die "[ERROR] Cannot open $opt{locus_organism}\n";
while (<$fh_map>) {
    chomp;
    next if /^\s*$/;
    next if /^#/;

    my @arr = split /\t/, $_;
    next unless @arr >= 2;

    my $acc  = $arr[0];
    my $name = $arr[1];

    $acc =~ s/\.\d+$//;  # remove version
    $acc2name{$acc} = $name;
}
close $fh_map;

#--------------------------------------------------
# parse blast output
#--------------------------------------------------
open my $fh_raw, ">", $blast_raw or die "[ERROR] Cannot write $blast_raw\n";
print $fh_raw join("\t",
    qw(query subject identity alignment_length mismatches gap_opens q.start q.end s.start s.end evalue bit_score)
), "\n";

open my $fh_out, ">", $outfile or die "[ERROR] Cannot write $outfile\n";
print $fh_out join("\t",
    qw(sample query subject identity alignment_length evalue viral_name)
), "\n";

open my $fh_in, "<", $opt{blast_out} or die "[ERROR] Cannot open $opt{blast_out}\n";

my %hits;

while (<$fh_in>) {
    chomp;
    next if /^\s*$/;
    next if /^#/;

    my @val = split /\t/, $_;
    next unless @val >= 12;

    my $query            = $val[0];
    my $subject_raw      = $val[1];
    my $identity         = $val[2];
    my $alignment_length = $val[3];
    my $mismatches       = $val[4];
    my $gap_opens        = $val[5];
    my $q_start          = $val[6];
    my $q_end            = $val[7];
    my $s_start          = $val[8];
    my $s_end            = $val[9];
    my $evalue           = $val[10];
    my $bit_score        = $val[11];

    # 只保留常见 RefSeq accession
    next unless $subject_raw =~ /^(NC|AC|NG)_/;

    print $fh_raw join("\t",
        $query, $subject_raw, $identity, $alignment_length, $mismatches,
        $gap_opens, $q_start, $q_end, $s_start, $s_end, $evalue, $bit_score
    ), "\n";

    my $subject = $subject_raw;
    $subject =~ s/\.\d+$//;   # remove version

    if (
        $identity >= $opt{identity_cutoff}
        && $alignment_length >= $opt{alignment_length_cutoff}
        && blast_evalue_le($evalue, $opt{evalue_cutoff})
    ) {
        push @{$hits{$query}}, join("\t", $subject, $identity, $alignment_length, $evalue);
    }
}
close $fh_in;
close $fh_raw;

#--------------------------------------------------
# keep one representative hit per query
# current logic follows your earlier script:
#   if >=2 hits exist, keep the first one
#   if only 1 hit exists, keep that one
#--------------------------------------------------
my %uniq_hit;
for my $query (keys %hits) {
    my $hit_num = scalar @{$hits{$query}};
    if ($hit_num >= 1) {
        $uniq_hit{$query} = $hits{$query}->[0];
    }
}

for my $query (sort keys %uniq_hit) {
    my @info = split /\t/, $uniq_hit{$query};
    my ($hit, $identity, $alignment_length, $evalue) = @info;

    my $viral_name = exists $acc2name{$hit} ? $acc2name{$hit} : "NA";

    print $fh_out join("\t",
        $sample,
        $query,
        $hit,
        $identity,
        $alignment_length,
        $evalue,
        $viral_name
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
  perl VSN_05_parse_blast.pl \
    --sample sample_id \
    --blast-out sample.blast.out \
    --outdir sample_dir \
    --locus-organism locus_organism.tsv \
    --identity-cutoff 97 \
    --alignment-length-cutoff 35 \
    --evalue-cutoff 1e-5

Description:
  Parse one sample BLAST output and generate:
    1. <sample>_blast_raw_result.xls
    2. <sample>_blast_out_parsed.txt

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

sub blast_evalue_le {
    my ($a, $b) = @_;
    # force numeric comparison, works for strings like 1e-5, 2.3e-10, 0.00001
    return ($a + 0) <= ($b + 0);
}
