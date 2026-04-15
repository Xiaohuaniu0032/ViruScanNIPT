#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_07_annotate_human_dna_virus.pl
#
# Description:
#   Annotate merged BLAST parsed results using a human-host DNA virus database.
#
#   Main functions:
#     1. read merged result file
#     2. read human-host DNA virus database
#     3. remove accession version suffix before matching
#     4. append Groupname column
#     5. if one sample has no matched human-host DNA virus record,
#        output one NA row for that sample
#
# Input:
#   1. merge_blast_parsed_res.xls
#   2. human DNA virus database (accession<TAB>Groupname)
#
# Output:
#   Annotated_Human_Virus_Results.xls
#==================================================

my %opt;
GetOptions(
    "input=s"    => \$opt{input},
    "outdir=s"   => \$opt{outdir},
    "virus-db=s" => \$opt{virus_db},
    "help"       => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/input outdir virus_db/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{input}    = abs_path($opt{input});
$opt{outdir}   = normalize_path($opt{outdir});
$opt{virus_db} = abs_path($opt{virus_db});

-f $opt{input}    or die "[ERROR] input file not found: $opt{input}\n";
-f $opt{virus_db} or die "[ERROR] virus database file not found: $opt{virus_db}\n";
mkdir_if_not_exists($opt{outdir});

my $outfile = "$opt{outdir}/Annotated_Human_Virus_Results.xls";

#--------------------------------------------------
# load virus database
# format:
#   accession<TAB>Groupname
#--------------------------------------------------
my %db;
open my $fh_db, "<", $opt{virus_db} or die "[ERROR] Cannot open $opt{virus_db}\n";
while (<$fh_db>) {
    chomp;
    next if /^\s*$/;
    next if /^#/;

    my @arr = split /\t/, $_;
    next unless @arr >= 2;

    my $accession = $arr[0];
    my $groupname = $arr[1];

    $accession =~ s/\.\d+$//;   # remove version suffix
    $db{$accession} = $groupname;
}
close $fh_db;

#--------------------------------------------------
# read merged input and annotate
#--------------------------------------------------
open my $fh_in,  "<", $opt{input}  or die "[ERROR] Cannot open $opt{input}\n";
open my $fh_out, ">", $outfile     or die "[ERROR] Cannot write $outfile\n";

my $header = <$fh_in>;
defined $header or die "[ERROR] input file is empty: $opt{input}\n";
chomp $header;
print $fh_out "$header\tGroupname\n";

my %sample_map;
my @sample_order;
my %seen_sample;

while (<$fh_in>) {
    chomp;
    next if /^\s*$/;

    my @fields = split /\t/, $_;
    next unless @fields >= 7;

    my $sample  = $fields[0];
    my $subject = $fields[2];

    if (!$seen_sample{$sample}) {
        push @sample_order, $sample;
        $seen_sample{$sample} = 1;
        $sample_map{$sample} = [];
    }

    # NA line from previous step
    if ($subject eq 'NA') {
        next;
    }

    my $subject_clean = $subject;
    $subject_clean =~ s/\.\d+$//;

    if (exists $db{$subject_clean}) {
        push @{$sample_map{$sample}}, join("\t", @fields) . "\t" . $db{$subject_clean};
    }
}
close $fh_in;

#--------------------------------------------------
# output by original sample order
# if one sample has no matched records, output one NA row
# header now has 8 columns:
#   sample query subject identity alignment_length evalue viral_name Groupname
#--------------------------------------------------
for my $sample (@sample_order) {
    if (@{$sample_map{$sample}}) {
        for my $line (@{$sample_map{$sample}}) {
            print $fh_out "$line\n";
        }
    } else {
        print $fh_out join("\t", $sample, ("NA") x 7), "\n";
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
  perl VSN_07_annotate_human_dna_virus.pl \
    --input merge_blast_parsed_res.xls \
    --outdir 03_summary_analysis/02_human_dna_virus_annotation \
    --virus-db 224_human_host_DNA_viruses.txt

Description:
  Annotate merged sample results using a human-host DNA virus database
  and generate:
    Annotated_Human_Virus_Results.xls

Input virus-db format:
  accession<TAB>Groupname

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
