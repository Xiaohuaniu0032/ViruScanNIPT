#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename dirname);
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_06_merge_sample_results.pl
#
# Description:
#   Merge parsed BLAST results from all sample directories.
#
#   Input directory structure:
#     sample_dir/
#       ├── sampleA/
#       │   └── sampleA_blast_out_parsed.txt
#       ├── sampleB/
#       │   └── sampleB_blast_out_parsed.txt
#       └── ...
#
#   Main functions:
#     1. detect all sample subdirectories
#     2. find *_blast_out_parsed.txt for each sample
#     3. filter out records with viral_name containing "phage"
#     4. merge all valid records into one summary table
#     5. if a sample has no valid record, output one NA row
#
# Output:
#   merge_blast_parsed_res.xls
#==================================================

my %opt;
GetOptions(
    "sample-dir=s" => \$opt{sample_dir},
    "outdir=s"     => \$opt{outdir},
    "help"         => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/sample_dir outdir/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{sample_dir} = normalize_path($opt{sample_dir});
$opt{outdir}     = normalize_path($opt{outdir});

-d $opt{sample_dir} or die "[ERROR] sample-dir not found: $opt{sample_dir}\n";
mkdir_if_not_exists($opt{outdir});

my $outfile = "$opt{outdir}/merge_blast_parsed_res.xls";

open my $fh_out, ">", $outfile or die "[ERROR] Cannot write $outfile\n";
print $fh_out join("\t", qw(sample query subject identity alignment_length evalue viral_name)), "\n";

my @sample_dirs = sort grep { -d $_ } glob("$opt{sample_dir}/*");

for my $sdir (@sample_dirs) {
    my $sample = basename($sdir);

    my @parsed_files = glob("$sdir/*_blast_out_parsed.txt");

    if (!@parsed_files) {
        print $fh_out join("\t", $sample, ("NA") x 6), "\n";
        next;
    }

    my $file = $parsed_files[0];
    my $eff_line = 0;

    open my $fh_in, "<", $file or die "[ERROR] Cannot open $file\n";

    my $header = <$fh_in>;
    defined $header or do {
        close $fh_in;
        print $fh_out join("\t", $sample, ("NA") x 6), "\n";
        next;
    };

    while (<$fh_in>) {
        chomp;
        next if /^\s*$/;

        my @arr = split /\t/, $_;
        next unless @arr >= 7;

        my $viral_name = $arr[-1];

        if ($viral_name !~ /phage/i) {
            print $fh_out "$_\n";
            $eff_line++;
        }
    }
    close $fh_in;

    if ($eff_line == 0) {
        print $fh_out join("\t", $sample, ("NA") x 6), "\n";
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
  perl VSN_06_merge_sample_results.pl \
    --sample-dir 01_sample_shells \
    --outdir 03_summary_analysis/01_merge

Description:
  Merge all *_blast_out_parsed.txt files under sample directories.
  Records with viral_name containing "phage" will be excluded.
  If a sample has no valid record, one NA row will be added.

Output:
  merge_blast_parsed_res.xls

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

