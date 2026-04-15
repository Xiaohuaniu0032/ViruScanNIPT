#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   VSN_01_make_sample_shells.pl
#
# Description:
#   Read fq.list and generate one executable shell script (*.sh)
#   for each sample.
#
#   Each generated shell script includes:
#     1. bwa mem to hg38
#     2. SAM -> BAM
#     3. BAM sort/index
#     4. extract unmapped reads
#     5. unmapped BAM -> FASTA
#     6. BLAST against viral database
#     7. parse BLAST output
#
#   This script does NOT execute the shell scripts.
#==================================================

my %opt;
GetOptions(
    "fq-list=s"                 => \$opt{fq_list},
    "outdir=s"                  => \$opt{outdir},
    "bin-dir=s"                 => \$opt{bin_dir},
    "bwa=s"                     => \$opt{bwa},
    "samtools=s"                => \$opt{samtools},
    "blastn=s"                  => \$opt{blastn},
    "ref=s"                     => \$opt{ref},
    "viral-db=s"                => \$opt{viral_db},
    "locus-organism=s"          => \$opt{locus_organism},
    "bwa-threads=i"             => \$opt{bwa_threads},
    "blast-threads=i"           => \$opt{blast_threads},
    "max-target-seqs=i"         => \$opt{max_target_seqs},
    "max-hsps=i"                => \$opt{max_hsps},
    "identity-cutoff=f"         => \$opt{identity_cutoff},
    "alignment-length-cutoff=i" => \$opt{alignment_length_cutoff},
    "evalue-cutoff=s"           => \$opt{evalue_cutoff},
    "keep-intermediate=i"       => \$opt{keep_intermediate},
    "help"                      => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (
    qw/
      fq_list
      outdir
      bin_dir
      bwa
      samtools
      blastn
      ref
      viral_db
      locus_organism
      bwa_threads
      blast_threads
      max_target_seqs
      max_hsps
      identity_cutoff
      alignment_length_cutoff
      evalue_cutoff
      keep_intermediate
    /
) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{fq_list} = abs_path($opt{fq_list});
$opt{outdir}  = normalize_path($opt{outdir});
$opt{bin_dir} = normalize_path($opt{bin_dir});

-f $opt{fq_list} or die "[ERROR] fq.list not found: $opt{fq_list}\n";
-d $opt{bin_dir} or die "[ERROR] bin-dir not found: $opt{bin_dir}\n";

my $parse_script = "$opt{bin_dir}/VSN_05_parse_blast.pl";

mkdir_if_not_exists($opt{outdir});

open my $fh, "<", $opt{fq_list} or die "[ERROR] Cannot open fq.list: $opt{fq_list}\n";

while (<$fh>) {
    chomp;
    s/^\s+//;
    s/\s+$//;
    next if $_ eq '';
    next if /^#/;

    my $fq = $_;
    my $sample = get_sample_name($fq);

    my $sample_dir   = "$opt{outdir}/$sample";
    my $sh_file      = "$sample_dir/$sample.run_viral_scan.sh";
    my $parse_log    = "$sample_dir/$sample.parse_blast.log";

    my $sam          = "$sample_dir/$sample.sam";
    my $bam          = "$sample_dir/$sample.bam";
    my $sort_bam     = "$sample_dir/$sample.sort.bam";
    my $unmapped_bam = "$sample_dir/$sample.unmapped.bam";
    my $unmapped_fa  = "$sample_dir/$sample.unmapped.fa";
    my $blast_out    = "$sample_dir/$sample.blast.out";

    mkdir_if_not_exists($sample_dir);

    open my $sh, ">", $sh_file or die "[ERROR] Cannot write shell file: $sh_file\n";

    print $sh "#!/bin/bash\n";
    print $sh "set -euo pipefail\n\n";
    print $sh "# sample: $sample\n";
    print $sh "# fastq : $fq\n\n";

    my $cmd;

    # 1. bwa mem
    $cmd = qq{$opt{bwa} mem -t $opt{bwa_threads} $opt{ref} $fq > $sam};
    print $sh "$cmd\n";

    # 2. sam -> bam
    $cmd = qq{$opt{samtools} view -b -o $bam $sam};
    print $sh "$cmd\n";

    # 3. sort bam
    $cmd = qq{$opt{samtools} sort -o $sort_bam $bam};
    print $sh "$cmd\n";

    # 4. index bam
    $cmd = qq{$opt{samtools} index $sort_bam};
    print $sh "$cmd\n";

    # 5. extract unmapped reads
    $cmd = qq{$opt{samtools} view -b -f 4 -o $unmapped_bam $sort_bam};
    print $sh "$cmd\n";

    # 6. bam -> fasta
    $cmd = qq{$opt{samtools} fasta $unmapped_bam > $unmapped_fa};
    print $sh "$cmd\n";

    # 7. blastn
    $cmd = qq{$opt{blastn} -db $opt{viral_db} -query $unmapped_fa -outfmt 7 -max_target_seqs $opt{max_target_seqs} -max_hsps $opt{max_hsps} -num_threads $opt{blast_threads} > $blast_out};
    print $sh "$cmd\n";

    # 8. parse blast
    if (-f $parse_script) {
        $cmd = join(" ",
            "perl",
            shell_quote($parse_script),
            "--sample", shell_quote($sample),
            "--blast-out", shell_quote($blast_out),
            "--outdir", shell_quote($sample_dir),
            "--locus-organism", shell_quote($opt{locus_organism}),
            "--identity-cutoff", shell_quote($opt{identity_cutoff}),
            "--alignment-length-cutoff", shell_quote($opt{alignment_length_cutoff}),
            "--evalue-cutoff", shell_quote($opt{evalue_cutoff}),
            ">", shell_quote($parse_log),
            "2>&1"
        );
        print $sh "$cmd\n";
    } else {
        print $sh qq{echo "[ERROR] required script not found: $parse_script" >&2\n};
        print $sh "exit 1\n";
    }

    # remove selected intermediate files
    if ($opt{keep_intermediate} == 0) {
        print $sh "rm -f $sam\n";
        print $sh "rm -f $bam\n";
        print $sh "rm -f $sort_bam\n";
        print $sh "rm -f $sort_bam.bai\n";
    }

    close $sh;

    chmod 0755, $sh_file or die "[ERROR] Cannot chmod 755 for $sh_file\n";

    print "Generated shell for sample: $sample\n";
}
close $fh;

exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl VSN_01_make_sample_shells.pl \
    --fq-list fq.list \
    --outdir 01_sample_shells \
    --bin-dir /path/to/bin \
    --bwa /path/to/bwa \
    --samtools /path/to/samtools \
    --blastn /path/to/blastn \
    --ref /path/to/hg38.fa \
    --viral-db /path/to/viral \
    --locus-organism /path/to/locus_organism.tsv \
    --bwa-threads 4 \
    --blast-threads 4 \
    --max-target-seqs 5 \
    --max-hsps 1 \
    --identity-cutoff 97 \
    --alignment-length-cutoff 35 \
    --evalue-cutoff 1e-5 \
    --keep-intermediate 0

Description:
  Read fq.list, create one directory per sample, and generate one executable
  shell script (*.sh) for each sample. The shell script is NOT executed here.
  The generated shell includes BLAST parsing at the end.

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

sub get_sample_name {
    my ($fq) = @_;
    my $base = basename($fq);

    $base =~ s/\.fastq\.gz$//;
    $base =~ s/\.fq\.gz$//;
    $base =~ s/\.fastq$//;
    $base =~ s/\.fq$//;

    return $base;
}

sub shell_quote {
    my ($str) = @_;
    $str =~ s/'/'"'"'/g;
    return "'$str'";
}

