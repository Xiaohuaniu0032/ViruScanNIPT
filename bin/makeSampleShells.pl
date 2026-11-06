#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# Generate one sample directory and one executable shell directly under --outdir.
# No per-step subdirectories are created.

my %opt;
GetOptions(
    "fq-list=s"                  => \$opt{fq_list},
    "outdir=s"                   => \$opt{outdir},
    "bin-dir=s"                  => \$opt{bin_dir},
    "bwa=s"                      => \$opt{bwa},
    "samtools=s"                 => \$opt{samtools},
    "blastn=s"                   => \$opt{blastn},
    "ref=s"                      => \$opt{ref},
    "viral-db=s"                 => \$opt{viral_db},
    "locus-organism=s"           => \$opt{locus_organism},
    "bwa-threads=i"              => \$opt{bwa_threads},
    "blast-threads=i"            => \$opt{blast_threads},
    "max-target-seqs=i"          => \$opt{max_target_seqs},
    "max-hsps=i"                 => \$opt{max_hsps},
    "identity-cutoff=f"          => \$opt{identity_cutoff},
    "alignment-length-cutoff=i"  => \$opt{alignment_length_cutoff},
    "evalue-cutoff=s"            => \$opt{evalue_cutoff},
    "keep-intermediate=i"        => \$opt{keep_intermediate},
    "help"                       => \$opt{help},
) or usage(1);

usage(0) if $opt{help};

for my $key (qw(
    fq_list outdir bin_dir bwa samtools blastn ref viral_db locus_organism
    bwa_threads blast_threads max_target_seqs max_hsps identity_cutoff
    alignment_length_cutoff evalue_cutoff keep_intermediate
)) {
    die "[ERROR] --$key is required\n" unless defined $opt{$key};
}

die "[ERROR] --keep-intermediate must be 0 or 1\n"
    unless $opt{keep_intermediate} == 0 || $opt{keep_intermediate} == 1;

$opt{fq_list} = abs_path($opt{fq_list})
    or die "[ERROR] Cannot resolve fq.list: $opt{fq_list}\n";
$opt{outdir} = normalize_path($opt{outdir});
$opt{bin_dir} = normalize_path($opt{bin_dir});
$opt{ref} = abs_path($opt{ref})
    or die "[ERROR] Cannot resolve reference: $opt{ref}\n";
$opt{locus_organism} = abs_path($opt{locus_organism})
    or die "[ERROR] Cannot resolve locus-organism table: $opt{locus_organism}\n";

-f $opt{fq_list} or die "[ERROR] fq.list not found: $opt{fq_list}\n";
-d $opt{bin_dir} or die "[ERROR] bin directory not found: $opt{bin_dir}\n";
-f $opt{ref} or die "[ERROR] Reference file not found: $opt{ref}\n";
-f $opt{locus_organism} or die "[ERROR] locus-organism table not found: $opt{locus_organism}\n";
-f "$opt{bin_dir}/parseBlast.pl"
    or die "[ERROR] Required parser not found: $opt{bin_dir}/parseBlast.pl\n";

make_path($opt{outdir}) unless -d $opt{outdir};

my $fq_list_dir = dirname($opt{fq_list});
my %seen_sample;
my $generated = 0;

open my $LIST, '<', $opt{fq_list}
    or die "[ERROR] Cannot open fq.list: $opt{fq_list}: $!\n";

while (my $line = <$LIST>) {
    chomp $line;
    $line =~ s/\r$//;
    $line =~ s/^\s+|\s+$//g;
    next if $line eq '' || $line =~ /^#/;

    my $fq = resolve_fastq_path($line, $fq_list_dir);
    my $sample = get_sample_name($fq);

    die "[ERROR] Invalid sample name derived from FASTQ: $fq\n"
        unless $sample =~ /\A[A-Za-z0-9_.-]+\z/ && $sample ne '.' && $sample ne '..';
    die "[ERROR] Duplicate sample name '$sample' in fq.list\n"
        if $seen_sample{$sample}++;

    my $sample_dir = "$opt{outdir}/$sample";
    make_path($sample_dir) unless -d $sample_dir;

    my $shell_file = "$sample_dir/$sample.run.sh";
    my $parse_script = "$opt{bin_dir}/parseBlast.pl";

    my $sam          = "$sample_dir/$sample.sam";
    my $bam          = "$sample_dir/$sample.bam";
    my $sort_bam     = "$sample_dir/$sample.sort.bam";
    my $unmapped_bam = "$sample_dir/$sample.unmapped.bam";
    my $unmapped_fa  = "$sample_dir/$sample.unmapped.fa";
    my $blast_out    = "$sample_dir/$sample.blast.out";

    my $bwa_log      = "$sample_dir/$sample.bwa.log";
    my $samtools_log = "$sample_dir/$sample.samtools.log";
    my $blast_log    = "$sample_dir/$sample.blast.log";
    my $parse_log    = "$sample_dir/$sample.parseBlast.log";

    open my $SH, '>', $shell_file
        or die "[ERROR] Cannot write shell file $shell_file: $!\n";

    print {$SH} "#!/usr/bin/env bash\n";
    print {$SH} "set -euo pipefail\n\n";
    print {$SH} "# ViruScanNIPT sample-level workflow\n\n";

    # 1. Align reads to the human reference.
    print {$SH} qdouble($opt{bwa}), " mem -t ", qdouble($opt{bwa_threads}),
        " ", qdouble($opt{ref}), " ", qdouble($fq),
        " > ", qdouble($sam), " 2> ", qdouble($bwa_log), "\n\n";

    # 2. Convert, sort, and index BAM.
    print {$SH} qdouble($opt{samtools}), " view -b -o ",
        qdouble($bam), " ", qdouble($sam),
        " 2> ", qdouble($samtools_log), "\n";
    print {$SH} qdouble($opt{samtools}), " sort -o ",
        qdouble($sort_bam), " ", qdouble($bam),
        " 2>> ", qdouble($samtools_log), "\n";
    print {$SH} qdouble($opt{samtools}), " index ",
        qdouble($sort_bam), " 2>> ", qdouble($samtools_log), "\n\n";

    # 3. Extract unmapped reads.
    print {$SH} qdouble($opt{samtools}), " view -b -f 4 -o ",
        qdouble($unmapped_bam), " ", qdouble($sort_bam),
        " 2>> ", qdouble($samtools_log), "\n";
    print {$SH} qdouble($opt{samtools}), " fasta ",
        qdouble($unmapped_bam), " > ", qdouble($unmapped_fa),
        " 2>> ", qdouble($samtools_log), "\n\n";

    # 4. BLASTN against the viral reference database.
    print {$SH} qdouble($opt{blastn}),
        " -db ", qdouble($opt{viral_db}),
        " -query ", qdouble($unmapped_fa),
        " -outfmt 7",
        " -max_target_seqs ", qdouble($opt{max_target_seqs}),
        " -max_hsps ", qdouble($opt{max_hsps}),
        " -num_threads ", qdouble($opt{blast_threads}),
        " -out ", qdouble($blast_out),
        " 2> ", qdouble($blast_log), "\n\n";

    # 5. Parse BLAST output.
    print {$SH} "perl ", qdouble($parse_script),
        " --sample ", qdouble($sample),
        " --blast-out ", qdouble($blast_out),
        " --outdir ", qdouble($sample_dir),
        " --locus-organism ", qdouble($opt{locus_organism}),
        " --identity-cutoff ", qdouble($opt{identity_cutoff}),
        " --alignment-length-cutoff ", qdouble($opt{alignment_length_cutoff}),
        " --evalue-cutoff ", qdouble($opt{evalue_cutoff}),
        " > ", qdouble($parse_log), " 2>&1\n\n";

    # 6. Remove selected intermediate files if requested.
    if ($opt{keep_intermediate} == 0) {
        print {$SH} "rm -f ",
            join(" ", map { qdouble($_) } ($sam, $bam, $sort_bam, "$sort_bam.bai")),
            "\n";
    }

    close $SH or die "[ERROR] Cannot finish shell file $shell_file: $!\n";
    chmod 0755, $shell_file
        or die "[ERROR] Cannot chmod shell file $shell_file: $!\n";

    $generated++;
    print "[INFO] Generated: $shell_file\n";
}
close $LIST;

die "[ERROR] No valid FASTQ entries found in $opt{fq_list}\n" unless $generated;

print "[INFO] Samples generated: $generated\n";
print "[INFO] Output directory : $opt{outdir}\n";
exit 0;

sub resolve_fastq_path {
    my ($path, $list_dir) = @_;
    return abs_path($path) if -f $path;
    my $relative_to_list = "$list_dir/$path";
    return abs_path($relative_to_list) if -f $relative_to_list;
    die "[ERROR] FASTQ file not found: $path\n";
}

sub get_sample_name {
    my ($fq) = @_;
    my $base = basename($fq);
    $base =~ s/\.(?:fastq|fq)(?:\.gz)?$//i;
    return $base;
}

sub normalize_path {
    my ($path) = @_;
    return abs_path($path) if -d $path;
    if ($path =~ m{^/}) {
        return $path;
    }
    return getcwd() . "/$path";
}

# Bash double-quoted argument. Single quotes are deliberately prohibited
# to guarantee that generated sample shells contain no single-quote character.
sub qdouble {
    my ($value) = @_;
    $value = "$value";
    die "[ERROR] A shell argument contains a single quote, which is not allowed: $value\n"
        if index($value, "'") >= 0;
    die "[ERROR] A shell argument contains a newline: $value\n"
        if $value =~ /[\r\n]/;
    $value =~ s/\\/\\\\/g;
    $value =~ s/"/\\"/g;
    $value =~ s/\$/\\\$/g;
    $value =~ s/`/\\`/g;
    return '"' . $value . '"';
}

sub usage {
    my ($exit_code) = @_;
    print STDERR <<'USAGE';
Usage:
  perl bin/makeSampleShells.pl \
    --fq-list fq.list \
    --outdir ViruScanNIPT_result \
    --bin-dir bin \
    --bwa /path/to/bwa \
    --samtools /path/to/samtools \
    --blastn /path/to/blastn \
    --ref /path/to/hg38.fa \
    --viral-db /path/to/viral \
    --locus-organism database/annotation/locusOrganism.tsv \
    --bwa-threads 4 \
    --blast-threads 4 \
    --max-target-seqs 5 \
    --max-hsps 1 \
    --identity-cutoff 97 \
    --alignment-length-cutoff 35 \
    --evalue-cutoff 1e-5 \
    --keep-intermediate 0

Creates:
  ViruScanNIPT_result/<sample>/<sample>.run.sh

No additional per-step directories are created. The generated shell is
not executed by this script.
USAGE
    exit($exit_code);
}
