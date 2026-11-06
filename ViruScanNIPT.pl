#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ViruScanNIPT.pl is a framework generator.
# It creates one directory per FASTQ sample directly under --outdir
# and writes one executable <sample>.run.sh in each sample directory.
# Sample-level shells are not executed automatically.
# Cohort-level statistics are handled separately by:
#   bin/summarizeGroupnameResults.pl

my %opt;
GetOptions(
    "config=s"  => \$opt{config},
    "fq-list=s" => \$opt{fq_list},
    "outdir=s"  => \$opt{outdir},
    "help"      => \$opt{help},
) or usage(1);

usage(0) if $opt{help};

for my $key (qw(config fq_list outdir)) {
    die "[ERROR] --$key is required\n" unless defined $opt{$key};
}

$opt{config} = abs_path($opt{config})
    or die "[ERROR] Cannot resolve config file: $opt{config}\n";
$opt{fq_list} = abs_path($opt{fq_list})
    or die "[ERROR] Cannot resolve fq.list: $opt{fq_list}\n";
$opt{outdir} = normalize_path($opt{outdir});

-f $opt{config} or die "[ERROR] Config file not found: $opt{config}\n";
-f $opt{fq_list} or die "[ERROR] FASTQ list not found: $opt{fq_list}\n";

my $project_dir = abs_path(dirname($0));
my $bin_dir = "$project_dir/bin";
my $shell_generator = "$bin_dir/makeSampleShells.pl";

-f $shell_generator
    or die "[ERROR] Sample shell generator not found: $shell_generator\n";

my %conf = read_config($opt{config});
validate_config(\%conf);
check_required_paths(\%conf);

make_path($opt{outdir}) unless -d $opt{outdir};

my @cmd = (
    $^X,
    $shell_generator,
    "--fq-list", $opt{fq_list},
    "--outdir", $opt{outdir},
    "--bin-dir", $bin_dir,
    "--bwa", $conf{BWA},
    "--samtools", $conf{SAMTOOLS},
    "--blastn", $conf{BLASTN},
    "--ref", $conf{HG38_REF},
    "--viral-db", $conf{VIRAL_DB},
    "--locus-organism", $conf{LOCUS_ORGANISM},
    "--bwa-threads", $conf{BWA_THREADS},
    "--blast-threads", $conf{BLAST_THREADS},
    "--max-target-seqs", $conf{MAX_TARGET_SEQS},
    "--max-hsps", $conf{MAX_HSPS},
    "--identity-cutoff", $conf{IDENTITY_CUTOFF},
    "--alignment-length-cutoff", $conf{ALIGNMENT_LENGTH_CUTOFF},
    "--evalue-cutoff", $conf{EVALUE_CUTOFF},
    "--keep-intermediate", $conf{KEEP_INTERMEDIATE},
);

print "[INFO] Generating sample shells\n";
print "[INFO] FASTQ list      : $opt{fq_list}\n";
print "[INFO] Output directory: $opt{outdir}\n";

my $status = system @cmd;
if ($status == -1) {
    die "[ERROR] Cannot execute sample shell generator: $!\n";
}
if ($status & 127) {
    die "[ERROR] Sample shell generator terminated by signal " . ($status & 127) . "\n";
}
my $exit_code = $status >> 8;
die "[ERROR] Sample shell generation failed with exit code $exit_code\n"
    if $exit_code != 0;

print "[DONE] Sample shell generation completed.\n";
print "[INFO] Each sample directory is directly under: $opt{outdir}\n";
print "[INFO] Run the generated <sample>.run.sh files separately.\n";
my $group_table = $conf{GROUPNAME_TABLE}
    // $conf{HUMAN_DNA_VIRUS_DB};
my $positive_cutoff = $conf{POSITIVE_HIT_CUTOFF} // 2;
if (defined $group_table && $group_table ne "") {
    my $summary_script = "$bin_dir/summarizeGroupnameResults.pl";
    print "[INFO] Run Summary after all sample parsed files are ready:\n";
    print '       perl "' . $summary_script . '" --input "' . $opt{outdir} . '" '
        . '--group-table "' . $group_table . '" --outdir "' . $opt{outdir} . '/Summary" '
        . '--positive-hit-cutoff "' . $positive_cutoff . '"' . "\n";
} else {
    print "[INFO] Run Summary after all sample parsed files are ready; "
        . "provide the Groupname table with --group-table.\n";
}

exit 0;

sub read_config {
    my ($file) = @_;
    my %cfg;
    open my $fh, '<', $file or die "[ERROR] Cannot open config: $file: $!\n";
    while (my $line = <$fh>) {
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;
        my ($key, $value) = split(/\s*=\s*/, $line, 2);
        defined $key && defined $value
            or die "[ERROR] Invalid config line: $line\n";
        $key =~ s/^\s+|\s+$//g;
        $value =~ s/^\s+|\s+$//g;
        $cfg{$key} = $value;
    }
    close $fh;
    return %cfg;
}

sub validate_config {
    my ($cfg) = @_;
    my @required = qw(
        BWA SAMTOOLS BLASTN HG38_REF VIRAL_DB LOCUS_ORGANISM
        BWA_THREADS BLAST_THREADS MAX_TARGET_SEQS MAX_HSPS
        IDENTITY_CUTOFF ALIGNMENT_LENGTH_CUTOFF EVALUE_CUTOFF
        KEEP_INTERMEDIATE
    );
    for my $key (@required) {
        die "[ERROR] Missing required config key: $key\n"
            unless exists $cfg->{$key} && defined $cfg->{$key} && $cfg->{$key} ne '';
    }
    die "[ERROR] KEEP_INTERMEDIATE must be 0 or 1\n"
        unless $cfg->{KEEP_INTERMEDIATE} =~ /\A[01]\z/;
}

sub check_required_paths {
    my ($cfg) = @_;

    for my $exe (qw(BWA SAMTOOLS BLASTN)) {
        -x $cfg->{$exe}
            or die "[ERROR] Executable not found or not executable: $cfg->{$exe}\n";
    }

    for my $file_key (qw(HG38_REF LOCUS_ORGANISM)) {
        -f $cfg->{$file_key}
            or die "[ERROR] Required file not found: $cfg->{$file_key}\n";
    }

    my $group_table = $cfg->{GROUPNAME_TABLE}
        // $cfg->{HUMAN_DNA_VIRUS_DB};
    if (defined $group_table && $group_table ne "") {
        -f $group_table
            or die "[ERROR] Groupname table not found: $group_table\n";
    }

    my $db_prefix = $cfg->{VIRAL_DB};
    my $db_found = 0;
    for my $suffix (qw(.nhr .nin .nsq .ndb)) {
        if (-f $db_prefix . $suffix) {
            $db_found = 1;
            last;
        }
    }
    die "[ERROR] BLAST viral database not found for prefix: $db_prefix\n"
        unless $db_found;
}

sub normalize_path {
    my ($path) = @_;
    return abs_path($path) if -d $path;
    return $path if $path =~ m{\A/};
    return getcwd() . "/$path";
}

sub usage {
    my ($exit_code) = @_;
    print STDERR <<'USAGE';
Usage:
  perl ViruScanNIPT.pl \
    --config conf/ViruScanNIPT.conf \
    --fq-list fq.list \
    --outdir ViruScanNIPT_result

Required arguments:
  --config     Pipeline configuration file
  --fq-list    FASTQ list, one FASTQ path per line
  --outdir     Output directory

Behavior:
  - Creates one sample directory directly under --outdir.
  - Writes <sample>.run.sh inside each sample directory.
  - Does not create summary-shell, log, or intermediate directory trees.
  - Does not execute the generated sample shells.

Run cohort-level summary separately after sample analysis:
  perl bin/summarizeGroupnameResults.pl \
    --input ViruScanNIPT_result \
    --group-table database/annotation/viralGroupname.tsv \
    --outdir ViruScanNIPT_result/Summary \
    --positive-hit-cutoff 2
USAGE
    exit($exit_code);
}
