#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);

# ViruScanNIPT.pl generates one executable sample-level shell per FASTQ.
# It does not execute sample shells or run cohort Summary automatically.
# All configuration keys in ViruScanNIPT.conf are mandatory and must match
# the current conf/ViruScanNIPT.conf.example exactly.

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

die "[ERROR] Config file not found: $opt{config}\n" unless -f $opt{config};
die "[ERROR] FASTQ list not found: $opt{fq_list}\n" unless -f $opt{fq_list};
die "[ERROR] --outdir cannot contain a single quote or newline\n"
    if $opt{outdir} =~ /['\r\n]/;

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
print "[INFO] Configuration  : $opt{config}\n";
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
my $summary_script = "$bin_dir/summarizeGroupnameResults.pl";
print "[INFO] Run Summary after all sample parsed files are ready:\n";
print '       perl "' . $summary_script . '" --input "' . $opt{outdir} . '" '
    . '--group-table "' . $conf{GROUPNAME_TABLE} . '" --outdir "' . $opt{outdir} . '/Summary" '
    . '--positive-hit-cutoff "' . $conf{POSITIVE_HIT_CUTOFF} . '"' . "\n";

exit 0;

sub read_config {
    my ($file) = @_;
    my %cfg;
    my %allowed = map { $_ => 1 } qw(
        BWA SAMTOOLS BLASTN HG38_REF VIRAL_DB LOCUS_ORGANISM GROUPNAME_TABLE
        BWA_THREADS BLAST_THREADS MAX_TARGET_SEQS MAX_HSPS
        EVALUE_CUTOFF IDENTITY_CUTOFF ALIGNMENT_LENGTH_CUTOFF
        POSITIVE_HIT_CUTOFF KEEP_INTERMEDIATE
    );

    open my $fh, '<', $file
        or die "[ERROR] Cannot open config: $file: $!\n";

    my $line_number = 0;
    while (my $line = <$fh>) {
        ++$line_number;
        chomp $line;
        $line =~ s/\r$//;
        $line =~ s/^\s+|\s+$//g;
        next if $line eq '' || $line =~ /^#/;

        my ($key, $value) = $line =~ /\A([A-Z][A-Z0-9_]*)\s*=\s*(.*?)\s*\z/;
        defined $key
            or die "[ERROR] Invalid config syntax at line $line_number: $line\n";
        exists $allowed{$key}
            or die "[ERROR] Unknown config key '$key' at line $line_number. "
                . "Use only keys defined in the current configuration template.\n";
        exists $cfg{$key}
            and die "[ERROR] Duplicate config key '$key' at line $line_number\n";
        defined $value && length $value
            or die "[ERROR] Empty value for config key '$key' at line $line_number\n";

        $cfg{$key} = $value;
    }
    close $fh or die "[ERROR] Cannot close config file $file: $!\n";
    return %cfg;
}

sub validate_config {
    my ($cfg) = @_;
    my @required = qw(
        BWA SAMTOOLS BLASTN HG38_REF VIRAL_DB LOCUS_ORGANISM GROUPNAME_TABLE
        BWA_THREADS BLAST_THREADS MAX_TARGET_SEQS MAX_HSPS
        EVALUE_CUTOFF IDENTITY_CUTOFF ALIGNMENT_LENGTH_CUTOFF
        POSITIVE_HIT_CUTOFF KEEP_INTERMEDIATE
    );

    for my $key (@required) {
        die "[ERROR] Missing required config key: $key\n"
            unless exists $cfg->{$key} && defined $cfg->{$key} && $cfg->{$key} ne '';
    }

    for my $key (qw(
        BWA_THREADS BLAST_THREADS MAX_TARGET_SEQS MAX_HSPS
        ALIGNMENT_LENGTH_CUTOFF POSITIVE_HIT_CUTOFF
    )) {
        die "[ERROR] $key must be a positive integer\n"
            unless $cfg->{$key} =~ /\A[1-9][0-9]*\z/;
    }

    die "[ERROR] IDENTITY_CUTOFF must be a number greater than 0 and <= 100\n"
        unless $cfg->{IDENTITY_CUTOFF} =~ /\A(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)\z/
            && $cfg->{IDENTITY_CUTOFF} > 0
            && $cfg->{IDENTITY_CUTOFF} <= 100;

    die "[ERROR] EVALUE_CUTOFF must be a positive numeric value (for example, 1e-5)\n"
        unless $cfg->{EVALUE_CUTOFF} =~ /\A(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\z/
            && $cfg->{EVALUE_CUTOFF} > 0;

    die "[ERROR] KEEP_INTERMEDIATE must be 0 or 1\n"
        unless $cfg->{KEEP_INTERMEDIATE} =~ /\A[01]\z/;

    for my $key (qw(BWA SAMTOOLS BLASTN HG38_REF VIRAL_DB LOCUS_ORGANISM GROUPNAME_TABLE)) {
        die "[ERROR] $key must be an absolute path: $cfg->{$key}\n"
            unless $cfg->{$key} =~ m{\A/};
        die "[ERROR] $key cannot contain a single quote or newline\n"
            if $cfg->{$key} =~ /['\r\n]/;
    }
}

sub check_required_paths {
    my ($cfg) = @_;

    for my $exe (qw(BWA SAMTOOLS BLASTN)) {
        -x $cfg->{$exe}
            or die "[ERROR] Executable not found or not executable: $cfg->{$exe}\n";
    }

    for my $file_key (qw(HG38_REF LOCUS_ORGANISM GROUPNAME_TABLE)) {
        -s $cfg->{$file_key}
            or die "[ERROR] Required file not found or empty: $cfg->{$file_key}\n";
    }

    # BWA-MEM requires the five classic BWA index components beside the FASTA.
    for my $suffix (qw(amb ann bwt pac sa)) {
        -s "$cfg->{HG38_REF}.$suffix"
            or die "[ERROR] BWA index component is missing: $cfg->{HG38_REF}.$suffix\n";
    }

    # Accept a nucleotide BLAST database alias (.nal), or a complete classic
    # nucleotide database index set (.nhr/.nin/.nsq; also present in BLAST DB v5).
    my $db_prefix = $cfg->{VIRAL_DB};
    my $has_alias = -s "$db_prefix.nal";
    my $has_index = 1;
    for my $suffix (qw(nhr nin nsq)) {
        $has_index = 0 unless -s "$db_prefix.$suffix";
    }
    die "[ERROR] Viral nucleotide BLAST database not found or incomplete for prefix: $db_prefix\n"
        unless $has_alias || $has_index;

    open my $GROUP, '<', $cfg->{GROUPNAME_TABLE}
        or die "[ERROR] Cannot open Groupname table: $cfg->{GROUPNAME_TABLE}: $!\n";
    my $header = <$GROUP>;
    defined $header
        or die "[ERROR] Groupname table is empty: $cfg->{GROUPNAME_TABLE}\n";
    $header =~ s/\r?\n\z//;
    my @columns = split /\t/, $header, -1;
    my %column;
    for my $i (0 .. $#columns) {
        die "[ERROR] Duplicate Groupname table column '$columns[$i]'\n"
            if exists $column{$columns[$i]};
        $column{$columns[$i]} = $i;
    }
    for my $required_column (qw(ACCESSION Groupname)) {
        die "[ERROR] Groupname table must contain column '$required_column': $cfg->{GROUPNAME_TABLE}\n"
            unless exists $column{$required_column};
    }

    my %seen_accession;
    my $data_rows = 0;
    while (my $line = <$GROUP>) {
        chomp $line;
        $line =~ s/\r\z//;
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
        my @fields = split /\t/, $line, -1;
        my $accession = $fields[$column{ACCESSION}] // '';
        my $groupname = $fields[$column{Groupname}] // '';
        $accession =~ s/^\s+|\s+$//g;
        $groupname =~ s/^\s+|\s+$//g;
        die "[ERROR] Blank ACCESSION or Groupname in $cfg->{GROUPNAME_TABLE}\n"
            unless length $accession && length $groupname;
        die "[ERROR] Duplicate ACCESSION '$accession' in $cfg->{GROUPNAME_TABLE}\n"
            if $seen_accession{$accession}++;
        ++$data_rows;
    }
    close $GROUP or die "[ERROR] Cannot close Groupname table: $!\n";
    die "[ERROR] Groupname table contains no data rows: $cfg->{GROUPNAME_TABLE}\n"
        unless $data_rows;
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
  --config     Current-format pipeline configuration file (all keys required)
  --fq-list    FASTQ list, one single-end FASTQ path per line
  --outdir     Output directory

Behavior:
  - Accepts only keys defined in conf/ViruScanNIPT.conf.example.
  - Rejects unknown, duplicate, missing or empty configuration entries.
  - Requires absolute paths and checks BWA/BLAST database indexes.
  - Creates one sample directory directly under --outdir.
  - Writes <sample>.run.sh inside each sample directory.
  - Does not execute the generated sample shells.

Run cohort-level Summary separately after all parsed files are ready.
USAGE
    exit($exit_code);
}
