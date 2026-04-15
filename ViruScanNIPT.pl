#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(dirname);
use Cwd qw(abs_path getcwd);

#==================================================
# Script Name:
#   ViruScanNIPT.pl
#
# Description:
#   Main framework generator for ViruScanNIPT pipeline.
#
#   This script does NOT execute the analysis directly.
#   Instead, it generates:
#     1. per-sample shell scripts
#     2. summary-level shell scripts
#
#   After running this script, users can execute all shell scripts
#   step by step according to their own compute environment.
#
#   Sample-level shell:
#     bwa mem
#     -> SAM/BAM conversion
#     -> sort/index
#     -> unmapped read extraction
#     -> unmapped BAM to FASTA
#     -> BLAST to viral DB
#     -> parse BLAST result
#
#   Summary-level shell:
#     01 merge parsed results
#     02 annotate human DNA virus
#     03 call group status
#     04 group positive ranking
#     05 codetection analysis
#     06 group distribution stats
#     07 visualize group distribution
#
#==================================================

my %opt;
GetOptions(
    "config=s"  => \$opt{config},
    "fq-list=s" => \$opt{fq_list},
    "outdir=s"  => \$opt{outdir},
    "help"      => \$opt{help},
) or usage();

usage() if $opt{help};

for my $arg (qw/config fq_list outdir/) {
    die "[ERROR] --$arg is required\n" unless defined $opt{$arg};
}

$opt{config}  = abs_path($opt{config});
$opt{fq_list} = abs_path($opt{fq_list});
$opt{outdir}  = normalize_path($opt{outdir});

-f $opt{config}  or die "[ERROR] config file not found: $opt{config}\n";
-f $opt{fq_list} or die "[ERROR] fq.list not found: $opt{fq_list}\n";

my $script_dir = abs_path(dirname($0));
my $bin_dir    = "$script_dir/bin";
-d $bin_dir or die "[ERROR] bin directory not found: $bin_dir\n";

my %conf = read_config($opt{config});
validate_config(\%conf);
check_required_paths(\%conf);

my %step_script = (
    make_sample_shells         => "$bin_dir/VSN_01_make_sample_shells.pl",
    parse_blast                => "$bin_dir/VSN_05_parse_blast.pl",
    merge_sample_results       => "$bin_dir/VSN_06_merge_sample_results.pl",
    annotate_human_dna_virus   => "$bin_dir/VSN_07_annotate_human_dna_virus.pl",
    call_group_status          => "$bin_dir/VSN_08_call_group_status.pl",
    group_positive_ranking     => "$bin_dir/VSN_09_group_positive_ranking.pl",
    codetection_analysis       => "$bin_dir/VSN_10_codetection_analysis.pl",
    group_distribution_stats   => "$bin_dir/VSN_11_group_distribution_stats.pl",
    virostat_visualizer        => "$bin_dir/ViroStat_Visualizer.py",
);


create_dir_structure($opt{outdir});
my %dir = build_dir_map($opt{outdir});

save_runtime_files(
    outdir     => $opt{outdir},
    fq_list    => $opt{fq_list},
    config     => $opt{config},
    conf_ref   => \%conf,
    dir_ref    => \%dir,
);

log_message("$dir{logs}/ViruScanNIPT.pipeline.log", "Pipeline framework generation started");
log_message("$dir{logs}/ViruScanNIPT.pipeline.log", "Config: $opt{config}");
log_message("$dir{logs}/ViruScanNIPT.pipeline.log", "fq.list: $opt{fq_list}");
log_message("$dir{logs}/ViruScanNIPT.pipeline.log", "outdir: $opt{outdir}");

run_make_sample_shells(
    script                  => $step_script{make_sample_shells},
    fq_list                 => $opt{fq_list},
    outdir                  => $dir{sample_shells},
    bin_dir                 => $bin_dir,
    bwa                     => $conf{BWA},
    samtools                => $conf{SAMTOOLS},
    blastn                  => $conf{BLASTN},
    ref                     => $conf{HG38_REF},
    viral_db                => $conf{VIRAL_DB},
    locus_organism          => $conf{LOCUS_ORGANISM},
    bwa_threads             => $conf{BWA_THREADS},
    blast_threads           => $conf{BLAST_THREADS},
    max_target_seqs         => $conf{MAX_TARGET_SEQS},
    max_hsps                => $conf{MAX_HSPS},
    identity_cutoff         => $conf{IDENTITY_CUTOFF},
    alignment_length_cutoff => $conf{ALIGNMENT_LENGTH_CUTOFF},
    evalue_cutoff           => $conf{EVALUE_CUTOFF},
    keep_intermediate       => $conf{KEEP_INTERMEDIATE},
    log_file                => "$dir{logs}/make_sample_shells.log",
);

generate_summary_shells(
    outdir                    => $dir{summary_shells},
    logs_dir                  => $dir{logs},
    summary_merge_dir         => $dir{summary_merge},
    summary_annot_dir         => $dir{summary_annot},
    summary_group_status_dir  => $dir{summary_group_status},
    summary_group_ranking_dir => $dir{summary_group_ranking},
    summary_codetection_dir   => $dir{summary_codetection},
    summary_group_dist_dir    => $dir{summary_group_distribution},
    sample_shells_dir         => $dir{sample_shells},
    merge_script              => $step_script{merge_sample_results},
    annotate_script           => $step_script{annotate_human_dna_virus},
    group_status_script       => $step_script{call_group_status},
    ranking_script            => $step_script{group_positive_ranking},
    codetection_script        => $step_script{codetection_analysis},
    distribution_script       => $step_script{group_distribution_stats},
    visualizer_script         => $step_script{virostat_visualizer},
    python3                   => $conf{PYTHON3},
    human_dna_virus_db        => $conf{HUMAN_DNA_VIRUS_DB},
    positive_hit_cutoff       => $conf{POSITIVE_HIT_CUTOFF},
);

log_message("$dir{logs}/ViruScanNIPT.pipeline.log", "Pipeline framework generation finished successfully");
exit(0);

#==================================================
# subroutines
#==================================================

sub usage {
    print STDERR <<'USAGE';
Usage:
  perl ViruScanNIPT.pl --config ViruScanNIPT.conf --fq-list fq.list --outdir result_dir

Required arguments:
  --config     ViruScanNIPT config file
  --fq-list    Input fastq list, one fastq path per line
  --outdir     Output directory

Optional:
  --help       Print this help message

Notes:
  1. This script only generates shell scripts and directory framework.
  2. It does NOT execute sample-level or summary-level analysis.
  3. Users should run generated shell scripts step by step.

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

sub read_config {
    my ($file) = @_;
    my %cfg;

    open my $fh, "<", $file or die "[ERROR] Cannot open config file: $file\n";
    while (<$fh>) {
        chomp;
        s/^\s+//;
        s/\s+$//;
        next if /^$/;
        next if /^#/;

        my ($k, $v) = split(/\s*=\s*/, $_, 2);
        die "[ERROR] Bad config line: $_\n" unless defined $k && defined $v;

        $k =~ s/^\s+//;
        $k =~ s/\s+$//;
        $v =~ s/^\s+//;
        $v =~ s/\s+$//;

        $cfg{$k} = $v;
    }
    close $fh;

    return %cfg;
}

sub validate_config {
    my ($conf_ref) = @_;

    my @required = qw(
        BWA
        SAMTOOLS
        BLASTN
        PYTHON3
        HG38_REF
        VIRAL_DB
        LOCUS_ORGANISM
        HUMAN_DNA_VIRUS_DB
        BWA_THREADS
        BLAST_THREADS
        EVALUE_CUTOFF
        IDENTITY_CUTOFF
        ALIGNMENT_LENGTH_CUTOFF
        POSITIVE_HIT_CUTOFF
        MAX_TARGET_SEQS
        MAX_HSPS
        KEEP_INTERMEDIATE
    );

    for my $k (@required) {
        die "[ERROR] Missing required config key: $k\n" unless exists $conf_ref->{$k};
    }
}

sub check_required_paths {
    my ($conf_ref) = @_;

    for my $exe (qw/BWA SAMTOOLS BLASTN/) {
        -x $conf_ref->{$exe} or die "[ERROR] Executable not found or not executable: $conf_ref->{$exe}\n";
    }

    if ($conf_ref->{PYTHON3} =~ m{/}) {
        -x $conf_ref->{PYTHON3} or die "[ERROR] Executable not found or not executable: $conf_ref->{PYTHON3}\n";
    }

    for my $f (qw/HG38_REF LOCUS_ORGANISM HUMAN_DNA_VIRUS_DB/) {
        -f $conf_ref->{$f} or die "[ERROR] Required file not found: $conf_ref->{$f}\n";
    }

    my $viral_db_found = 0;
    for my $suffix (qw/.nhr .nin .nsq/) {
        if (-f $conf_ref->{VIRAL_DB} . $suffix) {
            $viral_db_found = 1;
            last;
        }
    }
    die "[ERROR] BLAST viral database not found: $conf_ref->{VIRAL_DB}\n" unless $viral_db_found;
}

sub create_dir_structure {
    my ($outdir) = @_;

    my @dirs = (
        $outdir,
        "$outdir/00_run_config",
        "$outdir/01_sample_shells",
        "$outdir/02_summary_shells",
        "$outdir/03_summary_analysis",
        "$outdir/03_summary_analysis/01_merge",
        "$outdir/03_summary_analysis/02_human_dna_virus_annotation",
        "$outdir/03_summary_analysis/03_group_status",
        "$outdir/03_summary_analysis/04_group_ranking",
        "$outdir/03_summary_analysis/05_codetection",
        "$outdir/03_summary_analysis/06_group_distribution",
        "$outdir/04_logs",
    );

    for my $d (@dirs) {
        mkdir_if_not_exists($d);
    }
}

sub build_dir_map {
    my ($outdir) = @_;

    my %dir = (
        root                       => $outdir,
        run_config                 => "$outdir/00_run_config",
        sample_shells              => "$outdir/01_sample_shells",
        summary_shells             => "$outdir/02_summary_shells",
        summary_root               => "$outdir/03_summary_analysis",
        summary_merge              => "$outdir/03_summary_analysis/01_merge",
        summary_annot              => "$outdir/03_summary_analysis/02_human_dna_virus_annotation",
        summary_group_status       => "$outdir/03_summary_analysis/03_group_status",
        summary_group_ranking      => "$outdir/03_summary_analysis/04_group_ranking",
        summary_codetection        => "$outdir/03_summary_analysis/05_codetection",
        summary_group_distribution => "$outdir/03_summary_analysis/06_group_distribution",
        logs                       => "$outdir/04_logs",
    );

    return %dir;
}

sub save_runtime_files {
    my %args = @_;
    my $fq_list  = $args{fq_list};
    my $config   = $args{config};
    my %conf     = %{$args{conf_ref}};
    my %dir      = %{$args{dir_ref}};

    system("cp", $fq_list, "$dir{run_config}/fq.list") == 0
        or die "[ERROR] failed to copy fq.list into run_config\n";

    open my $cmd_fh, ">", "$dir{run_config}/ViruScanNIPT.command.txt"
        or die "[ERROR] cannot write command record\n";
    print $cmd_fh join(" ", $0, @ARGV), "\n";
    close $cmd_fh;

    open my $run_fh, ">", "$dir{run_config}/ViruScanNIPT.run.config.txt"
        or die "[ERROR] cannot write runtime config file\n";

    print $run_fh "PROJECT_NAME=ViruScanNIPT\n";
    print $run_fh "CONFIG_FILE=$config\n";
    print $run_fh "FQ_LIST=$fq_list\n";
    print $run_fh "OUTDIR=$args{outdir}\n";

    for my $k (sort keys %conf) {
        print $run_fh "$k=$conf{$k}\n";
    }
    close $run_fh;
}

sub run_make_sample_shells {
    my %args = @_;

    my $cmd = join(" ",
        "perl",
        shell_quote($args{script}),
        "--fq-list", shell_quote($args{fq_list}),
        "--outdir", shell_quote($args{outdir}),
        "--bin-dir", shell_quote($args{bin_dir}),
        "--bwa", shell_quote($args{bwa}),
        "--samtools", shell_quote($args{samtools}),
        "--blastn", shell_quote($args{blastn}),
        "--ref", shell_quote($args{ref}),
        "--viral-db", shell_quote($args{viral_db}),
        "--locus-organism", shell_quote($args{locus_organism}),
        "--bwa-threads", shell_quote($args{bwa_threads}),
        "--blast-threads", shell_quote($args{blast_threads}),
        "--max-target-seqs", shell_quote($args{max_target_seqs}),
        "--max-hsps", shell_quote($args{max_hsps}),
        "--identity-cutoff", shell_quote($args{identity_cutoff}),
        "--alignment-length-cutoff", shell_quote($args{alignment_length_cutoff}),
        "--evalue-cutoff", shell_quote($args{evalue_cutoff}),
        "--keep-intermediate", shell_quote($args{keep_intermediate}),
        ">", shell_quote($args{log_file}),
        "2>&1"
    );
    run_cmd($cmd);
}

sub generate_summary_shells {
    my %args = @_;

    my @jobs = (
        {
            file => "$args{outdir}/01_merge_sample_results.sh",
            need => $args{merge_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{merge_script}),
                "--sample-dir", shell_quote($args{sample_shells_dir}),
                "--outdir", shell_quote($args{summary_merge_dir}),
                ">", shell_quote("$args{logs_dir}/summary.merge.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/02_annotate_human_dna_virus.sh",
            need => $args{annotate_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{annotate_script}),
                "--input", shell_quote("$args{summary_merge_dir}/merge_blast_parsed_res.xls"),
                "--outdir", shell_quote($args{summary_annot_dir}),
                "--virus-db", shell_quote($args{human_dna_virus_db}),
                ">", shell_quote("$args{logs_dir}/summary.annotate_human_dna_virus.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/03_call_group_status.sh",
            need => $args{group_status_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{group_status_script}),
                "--input", shell_quote("$args{summary_annot_dir}/Annotated_Human_Virus_Results.xls"),
                "--outdir", shell_quote($args{summary_group_status_dir}),
                "--positive-hit-cutoff", shell_quote($args{positive_hit_cutoff}),
                ">", shell_quote("$args{logs_dir}/summary.group_status.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/04_group_positive_ranking.sh",
            need => $args{ranking_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{ranking_script}),
                "--input", shell_quote("$args{summary_group_status_dir}/Virus_Group_Statistics_Final.xls"),
                "--outdir", shell_quote($args{summary_group_ranking_dir}),
                ">", shell_quote("$args{logs_dir}/summary.group_ranking.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/05_codetection_analysis.sh",
            need => $args{codetection_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{codetection_script}),
                "--input", shell_quote("$args{summary_group_status_dir}/Virus_Group_Statistics_Final.xls"),
                "--outdir", shell_quote($args{summary_codetection_dir}),
                ">", shell_quote("$args{logs_dir}/summary.codetection.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/06_group_distribution_stats.sh",
            need => $args{distribution_script},
            cmd  => join(" ",
                "perl",
                shell_quote($args{distribution_script}),
                "--input", shell_quote("$args{summary_group_status_dir}/Virus_Group_Statistics_Final.xls"),
                "--outdir", shell_quote($args{summary_group_dist_dir}),
                ">", shell_quote("$args{logs_dir}/summary.group_distribution.log"),
                "2>&1"
            ),
        },
        {
            file => "$args{outdir}/07_visualize_group_distribution.sh",
            need => $args{visualizer_script},
            cmd  => join(" ",
                shell_quote($args{python3}),
                shell_quote($args{visualizer_script}),
                shell_quote("$args{summary_group_dist_dir}/Groupname_Sample_Count_Detail.xls"),
                ">", shell_quote("$args{logs_dir}/summary.visualize_group_distribution.log"),
                "2>&1"
            ),
        },
    );

    for my $job (@jobs) {
        open my $sh, ">", $job->{file} or die "[ERROR] cannot write shell file: $job->{file}\n";
        print $sh "#!/bin/bash\n";
        print $sh "set -euo pipefail\n\n";

        if (defined $job->{need} && -f $job->{need}) {
            print $sh "$job->{cmd}\n";
        } else {
            print $sh "echo \"[ERROR] required script not found: $job->{need}\" >&2\n";
            print $sh "exit 1\n";
        }

        close $sh;
        chmod 0755, $job->{file} or die "[ERROR] cannot chmod 755 for $job->{file}\n";
    }
}

sub run_cmd {
    my ($cmd) = @_;
    my $ret = system($cmd);
    if ($ret != 0) {
        die "[ERROR] command failed:\n$cmd\n";
    }
}

sub shell_quote {
    my ($str) = @_;
    $str =~ s/'/'"'"'/g;
    return "'$str'";
}

sub mkdir_if_not_exists {
    my ($dir) = @_;
    return if -d $dir;
    mkdir $dir or die "[ERROR] cannot create directory: $dir\n";
}

sub log_message {
    my ($logfile, $msg) = @_;
    open my $fh, ">>", $logfile or die "[ERROR] cannot write log: $logfile\n";
    my $time = localtime();
    print $fh "[$time] $msg\n";
    close $fh;
}

