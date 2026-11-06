#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Find qw(find);
use File::Path qw(make_path);
use Cwd qw(abs_path getcwd);
use POSIX qw(strftime);

# ============================================================
# summarizeGroupnameResults.pl
#
# Summarize sample-level viral Groupname detections directly
# from one or more *_blast_out_parsed.txt files.
#
# The parsed BLAST file is joined to a curated Groupname table
# by subject accession. Accession version suffixes are ignored.
#
# Counting unit:
#   One distinct query per sample per Groupname.
#
# Positive definition:
#   Distinct query count >= positive-hit-cutoff (default: 2).
#
# Supported Groupname table headers:
#   ACCESSION / Reference_accession + Groupname
#
# Output:
#   sampleGroupCounts.tsv
#   sampleSummary.tsv
#   groupSummary.tsv
#   groupReadCounts.tsv
#   groupCoDetection.tsv
#   unmappedAccessions.tsv
#   summary.log
# ============================================================

my %opt;
GetOptions(
    'input=s'                  => \$opt{input},
    'group-table=s'            => \$opt{group_table},
    'outdir=s'                 => \$opt{outdir},
    'positive-hit-cutoff|positive-cutoff=i' => \$opt{positive_hit_cutoff},
    'help'                     => \$opt{help},
) or usage(1);

usage(0) if $opt{help};

for my $key (qw(input group_table outdir)) {
    die "[ERROR] --$key is required\n" unless defined $opt{$key};
}

$opt{positive_hit_cutoff} = 2 unless defined $opt{positive_hit_cutoff};
die "[ERROR] --positive-hit-cutoff must be >= 1\n"
    unless $opt{positive_hit_cutoff} =~ /\A[1-9][0-9]*\z/;

$opt{input} = abs_path($opt{input})
    or die "[ERROR] Cannot resolve input: $opt{input}\n";
$opt{group_table} = abs_path($opt{group_table})
    or die "[ERROR] Cannot resolve Groupname table: $opt{group_table}\n";
$opt{outdir} = normalize_output_path($opt{outdir});

die "[ERROR] Input does not exist: $opt{input}\n" unless -e $opt{input};
die "[ERROR] Groupname table not found: $opt{group_table}\n"
    unless -f $opt{group_table};

make_path($opt{outdir}) unless -d $opt{outdir};

my ($accession_to_group_ref, $group_set_ref) = read_group_table($opt{group_table});
my %accession_to_group = %{$accession_to_group_ref};
my %group_set = %{$group_set_ref};
my @groups = sort keys %group_set;

my @parsed_files = collect_parsed_files($opt{input});
die "[ERROR] No *_blast_out_parsed.txt files found under: $opt{input}\n"
    unless @parsed_files;

check_incomplete_sample_directories($opt{input}, \@parsed_files)
    if -d $opt{input};

my (%sample_seen, %sample_queries);
my (%unmapped_record_count, %unmapped_pairs);
my ($parsed_records, $mapped_records, $unmapped_records, $phage_filtered_records) = (0, 0, 0, 0);
my %sample_source_file_count;

for my $file (@parsed_files) {
    my ($column_ref, $fallback_sample) = read_parsed_header($file);
    my %column = %{$column_ref};

    open my $IN, '<', $file
        or die "[ERROR] Cannot open parsed file $file: $!\n";
    my $header = <$IN>;
    defined $header or die "[ERROR] Parsed file is empty: $file\n";

    my $file_has_data_sample = 0;

    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/;

        my @field = split /\t/, $line, -1;
        my $sample = field_value(\@field, $column{sample});
        my $query = field_value(\@field, $column{query});
        my $subject = field_value(\@field, $column{subject});
        my $viral_name = field_value(\@field, $column{viral_name});

        $sample = $fallback_sample
            unless defined $sample && $sample ne '' && $sample ne 'NA';
        die "[ERROR] Sample name is missing in $file\n"
            unless defined $sample && $sample ne '';

        $sample_seen{$sample} = 1;
        $file_has_data_sample = 1;
        $sample_source_file_count{$sample}{$file} = 1;

        next unless defined $query && $query ne '' && $query ne 'NA';

        $parsed_records++;

        # Retain the previous workflow's explicit phage exclusion.
        if (defined $viral_name && $viral_name =~ /phage/i) {
            $phage_filtered_records++;
            next;
        }

        my $accession = normalize_accession($subject);

        if (!length($accession) || !exists $accession_to_group{$accession}) {
            $unmapped_records++;
            if (length($accession)) {
                $unmapped_record_count{$accession}++;
                $unmapped_pairs{$accession}{$sample}{$query} = 1;
            }
            next;
        }

        my $group = $accession_to_group{$accession};
        $mapped_records++;

        # Multiple reference accessions within one Groupname must not
        # inflate the read count for the same sample-query pair.
        $sample_queries{$sample}{$group}{$query} = 1;
    }
    close $IN;

    # The parser always writes a header, including when no hit passes filters.
    # In that case, infer the sample ID from the parsed filename.
    if (!$file_has_data_sample) {
        $sample_seen{$fallback_sample} = 1;
        $sample_source_file_count{$fallback_sample}{$file} = 1;
    }
}

my @samples = sort keys %sample_seen;
my $total_samples = scalar @samples;
die "[ERROR] No sample IDs could be recovered from parsed files\n"
    unless $total_samples;

# ------------------------------------------------------------
# Sample x Groupname count matrix in long-table format
# ------------------------------------------------------------
my (%group_positive_samples, %sample_positive_groups);
my (%group_counts_by_sample, %group_total_queries);
my (%sample_total_distinct_queries, %sample_detected_groups);

my $sample_group_file = "$opt{outdir}/sampleGroupCounts.tsv";
open my $SG, '>', $sample_group_file
    or die "[ERROR] Cannot write $sample_group_file: $!\n";
print {$SG} join("\t", qw(Sample Groupname DistinctQueryCount Status)), "\n";

for my $sample (@samples) {
    my %sample_query_union;

    for my $group (@groups) {
        my $count = 0;

        if (exists $sample_queries{$sample}{$group}) {
            $count = scalar keys %{$sample_queries{$sample}{$group}};
            $group_total_queries{$group} += $count;
            $sample_query_union{$_} = 1
                for keys %{$sample_queries{$sample}{$group}};
            $sample_detected_groups{$sample}{$group} = 1 if $count > 0;
        }

        $group_counts_by_sample{$group}{$sample} = $count;

        my $status = $count >= $opt{positive_hit_cutoff}
            ? 'Positive' : 'Negative';

        if ($status eq 'Positive') {
            $group_positive_samples{$group}{$sample} = 1;
            $sample_positive_groups{$sample}{$group} = 1;
        }

        print {$SG} join("\t", $sample, $group, $count, $status), "\n";
    }

    $sample_total_distinct_queries{$sample} = scalar keys %sample_query_union;
}
close $SG;

# ------------------------------------------------------------
# Per-sample summary
# ------------------------------------------------------------
my $sample_summary_file = "$opt{outdir}/sampleSummary.tsv";
open my $SS, '>', $sample_summary_file
    or die "[ERROR] Cannot write $sample_summary_file: $!\n";
print {$SS} join("\t", qw(
    Sample TotalDistinctAssignedQueries
    GroupnamesWithAtLeastOneRead PositiveGroupCount PositiveGroupnames
)), "\n";

for my $sample (@samples) {
    my @detected = sort keys %{$sample_detected_groups{$sample} || {}};
    my @positive = sort keys %{$sample_positive_groups{$sample} || {}};

    print {$SS} join("\t",
        $sample,
        ($sample_total_distinct_queries{$sample} || 0),
        scalar(@detected),
        scalar(@positive),
        (@positive ? join(';', @positive) : 'NA')
    ), "\n";
}
close $SS;

# ------------------------------------------------------------
# Group-level summary and positive sample read-count detail
# ------------------------------------------------------------
my $group_summary_file = "$opt{outdir}/groupSummary.tsv";
my $group_read_file = "$opt{outdir}/groupReadCounts.tsv";

open my $GS, '>', $group_summary_file
    or die "[ERROR] Cannot write $group_summary_file: $!\n";
print {$GS} join("\t", qw(
    Groupname PositiveSampleCount TotalSampleCount PositiveRatePercent
    TotalDistinctQueries MeanReadsPerPositiveSample MedianReadsPerPositiveSample
    Q1 Q3 Min Max
)), "\n";

open my $GR, '>', $group_read_file
    or die "[ERROR] Cannot write $group_read_file: $!\n";
print {$GR} join("\t", qw(Groupname Sample DistinctQueryCount)), "\n";

for my $group (@groups) {
    my $positive_n = scalar keys %{$group_positive_samples{$group} || {}};
    my $rate = sprintf('%.4f', $positive_n / $total_samples * 100);

    my @positive_counts;
    for my $sample (@samples) {
        my $count = $group_counts_by_sample{$group}{$sample} || 0;
        next unless $count >= $opt{positive_hit_cutoff};
        push @positive_counts, $count;
        print {$GR} join("\t", $group, $sample, $count), "\n";
    }

    @positive_counts = sort { $a <=> $b } @positive_counts;
    my ($mean, $median, $q1, $q3, $min, $max) =
        ('NA', 'NA', 'NA', 'NA', 'NA', 'NA');

    if (@positive_counts) {
        my $sum = 0;
        $sum += $_ for @positive_counts;
        $mean = sprintf('%.2f', $sum / @positive_counts);
        $median = quantile(\@positive_counts, 0.50);
        $q1 = quantile(\@positive_counts, 0.25);
        $q3 = quantile(\@positive_counts, 0.75);
        $min = $positive_counts[0];
        $max = $positive_counts[-1];
    }

    print {$GS} join("\t",
        $group,
        $positive_n,
        $total_samples,
        $rate,
        ($group_total_queries{$group} || 0),
        $mean, $median, $q1, $q3, $min, $max
    ), "\n";
}
close $GS;
close $GR;

# ------------------------------------------------------------
# Pairwise Groupname co-detection
# ------------------------------------------------------------
my %pair_count;
for my $sample (@samples) {
    my @positive = sort keys %{$sample_positive_groups{$sample} || {}};
    for (my $i = 0; $i < @positive; $i++) {
        for (my $j = $i + 1; $j < @positive; $j++) {
            $pair_count{$positive[$i]}{$positive[$j]}++;
        }
    }
}

my $codetect_file = "$opt{outdir}/groupCoDetection.tsv";
open my $CO, '>', $codetect_file
    or die "[ERROR] Cannot write $codetect_file: $!\n";
print {$CO} join("\t", qw(
    Groupname1 Groupname2 CoDetectedSamples CohortPercent
)), "\n";

for my $g1 (@groups) {
    for my $g2 (@groups) {
        next unless $g1 lt $g2;
        my $n = $pair_count{$g1}{$g2} || 0;
        next unless $n;
        my $pct = sprintf('%.4f', $n / $total_samples * 100);
        print {$CO} join("\t", $g1, $g2, $n, $pct), "\n";
    }
}
close $CO;

# ------------------------------------------------------------
# Exact positive-Groupname combinations (legacy combination summary)
# ------------------------------------------------------------
my %combination_count;
my $sample_codetect_file = "$opt{outdir}/sampleCoDetection.tsv";
open my $SCD, '>', $sample_codetect_file
    or die "[ERROR] Cannot write $sample_codetect_file: $!\n";
print {$SCD} join("\t", qw(
    Sample PositiveGroupCount CoDetectionCombination
)), "\n";

for my $sample (@samples) {
    my @positive = sort keys %{$sample_positive_groups{$sample} || {}};
    next unless @positive >= 2;
    my $combination = join(" + ", @positive);
    print {$SCD} join("\t", $sample, scalar(@positive), $combination), "\n";
    $combination_count{$combination}++;
}
close $SCD;

my $combination_summary_file = "$opt{outdir}/groupCombinationSummary.tsv";
open my $GCS, '>', $combination_summary_file
    or die "[ERROR] Cannot write $combination_summary_file: $!\n";
print {$GCS} join("\t", qw(
    CoDetectionCombination PositiveSampleCount
)), "\n";

for my $combination (
    sort {
        $combination_count{$b} <=> $combination_count{$a}
            || $a cmp $b
    } keys %combination_count
) {
    print {$GCS} join("\t", $combination, $combination_count{$combination}), "\n";
}
close $GCS;

# ------------------------------------------------------------
# Unmapped accession report for Groupname-table quality control
# ------------------------------------------------------------
my $unmapped_file = "$opt{outdir}/unmappedAccessions.tsv";
open my $UM, '>', $unmapped_file
    or die "[ERROR] Cannot write $unmapped_file: $!\n";
print {$UM} join("\t", qw(
    ACCESSION UnmappedRecordCount DistinctSampleQueryPairs
)), "\n";

for my $accession (sort keys %unmapped_record_count) {
    my %pairs;
    for my $sample (keys %{$unmapped_pairs{$accession} || {}}) {
        $pairs{"$sample\t$_"} = 1
            for keys %{$unmapped_pairs{$accession}{$sample}};
    }
    print {$UM} join("\t",
        $accession,
        $unmapped_record_count{$accession},
        scalar(keys %pairs)
    ), "\n";
}
close $UM;

# ------------------------------------------------------------
# Summary log
# ------------------------------------------------------------
my $log_file = "$opt{outdir}/summary.log";
open my $LOG, '>', $log_file
    or die "[ERROR] Cannot write $log_file: $!\n";
print {$LOG} "ViruScanNIPT Groupname summary\n";
print {$LOG} "Date: ", strftime('%Y-%m-%d %H:%M:%S', localtime), "\n";
print {$LOG} "Input: $opt{input}\n";
print {$LOG} "Groupname table: $opt{group_table}\n";
print {$LOG} "Positive cutoff: >= $opt{positive_hit_cutoff} distinct query reads per sample per Groupname\n";
print {$LOG} "Parsed files: ", scalar(@parsed_files), "\n";
print {$LOG} "Unique samples: $total_samples\n";
print {$LOG} "Parsed records: $parsed_records\n";
print {$LOG} "Mapped records: $mapped_records\n";
print {$LOG} "Unmapped records: $unmapped_records\n";
print {$LOG} "Phage-filtered records: $phage_filtered_records\n";
print {$LOG} "Groupnames in table: ", scalar(@groups), "\n";
print {$LOG} "Samples positive for at least one Groupname: ",
    scalar(grep { scalar(keys %{$sample_positive_groups{$_} || {}}) > 0 } @samples), "\n";
print {$LOG} "Samples with multiple parsed files: ",
    scalar(grep { scalar(keys %{$sample_source_file_count{$_}}) > 1 } @samples), "\n";
print {$LOG} "Output directory: $opt{outdir}\n";
close $LOG;

print "[INFO] Summary completed successfully\n";
print "[INFO] Parsed files       : ", scalar(@parsed_files), "\n";
print "[INFO] Unique samples     : $total_samples\n";
print "[INFO] Parsed records     : $parsed_records\n";
print "[INFO] Mapped records     : $mapped_records\n";
print "[INFO] Unmapped records   : $unmapped_records\n";
print "[INFO] Phage-filtered     : $phage_filtered_records\n";
print "[INFO] Groupnames         : ", scalar(@groups), "\n";
print "[INFO] Output directory   : $opt{outdir}\n";
exit 0;

# ============================================================
# Subroutines
# ============================================================

sub read_group_table {
    my ($file) = @_;
    open my $IN, '<', $file
        or die "[ERROR] Cannot open Groupname table $file: $!\n";

    my ($header_line, @first_data);
    while (my $line = <$IN>) {
        chomp $line;
        $line =~ s/\r$//;
        next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
        $header_line = $line;
        last;
    }
    defined $header_line
        or die "[ERROR] Groupname table is empty: $file\n";

    my @header = split /\t/, $header_line, -1;
    my %column;
    for my $i (0 .. $#header) {
        my $name = lc($header[$i]);
        $name =~ s/^\s+|\s+$//g;
        $column{$name} = $i;
    }

    my $accession_idx;
    for my $candidate (qw(accession reference_accession refseq_accession)) {
        if (exists $column{$candidate}) {
            $accession_idx = $column{$candidate};
            last;
        }
    }

    my $group_idx = $column{groupname};

    my @rows;
    if (defined $accession_idx && defined $group_idx) {
        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/\r$//;
            next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
            push @rows, [split(/\t/, $line, -1)];
        }
    } else {
        # Headerless two-column legacy table: accession<TAB>Groupname.
        push @rows, \@header;
        while (my $line = <$IN>) {
            chomp $line;
            $line =~ s/\r$//;
            next if $line =~ /^\s*$/ || $line =~ /^\s*#/;
            push @rows, [split(/\t/, $line, -1)];
        }
        $accession_idx = 0;
        $group_idx = 1;
    }
    close $IN;

    my (%map, %groups);
    my $line_number = 1;
    for my $row (@rows) {
        $line_number++;
        next unless defined $row->[$accession_idx] && defined $row->[$group_idx];

        my $accession = normalize_accession($row->[$accession_idx]);
        my $group = $row->[$group_idx];
        $group =~ s/^\s+|\s+$//g if defined $group;

        next unless length $accession;
        next unless defined $group && length $group;
        next if uc($group) eq 'NA' || uc($group) eq 'N/A';

        if (exists $map{$accession} && $map{$accession} ne $group) {
            die "[ERROR] Conflicting Groupname assignments for $accession in $file\n";
        }

        $map{$accession} = $group;
        $groups{$group} = 1;
    }

    die "[ERROR] No valid ACCESSION-to-Groupname rows found in $file\n"
        unless keys %map;

    return (\%map, \%groups);
}

sub collect_parsed_files {
    my ($path) = @_;
    return ($path) if -f $path;

    my @files;
    find(
        {
            wanted => sub {
                return unless -f $_;
                return unless $_ =~ /_blast_out_parsed\.txt\z/;
                push @files, $File::Find::name;
            },
            no_chdir => 1,
        },
        $path
    );
    return sort @files;
}

sub check_incomplete_sample_directories {
    my ($root, $parsed_ref) = @_;
    my %parsed_parent;
    for my $file (@{$parsed_ref}) {
        $parsed_parent{dirname_for_file($file)} = 1;
    }

    opendir my $DH, $root or die "[ERROR] Cannot open input directory $root: $!\n";
    my @children = grep { $_ ne '.' && $_ ne '..' && -d "$root/$_" } readdir $DH;
    closedir $DH;

    my @missing;
    for my $child (@children) {
        my $dir = "$root/$child";
        my @shells = (
            glob("$dir/*.run.sh"),
            glob("$dir/*.run_viral_scan.sh")
        );
        next unless @shells;
        push @missing, $child unless $parsed_parent{$dir};
    }

    if (@missing) {
        die "[ERROR] Parsed result is missing for sample directories: "
            . join(', ', sort @missing)
            . "\nRun the corresponding sample shells before Summary.\n";
    }
}

sub dirname_for_file {
    my ($file) = @_;
    require File::Basename;
    return File::Basename::dirname($file);
}

sub read_parsed_header {
    my ($file) = @_;
    open my $IN, '<', $file
        or die "[ERROR] Cannot open parsed file $file: $!\n";
    my $header = <$IN>;
    close $IN;
    defined $header
        or die "[ERROR] Parsed file has no header: $file\n";

    chomp $header;
    $header =~ s/\r$//;
    my @header = split /\t/, $header, -1;
    my %column;

    for my $i (0 .. $#header) {
        my $name = lc($header[$i]);
        $name =~ s/^\s+|\s+$//g;
        $column{$name} = $i;
    }

    for my $required (qw(sample query subject)) {
        die "[ERROR] Required column '$required' not found in $file\n"
            unless exists $column{$required};
    }

    my $sample = basename($file);
    $sample =~ s/_blast_out_parsed\.txt\z//;
    return (\%column, $sample);
}

sub field_value {
    my ($fields, $index) = @_;
    return undef unless defined $index && $index <= $#{$fields};
    return $fields->[$index];
}

sub normalize_accession {
    my ($accession) = @_;
    return '' unless defined $accession;
    $accession =~ s/^\s+|\s+$//g;
    $accession =~ s/\.\d+\z//;
    return $accession;
}

sub quantile {
    my ($values, $probability) = @_;
    my @v = @{$values};
    return 'NA' unless @v;
    return $v[0] if @v == 1;

    my $position = ($#v) * $probability;
    my $lower = int($position);
    my $upper = $lower + 1;
    return $v[$lower] if $upper > $#v;

    my $fraction = $position - $lower;
    return sprintf('%.2f', $v[$lower] + ($v[$upper] - $v[$lower]) * $fraction);
}

sub normalize_output_path {
    my ($path) = @_;
    return abs_path($path) if -d $path;
    return $path if $path =~ m{\A/};
    return getcwd() . "/$path";
}

sub usage {
    my ($exit_code) = @_;
    print STDERR <<'USAGE';
Usage:
  perl bin/summarizeGroupnameResults.pl \
    --input ViruScanNIPT_result \
    --group-table database/annotation/viralGroupname.tsv \
    --outdir ViruScanNIPT_result/Summary \
    --positive-hit-cutoff 2

Required:
  --input             One parsed file or a directory containing parsed files.
  --group-table       TSV with ACCESSION (or Reference_accession) and Groupname.
  --outdir            Summary output directory.

Optional:
  --positive-hit-cutoff  Minimum distinct query reads per sample per Groupname.
                         Default: 2.
  --help                 Show this help.

The input directory is searched recursively for *_blast_out_parsed.txt.
Accession version suffixes are ignored during mapping.
Each query is counted once per sample per Groupname.
When sample run shells are present but parsed results are missing, the script
stops rather than silently reducing the cohort denominator.
USAGE
    exit($exit_code);
}
