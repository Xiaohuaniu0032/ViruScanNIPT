#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename qw(basename);
use File::Temp qw(tempfile);

# =============================================================================
# Parameters
# =============================================================================

my $input;
my $output;
my $sample;

my $taxonomy_dir = '/ehpcdata/fulongfei/database/taxonomy';

my ($acc2tax_file, $nodes_file, $names_file, $merged_file);

my $help = 0;

GetOptions(
    'input|i=s'         => \$input,
    'output|o=s'        => \$output,
    'sample|s=s'        => \$sample,
    'taxonomy-dir|t=s'  => \$taxonomy_dir,

    'accession2taxid=s' => \$acc2tax_file,
    'nodes=s'           => \$nodes_file,
    'names=s'           => \$names_file,
    'merged=s'          => \$merged_file,

    'help|h'            => \$help,
) or usage(1);

usage(0) if $help;
usage(1) unless defined $input && defined $output;

# =============================================================================
# Taxonomy files
# =============================================================================

$acc2tax_file ||= "$taxonomy_dir/nucl_gb.accession2taxid.gz";
$nodes_file   ||= "$taxonomy_dir/nodes.dmp";
$names_file   ||= "$taxonomy_dir/names.dmp";
$merged_file  ||= "$taxonomy_dir/merged.dmp";

for my $f ($input, $acc2tax_file, $nodes_file, $names_file) {
    die "[ERROR] File not found: $f\n" unless -e $f;
}

# =============================================================================
# Infer sample name
# =============================================================================

if (!defined $sample || $sample eq '') {
    $sample = basename($input);

    $sample =~ s/\.blast\.out$//;
    $sample =~ s/\.out$//;
    $sample =~ s/\.txt$//;
}

print STDERR "[INFO] Input              : $input\n";
print STDERR "[INFO] Output             : $output\n";
print STDERR "[INFO] Sample             : $sample\n";
print STDERR "[INFO] accession2taxid    : $acc2tax_file\n";
print STDERR "[INFO] nodes.dmp          : $nodes_file\n";
print STDERR "[INFO] names.dmp          : $names_file\n";
print STDERR "[INFO] Final rank         : genus\n";
print STDERR "[INFO] Best-hit priority  : bit_score > E-value > identity > alignment_length\n";

# =============================================================================
# Pass 1
#
# Parse BLAST output.
#
# For each Query:
#   1. Collapse multiple HSPs from the same subject
#   2. Select best subject(s)
#   3. Keep all completely tied best subjects
#
# Store results in temporary file.
# =============================================================================

my ($best_fh, $best_tmp) = tempfile(
    'blast_best_XXXXXX',
    TMPDIR => 1,
    UNLINK => 1
);

my %wanted_acc;

my $query_count          = 0;
my $query_with_hit_count = 0;

open my $IN, '<', $input
    or die "[ERROR] Cannot open $input: $!\n";

my $current_query;
my @current_hits;

while (my $line = <$IN>) {

    chomp $line;
    $line =~ s/\r$//;

    # -------------------------------------------------------------------------
    # BLAST outfmt 7
    # -------------------------------------------------------------------------

    if ($line =~ /^#\s*Query:\s*(\S+)/) {

        if (defined $current_query) {

            flush_query(
                $current_query,
                \@current_hits,
                $best_fh,
                \%wanted_acc,
                \$query_count,
                \$query_with_hit_count
            );
        }

        $current_query = $1;
        @current_hits  = ();

        next;
    }

    next if $line =~ /^#/;
    next if $line =~ /^\s*$/;

    # -------------------------------------------------------------------------
    # Data line
    # -------------------------------------------------------------------------

    my @f = split /\t/, $line;

    # fallback for whitespace-delimited input
    @f = split /\s+/, $line if @f < 12;

    die "[ERROR] BLAST line has fewer than 12 columns:\n$line\n"
        if @f < 12;

    my (
        $q,
        $subject,
        $pident,
        $aln_len,
        $mismatch,
        $gapopen,
        $qstart,
        $qend,
        $sstart,
        $send,
        $evalue,
        $bitscore
    ) = @f[0..11];

    # -------------------------------------------------------------------------
    # Also support BLAST outfmt 6
    # -------------------------------------------------------------------------

    if (!defined $current_query) {

        $current_query = $q;

    }
    elsif ($q ne $current_query) {

        flush_query(
            $current_query,
            \@current_hits,
            $best_fh,
            \%wanted_acc,
            \$query_count,
            \$query_with_hit_count
        );

        $current_query = $q;
        @current_hits  = ();
    }

    push @current_hits, {

        subject  => $subject,
        pident   => 0 + $pident,
        aln_len  => 0 + $aln_len,
        evalue   => 0 + $evalue,
        bitscore => 0 + $bitscore,
    };
}

# Last Query
if (defined $current_query) {

    flush_query(
        $current_query,
        \@current_hits,
        $best_fh,
        \%wanted_acc,
        \$query_count,
        \$query_with_hit_count
    );
}

close $IN;
close $best_fh;

my $wanted_count = scalar keys %wanted_acc;

print STDERR "[INFO] Queries parsed     : $query_count\n";
print STDERR "[INFO] Queries with hits  : $query_with_hit_count\n";
print STDERR "[INFO] Best-hit accessions: $wanted_count\n";

# =============================================================================
# accession -> TaxID
#
# nucl_gb.accession2taxid.gz is ~2.5 GB.
#
# DO NOT load the whole file into RAM.
#
# Stream through it and only keep accessions actually found in the BLAST
# best-hit results.
# =============================================================================

my %acc2tax;

if ($wanted_count > 0) {

    print STDERR "[INFO] Scanning accession2taxid ...\n";

    my $ACC = open_text_or_gzip($acc2tax_file);

    my $mapped = 0;

    while (my $line = <$ACC>) {

        chomp $line;
        $line =~ s/\r$//;

        next if $line =~ /^accession\taccession\.version\ttaxid\tgi/;

        my (
            $acc,
            $accver,
            $taxid
        ) = split /\t/, $line, 4;

        next unless defined $taxid;

        # accession.version
        if (
            exists $wanted_acc{$accver}
            &&
            !exists $acc2tax{$accver}
        ) {

            $acc2tax{$accver} = $taxid;
            $mapped++;
        }

        # accession without version
        if (
            exists $wanted_acc{$acc}
            &&
            !exists $acc2tax{$acc}
        ) {

            $acc2tax{$acc} = $taxid;
            $mapped++;
        }

        # Stop early if all requested accessions have been found
        last if $mapped >= $wanted_count;
    }

    close $ACC;

    print STDERR "[INFO] Accessions mapped  : "
        . scalar(keys %acc2tax)
        . "/$wanted_count\n";
}

# =============================================================================
# Load merged.dmp
#
# Used when old TaxIDs have been merged into newer TaxIDs.
# =============================================================================

my %merged;

if (
    defined $merged_file
    &&
    -e $merged_file
) {

    open my $MERGED, '<', $merged_file
        or die "[ERROR] Cannot open $merged_file: $!\n";

    while (my $line = <$MERGED>) {

        my @f = split /\t\|\t/, $line;

        next unless @f >= 2;

        my $old = trim_tax_field($f[0]);
        my $new = trim_tax_field($f[1]);

        $merged{$old} = $new
            if $old ne '' && $new ne '';
    }

    close $MERGED;
}

# =============================================================================
# Load nodes.dmp
#
# For memory efficiency only store:
#
#   TaxID -> parent_taxid + whether rank == genus
#
# We do NOT need all taxonomy rank strings.
# =============================================================================

print STDERR "[INFO] Loading nodes.dmp ...\n";

my %node;

open my $NODES, '<', $nodes_file
    or die "[ERROR] Cannot open $nodes_file: $!\n";

while (my $line = <$NODES>) {

    my @f = split /\t\|\t/, $line;

    next unless @f >= 3;

    my $taxid  = trim_tax_field($f[0]);
    my $parent = trim_tax_field($f[1]);
    my $rank   = trim_tax_field($f[2]);

    $node{$taxid} =
        $parent
        . "\t"
        . ($rank eq 'genus' ? 1 : 0);
}

close $NODES;

print STDERR "[INFO] Taxonomy nodes      : "
    . scalar(keys %node)
    . "\n";

# =============================================================================
# Convert TaxID -> Genus TaxID
# =============================================================================

my %tax2genus;
my %genus_cache;
my %wanted_genus;

for my $acc (keys %acc2tax) {

    my $taxid = $acc2tax{$acc};

    my $genus_taxid = find_genus_taxid(
        $taxid,
        \%node,
        \%merged,
        \%genus_cache
    );

    if (defined $genus_taxid) {

        $tax2genus{$taxid} = $genus_taxid;

        $wanted_genus{$genus_taxid} = 1;
    }
}

# =============================================================================
# Load scientific names only for genus TaxIDs we actually need
# =============================================================================

my %genus_name;

if (%wanted_genus) {

    print STDERR "[INFO] Loading required genus names ...\n";

    open my $NAMES, '<', $names_file
        or die "[ERROR] Cannot open $names_file: $!\n";

    while (my $line = <$NAMES>) {

        my @f = split /\t\|\t/, $line;

        next unless @f >= 4;

        my $taxid = trim_tax_field($f[0]);

        next unless exists $wanted_genus{$taxid};

        my $name       = trim_tax_field($f[1]);
        my $name_class = trim_tax_field($f[3]);

        if ($name_class eq 'scientific name') {

            $genus_name{$taxid} = $name;
        }
    }

    close $NAMES;
}

# =============================================================================
# Final output
# =============================================================================

open my $OUT, '>', $output
    or die "[ERROR] Cannot write $output: $!\n";

print $OUT join(
    "\t",
    qw(
        Sample
        Query
        Genus
        Genus_taxid
        Assignment
        Best_hit_count
        Best_subjects
        Best_taxids
        Candidate_genera
        Pident
        Alignment_length
        Evalue
        Bit_score
    )
), "\n";

open my $BEST, '<', $best_tmp
    or die "[ERROR] Cannot reopen temporary file $best_tmp: $!\n";

my %stats = (

    No_hit                => 0,
    Unique_best           => 0,
    Multi_best_same_genus => 0,
    Ambiguous_genus       => 0,
    Unclassified_genus    => 0,
);

while (my $line = <$BEST>) {

    chomp $line;

    my (
        $query,
        $best_count,
        $subjects_str,
        $pident,
        $aln_len,
        $evalue,
        $bitscore
    ) = split /\t/, $line, -1;

    # =========================================================================
    # No BLAST hit
    # =========================================================================

    if ($best_count == 0) {

        print $OUT join(
            "\t",

            $sample,
            $query,

            'No_hit',
            'NA',
            'No_hit',

            0,

            'NA',
            'NA',
            'NA',

            'NA',
            'NA',
            'NA',
            'NA'

        ), "\n";

        $stats{No_hit}++;

        next;
    }

    # =========================================================================
    # Best hits
    # =========================================================================

    my @subjects = split /,/, $subjects_str;

    my @best_taxids;

    my %candidate_genus_ids;

    my $unresolved = 0;

    for my $subject (@subjects) {

        my $taxid = $acc2tax{$subject};

        # ---------------------------------------------------------------------
        # accession -> TaxID failed
        # ---------------------------------------------------------------------

        if (!defined $taxid) {

            push @best_taxids, 'NA';

            $unresolved++;

            next;
        }

        push @best_taxids, $taxid;

        # ---------------------------------------------------------------------
        # TaxID -> genus
        # ---------------------------------------------------------------------

        my $genus_taxid = $tax2genus{$taxid};

        if (defined $genus_taxid) {

            $candidate_genus_ids{$genus_taxid} = 1;

        }
        else {

            $unresolved++;
        }
    }

    my @genus_ids =
        sort { $a <=> $b }
        keys %candidate_genus_ids;

    my @candidate_names = map {

        exists $genus_name{$_}
            ? $genus_name{$_}
            : "TaxID_$_"

    } @genus_ids;

    my (
        $final_genus,
        $final_genus_taxid,
        $assignment
    );

    # =========================================================================
    # Some best hits cannot be resolved to genus
    #
    # Conservative approach:
    # do NOT force assignment even if the other hits belong to one genus.
    # =========================================================================

    if ($unresolved > 0) {

        $final_genus       = 'Unclassified_genus';

        $final_genus_taxid = 'NA';

        $assignment        = 'Unclassified_genus';
    }

    # =========================================================================
    # All best hits belong to exactly one genus
    # =========================================================================

    elsif (@genus_ids == 1) {

        $final_genus_taxid = $genus_ids[0];

        $final_genus =
            $genus_name{$genus_ids[0]}
            // "TaxID_$genus_ids[0]";

        if ($best_count == 1) {

            $assignment = 'Unique_best';

        }
        else {

            $assignment = 'Multi_best_same_genus';
        }
    }

    # =========================================================================
    # Best hits span >1 genus
    # =========================================================================

    elsif (@genus_ids > 1) {

        $final_genus       = 'Ambiguous_genus';

        $final_genus_taxid = 'NA';

        $assignment        = 'Ambiguous_genus';
    }

    # =========================================================================
    # Nothing can be assigned
    # =========================================================================

    else {

        $final_genus       = 'Unclassified_genus';

        $final_genus_taxid = 'NA';

        $assignment        = 'Unclassified_genus';
    }

    my $candidate_genera =
        @candidate_names
        ? join(',', @candidate_names)
        : 'NA';

    print $OUT join(
        "\t",

        $sample,
        $query,

        $final_genus,
        $final_genus_taxid,

        $assignment,

        $best_count,

        $subjects_str,

        join(',', @best_taxids),

        $candidate_genera,

        $pident,
        $aln_len,
        $evalue,
        $bitscore

    ), "\n";

    $stats{$assignment}++;
}

close $BEST;
close $OUT;

# =============================================================================
# Summary
# =============================================================================

print STDERR "=" x 72, "\n";

print STDERR "[INFO] Output file                : $output\n";

print STDERR "[INFO] No_hit                     : "
    . $stats{No_hit}
    . "\n";

print STDERR "[INFO] Unique_best                : "
    . $stats{Unique_best}
    . "\n";

print STDERR "[INFO] Multi_best_same_genus      : "
    . $stats{Multi_best_same_genus}
    . "\n";

print STDERR "[INFO] Ambiguous_genus            : "
    . $stats{Ambiguous_genus}
    . "\n";

print STDERR "[INFO] Unclassified_genus         : "
    . $stats{Unclassified_genus}
    . "\n";

print STDERR "[INFO] Done.\n";

exit 0;

# =============================================================================
# Functions
# =============================================================================

sub flush_query {

    my (
        $query,
        $hits_ref,
        $fh,
        $wanted_ref,
        $query_count_ref,
        $with_hit_ref
    ) = @_;

    $$query_count_ref++;

    # =========================================================================
    # No hit
    # =========================================================================

    if (!@$hits_ref) {

        print $fh join(
            "\t",

            $query,
            0,

            'NA',
            'NA',
            'NA',
            'NA',
            'NA'

        ), "\n";

        return;
    }

    $$with_hit_ref++;

    # =========================================================================
    # STEP 1
    #
    # Collapse multiple HSPs from the same subject.
    #
    # If one read has multiple alignments to different locations of the same
    # viral genome, the subject should still count only once.
    # =========================================================================

    my %subject_best;

    for my $h (@$hits_ref) {

        my $subject = $h->{subject};

        if (
            !exists $subject_best{$subject}
            ||
            hit_is_better(
                $h,
                $subject_best{$subject}
            )
        ) {

            $subject_best{$subject} = $h;
        }
    }

    # =========================================================================
    # STEP 2
    #
    # Find global best subject(s).
    # =========================================================================

    my @best;

    for my $subject (sort keys %subject_best) {

        my $h = $subject_best{$subject};

        if (
            !@best
            ||
            hit_is_better(
                $h,
                $best[0]
            )
        ) {

            @best = ($h);

        }
        elsif (
            hits_are_tied(
                $h,
                $best[0]
            )
        ) {

            push @best, $h;
        }
    }

    # =========================================================================
    # Collect subjects
    # =========================================================================

    my @subjects =
        sort
        map {
            $_->{subject}
        } @best;

    $wanted_ref->{$_} = 1
        for @subjects;

    my $b = $best[0];

    print $fh join(
        "\t",

        $query,

        scalar(@subjects),

        join(',', @subjects),

        $b->{pident},
        $b->{aln_len},
        $b->{evalue},
        $b->{bitscore}

    ), "\n";
}

# =============================================================================
# Compare two hits
#
# Priority:
#
#   1. highest bit score
#   2. lowest E-value
#   3. highest identity
#   4. longest alignment
# =============================================================================

sub hit_is_better {

    my ($a, $b) = @_;

    # bit score
    return 1
        if $a->{bitscore} > $b->{bitscore};

    return 0
        if $a->{bitscore} < $b->{bitscore};

    # E-value
    return 1
        if $a->{evalue} < $b->{evalue};

    return 0
        if $a->{evalue} > $b->{evalue};

    # identity
    return 1
        if $a->{pident} > $b->{pident};

    return 0
        if $a->{pident} < $b->{pident};

    # alignment length
    return 1
        if $a->{aln_len} > $b->{aln_len};

    return 0;
}

# =============================================================================
# Check whether two hits are completely tied
# =============================================================================

sub hits_are_tied {

    my ($a, $b) = @_;

    return (

        $a->{bitscore}
            ==
        $b->{bitscore}

        &&

        $a->{evalue}
            ==
        $b->{evalue}

        &&

        $a->{pident}
            ==
        $b->{pident}

        &&

        $a->{aln_len}
            ==
        $b->{aln_len}
    );
}

# =============================================================================
# Find genus TaxID
#
# Walk:
#
# species -> genus -> ...
#
# until rank == genus.
# =============================================================================

sub find_genus_taxid {

    my (
        $start_taxid,
        $node_ref,
        $merged_ref,
        $cache_ref
    ) = @_;

    my $taxid = $start_taxid;

    my @path;

    my %seen;

    # =========================================================================
    # Resolve merged/obsolete TaxID
    # =========================================================================

    while (
        !exists $node_ref->{$taxid}
        &&
        exists $merged_ref->{$taxid}
    ) {

        $taxid = $merged_ref->{$taxid};
    }

    # =========================================================================
    # Walk toward root
    # =========================================================================

    while (
        defined $taxid
        &&
        $taxid ne ''
    ) {

        # cached result
        if (exists $cache_ref->{$taxid}) {

            my $g = $cache_ref->{$taxid};

            $cache_ref->{$_} = $g
                for @path;

            return $g eq 'NA'
                ? undef
                : $g;
        }

        last if $seen{$taxid}++;

        push @path, $taxid;

        my $node = $node_ref->{$taxid};

        last unless defined $node;

        my (
            $parent,
            $is_genus
        ) = split /\t/, $node, 2;

        # ---------------------------------------------------------------------
        # Found genus
        # ---------------------------------------------------------------------

        if ($is_genus) {

            $cache_ref->{$_} = $taxid
                for @path;

            return $taxid;
        }

        # root
        last
            if !defined $parent
            ||
            $parent eq ''
            ||
            $parent eq $taxid;

        $taxid = $parent;
    }

    # =========================================================================
    # No genus found
    # =========================================================================

    $cache_ref->{$_} = 'NA'
        for @path;

    return undef;
}

# =============================================================================
# Clean NCBI taxonomy field
# =============================================================================

sub trim_tax_field {

    my ($x) = @_;

    return ''
        unless defined $x;

    $x =~ s/^\s+//;
    $x =~ s/\s+$//;
    $x =~ s/\t?\|\s*$//;

    return $x;
}

# =============================================================================
# Open plain text or .gz
# =============================================================================

sub open_text_or_gzip {

    my ($file) = @_;

    my $fh;

    if ($file =~ /\.gz$/i) {

        open(
            $fh,
            '-|',
            'gzip',
            '-dc',
            $file
        )
        or die
            "[ERROR] Cannot open gzip stream for $file: $!\n";

    }
    else {

        open(
            $fh,
            '<',
            $file
        )
        or die
            "[ERROR] Cannot open $file: $!\n";
    }

    return $fh;
}

# =============================================================================
# Help
# =============================================================================

sub usage {

    my ($exit_code) = @_;

    print STDERR <<'USAGE';

Usage:

  perl parse_blast_genus.pl \
      -i SAMPLE.blast.out \
      -o SAMPLE.blast.genus.tsv \
      -t /ehpcdata/fulongfei/database/taxonomy


Required:

  -i, --input FILE

      BLAST outfmt 7 or outfmt 6 file.

      First 12 columns must be:

      qaccver
      saccver
      pident
      length
      mismatch
      gapopen
      qstart
      qend
      sstart
      send
      evalue
      bitscore


  -o, --output FILE

      Output TSV file.


Optional:

  -s, --sample STRING

      Sample name.

      Default:
      infer from input filename.


  -t, --taxonomy-dir DIR

      NCBI taxonomy directory.

      Default:

      /ehpcdata/fulongfei/database/taxonomy


  --accession2taxid FILE

      Override:

      nucl_gb.accession2taxid.gz


  --nodes FILE

      Override:

      nodes.dmp


  --names FILE

      Override:

      names.dmp


  --merged FILE

      Override:

      merged.dmp


  -h, --help

      Show help.


Best-hit rule:

  1. Highest bit score

  2. Lowest E-value

  3. Highest percent identity

  4. Longest alignment length

  5. If still tied:
     keep ALL tied best subjects


Genus assignment:

  One best subject
      ->
  Unique_best


  Multiple best subjects
  but all belong to same genus
      ->
  Multi_best_same_genus


  Multiple best subjects
  belonging to >1 genus
      ->
  Ambiguous_genus


  TaxID or genus cannot be resolved
      ->
  Unclassified_genus


  Query has no BLAST hit
      ->
  No_hit

USAGE

    exit $exit_code;
}

