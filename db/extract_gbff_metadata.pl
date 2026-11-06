#!/usr/bin/env perl
use strict;
use warnings;

# Usage: perl extract_gbff_metadata.pl input.gbff[.gz] > metadata.tsv
# Requires gzip only for .gz input; BioPerl is not required.
# ACCESSION is the primary accession, without the VERSION suffix.
# Output: ACCESSION, ORGANISM, mol_type, host, taxon (TAB-delimited).
# Metadata is taken from the primary biological source, selected by /focus or
# /transgenic, then by /organism matching the record's ORGANISM header.
# Among matching sources, prefer those spanning the entire LOCUS sequence.
# Wrapped ORGANISM names are resolved against complete source /organism values.
# taxon is the numeric ID from that source's /db_xref="taxon:ID" qualifier.

@ARGV == 1 or die "Usage: perl $0 input.gbff[.gz] > metadata.tsv\n";
my ($input) = @ARGV;
-f $input && -r $input or die "Cannot read input file: $input\n";

my $fh;
if ($input =~ /\.gz\z/i) {
    # List-form open passes the filename directly, without a shell.
    open $fh, '-|', 'gzip', '-dc', '--', $input
        or die "Cannot run gzip for $input: $!\n";
} else {
    open $fh, '<', $input or die "Cannot open $input: $!\n";
}

print join("\t", qw(ACCESSION ORGANISM mol_type host taxon)), "\n";

my $record = '';
my $in_sequence = 0;
while (my $line = <$fh>) {
    $line =~ s/\r\n\z/\n/;
    if ($line =~ m{^//\s*$}) {
        emit_record($record) if $record =~ /\S/;
        $record = '';
        $in_sequence = 0;
        next;
    }
    if ($line =~ /^(?:ORIGIN|CONTIG)\b/) {
        $in_sequence = 1;
        next;
    }
    # Keep annotations for one record; do not store sequence bases.
    $record .= $line unless $in_sequence;
}
close $fh or die "Input read/decompression failed for $input (status=$?, error=$!)\n";
emit_record($record) if $record =~ /\S/;
close STDOUT or die "Cannot finish writing output: $!\n";

sub emit_record {
    my ($text) = @_;
    my ($accession) = $text =~ /^ACCESSION[ \t]+(\S+)/m;
    defined $accession or die "A record has no ACCESSION field\n";
    my ($sequence_length) = $text =~ /^LOCUS[ \t]+\S+[ \t]+(\d+)[ \t]+(?:bp|aa)\b/m;
    my @sources;

    # Read only source features, stopping before the next feature.
    while ($text =~ /^ {5}source[ \t]+[^\n]*(?:\n {21}[^\n]*)*/mg) {
        my $source = $&;
        my %entry = (values => {});
        my ($location) = $source =~ /^ {5}source[ \t]+([^\n]*(?:\n {21}(?!\/)[^\n]*)*)/;
        $location =~ s/\s+//g;
        $location =~ s/^complement\((.*)\)$/$1/;
        $entry{covers_record} = defined $sequence_length
            && $location =~ /\A<?1\.\.>?(\d+)\z/ && $1 == $sequence_length;
        # Consume all qualifiers so quoted note text is not misread as host.
        # Parse flag qualifiers too; do not detect /focus inside a quoted note.
        while ($source =~ /^ {21}\/(\w+)(?:=(?:"((?:[^"]|"")*)"|([^\n]*)))?[ \t]*(?=\n|\z)/mg) {
            my ($key, $quoted, $unquoted) = ($1, $2, $3);
            if (($key eq 'focus' || $key eq 'transgenic')
                    && !defined $quoted && !defined $unquoted) {
                $entry{$key} = 1;
                next;
            }
            next unless $key eq 'mol_type' || $key eq 'host'
                || $key eq 'organism' || $key eq 'db_xref';
            next unless defined $quoted || defined $unquoted;
            my $value = defined $quoted ? $quoted : $unquoted;
            $value =~ s/""/"/g if defined $quoted;
            $value = clean_value($value);
            if ($key eq 'db_xref') {
                # Ignore other database cross-references, such as GeneID.
                push @{ $entry{values}{taxon} }, $1 if $value =~ /\Ataxon:\s*(\d+)\z/;
            } else {
                push @{ $entry{values}{$key} }, $value if length $value;
            }
        }
        push @sources, \%entry;
    }

    # ORGANISM names can wrap, and lineage lines use the same indentation.
    # Resolve the complete name before comparing it with source organisms.
    my $organism = read_organism_header($text, \@sources);

    # Do not mix viral metadata with source features for flanking host DNA.
    my @selected = select_sources($accession, $organism, \@sources);
    my %values;
    for my $source (@selected) {
        for my $key (qw(organism mol_type host taxon)) {
            push @{ $values{$key} }, @{ $source->{values}{$key} || [] };
        }
    }

    # The header is authoritative; /organism is a fallback if it is absent.
    $organism = defined $organism
        ? clean_value($organism) : joined_values($values{organism});
    $organism = 'NA' unless length $organism;
    my %seen_taxon;
    my @taxa = grep { !$seen_taxon{$_}++ } @{ $values{taxon} || [] };
    my $taxon = @taxa == 1 ? $taxa[0] : 'NA';
    # Conflicts now refer only to sources selected for the main organism.
    warn "[WARN] $accession has conflicting taxon IDs in selected primary sources ("
        . join(', ', @taxa) . "); taxon set to NA\n" if @taxa > 1;
    print join("\t", $accession, $organism,
        joined_values($values{mol_type}), joined_values($values{host}), $taxon), "\n";
}

sub read_organism_header {
    my ($text, $sources) = @_;
    my ($first, $continuation) = $text =~
        /^  ORGANISM[ \t]+([^\n]*)([\n] {12}[^\n]*(?:\n {12}[^\n]*)*)?/m;
    return undef unless defined $first;

    my %source_names;
    for my $source (@$sources) {
        $source_names{$_} = 1 for @{ $source->{values}{organism} || [] };
    }

    # Only accept an exact complete source name after joining whole header
    # lines. The longest exact match prevents a shorter name prefix winning.
    # This also keeps the following taxonomy lineage out of ORGANISM.
    my $prefix = '';
    my $matched;
    for my $part ($first, split /\n/, $continuation || '') {
        $part = clean_value($part);
        next unless length $part;
        $prefix .= (length($prefix) ? ' ' : '') . $part;
        $matched = $prefix if exists $source_names{$prefix};
    }
    return $matched if defined $matched;

    # For absent/inconsistent source names, read continuation lines until
    # the lineage begins. A lineage normally contains semicolons; a single
    # taxon lineage is a period-terminated line (for example, "Viruses.").
    my @parts = (clean_value($first));
    for my $part (split /\n/, $continuation || '') {
        $part = clean_value($part);
        next unless length $part;
        last if $part =~ /;/ || $part =~ /\.\z/;
        push @parts, $part;
    }
    return join(' ', @parts);
}

sub select_sources {
    my ($accession, $organism, $sources) = @_;
    return () unless @$sources;

    # INSDC /focus and /transgenic identify the main organism of the entry.
    my @marked = grep { $_->{focus} || $_->{transgenic} } @$sources;
    if (@marked == 1) {
        my $names = $marked[0]{values}{organism} || [];
        if (defined $organism && length $organism && @$names
                && !grep { $_ eq $organism } @$names) {
            warn "[WARN] $accession: primary source organism differs from ORGANISM; source fields set to NA\n";
            return ();
        }
        return @marked;
    }
    if (@marked > 1) {
        warn "[WARN] $accession: multiple primary source markers; source fields set to NA\n";
        return ();
    }

    if (defined $organism && length $organism) {
        my @matching = grep {
            my $source = $_;
            scalar grep { $_ eq $organism } @{ $source->{values}{organism} || [] }
        } @$sources;
        if (@matching) {
            my @full_length = grep { $_->{covers_record} } @matching;
            return @full_length if @full_length;
            return @matching;
        }
    }

    # A sole source is usable if a name is missing on either side.
    if (@$sources == 1 && (!defined $organism || !length $organism
            || !@{ $sources->[0]{values}{organism} || [] })) {
        return @$sources;
    }
    warn "[WARN] $accession: cannot identify primary source; source fields set to NA\n";
    return ();
}

sub clean_value {
    my ($value) = @_;
    $value =~ s/\s+/ /g;
    $value =~ s/^ +| +$//g;
    return $value;
}

sub joined_values {
    my ($items) = @_;
    return 'NA' unless $items && @$items;
    my %seen;
    return join('; ', grep { !$seen{$_}++ } @$items);
}


