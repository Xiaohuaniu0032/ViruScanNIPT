#!/usr/bin/env perl
use strict;
use warnings;
use Fcntl qw(O_WRONLY O_CREAT O_EXCL);

# Usage:
#   perl add_genus_to_metadata.pl metadata.tsv nodes.dmp names.dmp output.tsv
# Input is TAB-delimited and must have a column named taxon.
# Preserve all input columns and append Genus. Only taxon is used for lookup.
# Genus names come from names.dmp records of class "scientific name".
# Unresolved values become NA; an existing output file is never overwritten.
# No non-core Perl modules are required. This script does not access the web.

@ARGV == 4 or die "Usage: perl $0 metadata.tsv nodes.dmp names.dmp output.tsv\n";
my ($input, $nodes_file, $names_file, $output) = @ARGV;
for my $file ($input, $nodes_file, $names_file) {
    -f $file && -r $file or die "Cannot read file: $file\n";
}
(-e $output || -l $output)
    and die "Output already exists; choose a new filename: $output\n";

# Read and validate the metadata table before creating the output file.
open my $in, '<', $input or die "Cannot open $input: $!\n";
my $header = <$in>;
defined $header or die "Input table is empty: $input\n";
$header =~ s/\r?\n\z//;
$header =~ s/^\xEF\xBB\xBF//;
my @columns = split /\t/, $header, -1;
my @taxon_columns = grep { $columns[$_] eq 'taxon' } 0 .. $#columns;
@taxon_columns == 1
    or die "Input must be TAB-delimited with exactly one taxon column\n";
!grep { $_ eq 'Genus' } @columns
    or die "Input already contains a Genus column\n";
my $taxon_index = $taxon_columns[0];
my (@rows, %wanted);
my $line_number = 1;
while (my $line = <$in>) {
    ++$line_number;
    $line =~ s/\r?\n\z//;
    next unless $line =~ /\S/;
    my @fields = split /\t/, $line, -1;
    @fields == @columns
        or die "Wrong column count at $input line $line_number\n";
    my $taxid = trim($fields[$taxon_index]);
    push @rows, [$line, $taxid];
    $wanted{$taxid} = 1;
}
close $in or die "Cannot close $input: $!\n";

# nodes.dmp: TaxID, parent TaxID, rank, ...
# Keep only parent links and whether each node has the exact rank genus.
my (%parent, %is_genus, %lineage_cache);
print STDERR "[INFO] Reading nodes.dmp\n";
open my $nodes, '<', $nodes_file or die "Cannot open $nodes_file: $!\n";
while (my $line = <$nodes>) {
    next unless $line =~ /\S/;
    my ($id, $parent_id, $rank) = split /\t\|\t?/, $line, 4;
    defined $rank or die "Malformed nodes.dmp record at line $.\n";
    ($id, $parent_id, $rank) = map { trim($_) } ($id, $parent_id, $rank);
    $id =~ /\A[0-9]+\z/ && $parent_id =~ /\A[0-9]+\z/
        or die "Invalid TaxID in nodes.dmp at line $.\n";
    $parent{$id} = $parent_id;
    $is_genus{$id} = 1 if $rank eq 'genus';
}
close $nodes or die "Cannot close $nodes_file: $!\n";
%parent or die "No taxonomy nodes found in $nodes_file\n";

# Resolve each distinct input TaxID once. Never infer taxonomy from ORGANISM.
my (%result, %reason, %genus_for_taxon, %needed_genus);
for my $taxid (keys %wanted) {
    $result{$taxid} = 'NA';
    if (!length($taxid) || $taxid eq 'NA') {
        $reason{$taxid} = 'missing_taxon';
    } elsif ($taxid !~ /\A[1-9][0-9]*\z/) {
        $reason{$taxid} = 'invalid_taxon';
    } elsif (!exists $parent{$taxid}) {
        $reason{$taxid} = 'taxon_not_in_nodes';
    } else {
        my $genus_id = find_genus($taxid);
        if ($genus_id) {
            $genus_for_taxon{$taxid} = $genus_id;
            $needed_genus{$genus_id} = 1;
        } else {
            $reason{$taxid} = 'no_genus_or_incomplete_lineage';
        }
    }
}

# names.dmp: TaxID, name text, unique name, name class.
# Retain only scientific names of the genus nodes needed for this table.
my %genus_name;
print STDERR "[INFO] Reading names.dmp\n";
open my $names, '<', $names_file or die "Cannot open $names_file: $!\n";
while (my $line = <$names>) {
    next unless $line =~ /\S/;
    my ($id, $name, $unique, $class) = split /\t\|\t?/, $line, 5;
    defined $class or die "Malformed names.dmp record at line $.\n";
    $id = trim($id);
    next unless exists $needed_genus{$id};
    next unless trim($class) eq 'scientific name';
    $name = trim($name);
    $genus_name{$id} = $name if length $name;
}
close $names or die "Cannot close $names_file: $!\n";
for my $taxid (keys %genus_for_taxon) {
    my $genus_id = $genus_for_taxon{$taxid};
    if (exists $genus_name{$genus_id}) {
        $result{$taxid} = $genus_name{$genus_id};
    } else {
        $reason{$taxid} = 'genus_scientific_name_missing';
    }
}

# O_EXCL also prevents overwriting an input, symlink, or newly created output.
sysopen my $out, $output, O_WRONLY | O_CREAT | O_EXCL
    or die "Cannot create new output file $output: $!\n";
print {$out} $header, "\tGenus\n" or die "Cannot write $output: $!\n";
my $assigned = 0;
my %na_counts;
for my $row (@rows) {
    my ($line, $taxid) = @$row;
    print {$out} $line, "\t", $result{$taxid}, "\n"
        or die "Cannot write $output: $!\n";
    if (exists $reason{$taxid}) {
        ++$na_counts{$reason{$taxid}};
    } else {
        ++$assigned;
    }
}
close $out or die "Cannot finish writing $output: $!\n";
my $total = scalar @rows;
print STDERR "[DONE] Output: $output\n";
print STDERR "[DONE] Rows=$total; Genus_assigned=$assigned; NA=", $total - $assigned, "\n";
for my $reason (sort keys %na_counts) {
    print STDERR "[NA] $reason=$na_counts{$reason}\n";
}

# Follow parents, including the starting node itself, until rank=genus.
# Cache shared paths; roots, missing ancestors and cycles resolve to 0 (NA).
sub find_genus {
    my ($id) = @_;
    my (@path, %seen);
    my $answer = 0;
    while (exists $parent{$id}) {
        if (exists $lineage_cache{$id}) {
            $answer = $lineage_cache{$id};
            last;
        }
        last if $seen{$id}++;
        push @path, $id;
        if ($is_genus{$id}) {
            $answer = $id;
            last;
        }
        last if $parent{$id} eq $id;
        $id = $parent{$id};
    }
    $lineage_cache{$_} = $answer for @path;
    return $answer;
}

sub trim {
    my ($value) = @_;
    $value =~ s/\A\s+|\s+\z//g;
    return $value;
}
