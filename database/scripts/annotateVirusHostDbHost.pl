
#!/usr/bin/env perl
use strict;
use warnings;

# --------------------------------------------------
# Usage:
# perl annotate_virushostdb_host.pl input.tsv virushostdb.tsv output.tsv
# --------------------------------------------------

my ($input_file, $db_file, $output_file) = @ARGV;

die "Usage: perl $0 input.tsv virushostdb.tsv output.tsv\n"
    unless defined $input_file
        && defined $db_file
        && defined $output_file;

# --------------------------------------------------
# Step 1: Load VirusHostDB
# --------------------------------------------------

open my $DB, '<', $db_file
    or die "Cannot open $db_file: $!\n";

my $db_header = <$DB>;
die "VirusHostDB file is empty!\n" unless defined $db_header;

chomp $db_header;
$db_header =~ s/\r$//;

my @db_header = split /\t/, $db_header, -1;

my %db_col;
for my $i (0 .. $#db_header) {
    $db_col{$db_header[$i]} = $i;
}

for my $name (qw(virus_tax_id host_name)) {
    die "Missing required column in VirusHostDB: $name\n"
        unless exists $db_col{$name};
}

my $virus_tax_idx = $db_col{'virus_tax_id'};
my $host_name_idx = $db_col{'host_name'};

# Store virus_tax_id -> unique host names
my %virus_host;
my %seen_host;

while (my $line = <$DB>) {

    chomp $line;
    $line =~ s/\r$//;

    next if $line eq '';

    my @fields = split /\t/, $line, -1;

    die "Invalid number of columns in VirusHostDB at line $.\n"
        if @fields <= $virus_tax_idx || @fields <= $host_name_idx;

    my $tax_id    = $fields[$virus_tax_idx];
    my $host_name = $fields[$host_name_idx];

    next if $tax_id eq '' || $host_name eq '';

    # Remove duplicate host names while preserving original order
    next if $seen_host{$tax_id}{$host_name}++;

    push @{$virus_host{$tax_id}}, $host_name;
}

close $DB;

# --------------------------------------------------
# Step 2: Read input and prepare output
# --------------------------------------------------

open my $IN, '<', $input_file
    or die "Cannot open $input_file: $!\n";

open my $OUT, '>', $output_file
    or die "Cannot write $output_file: $!\n";

my $header = <$IN>;
die "Input file is empty!\n" unless defined $header;

chomp $header;
$header =~ s/\r$//;

my @header = split /\t/, $header, -1;

my %col;
for my $i (0 .. $#header) {
    $col{$header[$i]} = $i;
}

die "Missing required column: taxon\n"
    unless exists $col{'taxon'};

die "Column VirusHostDB_host_name already exists!\n"
    if exists $col{'VirusHostDB_host_name'};

my $taxon_idx = $col{'taxon'};

# Output header
print $OUT join("\t", @header, 'VirusHostDB_host_name'), "\n";

# --------------------------------------------------
# Step 3: Match taxon and expand multiple host names
# --------------------------------------------------

my $total_count     = 0;
my $matched_count   = 0;
my $unmatched_count = 0;
my $output_count    = 0;

while (my $line = <$IN>) {

    chomp $line;
    $line =~ s/\r$//;

    next if $line eq '';

    $total_count++;

    my @fields = split /\t/, $line, -1;

    die "Invalid number of columns in input at line $.\n"
        if @fields <= $taxon_idx;

    my $tax_id = $fields[$taxon_idx];

    if (exists $virus_host{$tax_id}) {

        # Output one row for each unique host name
        for my $host_name (@{$virus_host{$tax_id}}) {

            print $OUT join("\t", @fields, $host_name), "\n";
            $output_count++;
        }

        $matched_count++;

    } else {

        # Retain unmatched records with NA
        print $OUT join("\t", @fields, 'NA'), "\n";

        $unmatched_count++;
        $output_count++;
    }
}

close $IN;
close $OUT;

# --------------------------------------------------
# Step 4: Statistics
# --------------------------------------------------

print "Annotation completed!\n";
print "Input records:        $total_count\n";
print "Matched records:      $matched_count\n";
print "Unmatched records:    $unmatched_count\n";
print "Output records:       $output_count\n";
print "Output file:          $output_file\n";

