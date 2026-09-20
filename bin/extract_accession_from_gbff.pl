#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my ($input, $output, $help);

GetOptions(
    'input|i=s'  => \$input,
    'output|o=s' => \$output,
    'help|h'     => \$help,
) or die usage();

die usage() if $help;
die "ERROR: missing --input\n\n"  . usage() unless $input;
die "ERROR: missing --output\n\n" . usage() unless $output;

die "ERROR: input file not found: $input\n" unless -e $input;

print STDERR "[INFO] Input  : $input\n";
print STDERR "[INFO] Output : $output\n";

# ------------------------------------------------------------
# 打开输入
# ------------------------------------------------------------
my $IN;

if ($input =~ /\.gz$/) {
    $IN = IO::Uncompress::Gunzip->new($input)
        or die "ERROR: cannot open gzip file $input: $GunzipError\n";
}
else {
    open $IN, '<', $input
        or die "ERROR: cannot open $input: $!\n";
}

open my $OUT, '>', $output
    or die "ERROR: cannot write $output: $!\n";

my %seen;
my $count = 0;

# ------------------------------------------------------------
# 提取 ACCESSION
#
# GenBank 格式：
#
# ACCESSION   NC_029549
#
# 有些记录理论上可能有多个 accession：
#
# ACCESSION   ABC123 ABC124
#
# 因此全部提取。
# ------------------------------------------------------------
while (my $line = <$IN>) {

    next unless $line =~ /^ACCESSION\s+(.+)/;

    my $accession_text = $1;

    $accession_text =~ s/^\s+//;
    $accession_text =~ s/\s+$//;

    my @accessions = split /\s+/, $accession_text;

    for my $acc (@accessions) {

        next unless defined $acc && $acc ne '';

        # 理论上 ACCESSION 行通常没有版本号，
        # 这里仍做兼容处理
        $acc =~ s/\.\d+$//;

        next if $seen{$acc}++;

        print $OUT "$acc\n";
        $count++;
    }
}

close $IN;
close $OUT;

print STDERR "[INFO] Unique ACCESSIONs: $count\n";
print STDERR "[INFO] Done.\n";

exit 0;


sub usage {

    return <<'USAGE';

Usage:

perl extract_accession_from_gbff.pl \
    -i viral.1.genomic.gbff.gz \
    -o viral_accession.list

Options:

  -i, --input
      Input GenBank GBFF or GBFF.GZ file

  -o, --output
      Output ACCESSION list

  -h, --help
      Show help

USAGE
}

