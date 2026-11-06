#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

my ($accession_list, $taxid_file, $output, $help);

GetOptions(
    'accession|a=s' => \$accession_list,
    'taxid|t=s'     => \$taxid_file,
    'output|o=s'    => \$output,
    'help|h'        => \$help,
) or die usage();

die usage() if $help;

die "ERROR: missing --accession\n\n" . usage()
    unless $accession_list;

die "ERROR: missing --taxid\n\n" . usage()
    unless $taxid_file;

die "ERROR: missing --output\n\n" . usage()
    unless $output;

die "ERROR: file not found: $accession_list\n"
    unless -e $accession_list;

die "ERROR: file not found: $taxid_file\n"
    unless -e $taxid_file;

print STDERR "[INFO] ACCESSION list : $accession_list\n";
print STDERR "[INFO] TaxID database : $taxid_file\n";
print STDERR "[INFO] Output         : $output\n";

# ============================================================
# 1. 读取目标 ACCESSION
# ============================================================

my %wanted;

open my $LIST, '<', $accession_list
    or die "ERROR: cannot open $accession_list: $!\n";

while (<$LIST>) {

    chomp;
    s/\r$//;

    next if /^\s*$/;
    next if /^#/;

    my ($acc) = split /\s+/;

    next unless defined $acc;

    # 去除版本号
    # NC_029549.1 -> NC_029549
    $acc =~ s/\.\d+$//;

    $wanted{$acc} = 1;
}

close $LIST;

my $target_count = scalar keys %wanted;

print STDERR "[INFO] Target ACCESSIONs: $target_count\n";

# ============================================================
# 2. 扫描 accession2taxid
#
# 格式：
#
# accession   accession.version   taxid   gi
#
# ============================================================

my $IN;

if ($taxid_file =~ /\.gz$/) {

    $IN = IO::Uncompress::Gunzip->new($taxid_file)
        or die "ERROR: cannot open gzip file $taxid_file: $GunzipError\n";
}
else {

    open $IN, '<', $taxid_file
        or die "ERROR: cannot open $taxid_file: $!\n";
}

my %result;

my $line_count = 0;
my $found_count = 0;

while (my $line = <$IN>) {

    $line_count++;

    # 跳过表头
    next if $line_count == 1 && $line =~ /^accession\s+/;

    chomp $line;
    $line =~ s/\r$//;

    my ($accession, $accession_version, $taxid, $gi)
        = split /\t/, $line;

    next unless defined $accession;
    next unless defined $taxid;

    next unless exists $wanted{$accession};

    # 防止重复计数
    unless (exists $result{$accession}) {

        $result{$accession} = $taxid;
        $found_count++;

        # 所有 accession 都已经找到后，立即停止扫描
        last if $found_count == $target_count;
    }

    # 每扫描 1000 万行报告一次进度
    if ($line_count % 10_000_000 == 0) {

        print STDERR
            "[INFO] Scanned ",
            commify($line_count),
            " lines; found $found_count / $target_count\n";
    }
}

close $IN;

# ============================================================
# 3. 输出
# ============================================================

open my $OUT, '>', $output
    or die "ERROR: cannot write $output: $!\n";

print $OUT "ACCESSION\ttaxid\n";

# 按 ACCESSION 排序输出
for my $acc (sort keys %wanted) {

    my $taxid = exists $result{$acc}
              ? $result{$acc}
              : 'NA';

    print $OUT "$acc\t$taxid\n";
}

close $OUT;

# ============================================================
# 4. 统计
# ============================================================

my $missing_count = $target_count - $found_count;

print STDERR "\n";
print STDERR "=" x 60, "\n";
print STDERR "[INFO] Target ACCESSIONs : $target_count\n";
print STDERR "[INFO] Found             : $found_count\n";
print STDERR "[INFO] Missing           : $missing_count\n";
print STDERR "[INFO] Scanned lines     : ",
    commify($line_count), "\n";
print STDERR "[INFO] Output            : $output\n";
print STDERR "[INFO] Done.\n";

# 如果存在未找到 accession，额外输出列表
if ($missing_count > 0) {

    my $missing_file = $output . ".missing";

    open my $MISS, '>', $missing_file
        or die "ERROR: cannot write $missing_file: $!\n";

    for my $acc (sort keys %wanted) {

        next if exists $result{$acc};

        print $MISS "$acc\n";
    }

    close $MISS;

    print STDERR "[INFO] Missing list      : $missing_file\n";
}

exit 0;


sub commify {

    my $text = reverse $_[0];

    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;

    return scalar reverse $text;
}


sub usage {

    return <<'USAGE';

Usage:

perl extract_accession_taxid.pl \
    -a viral_accession.list \
    -t nucl_gb.accession2taxid.gz \
    -o viral_accession_taxid.tsv

Options:

  -a, --accession
      ACCESSION list generated from GBFF

  -t, --taxid
      NCBI accession2taxid file
      Example:
      nucl_gb.accession2taxid.gz

  -o, --output
      Output ACCESSION-taxid table

  -h, --help
      Show help

USAGE
}

