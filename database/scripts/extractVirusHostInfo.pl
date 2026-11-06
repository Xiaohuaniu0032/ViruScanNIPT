#!/usr/bin/perl
use strict;
use warnings;

my $input  = "virushostdb.tsv";
my $output = "virushostdb.virus_host.tsv";

open my $IN, "<", $input or die "Cannot open $input: $!\n";
open my $OUT, ">", $output or die "Cannot open $output: $!\n";

# 读取表头
my $header = <$IN>;
defined $header or die "Input file is empty!\n";
chomp $header;
$header =~ s/\r$//;

my @header = split /\t/, $header, -1;

# 获取目标列索引
my %col;
for my $i (0 .. $#header) {
    $col{$header[$i]} = $i;
}

my @target = (
    "virus tax id",
    "virus name",
    "host tax id",
    "host name"
);

for my $name (@target) {
    exists $col{$name}
        or die "Column '$name' not found in input file!\n";
}

# 输出表头：列名用下划线连接
print $OUT join("\t", map { s/ /_/gr } @target), "\n";

# 提取目标列
while (my $line = <$IN>) {
    chomp $line;
    $line =~ s/\r$//;

    next if $line =~ /^\s*$/;

    my @fields = split /\t/, $line, -1;

    my @values = map {
        defined $fields[$col{$_}] ? $fields[$col{$_}] : ""
    } @target;

    print $OUT join("\t", @values), "\n";
}

close $IN;
close $OUT;

print "Done! Output: $output\n";

