# Workflow

ViruScanNIPT has two separate execution phases. The main program prepares sample-level shell scripts; the cohort Summary is run independently after all sample-level analyses finish.

## 1. Configuration and preflight checks

Run `ViruScanNIPT.pl` with:

- a current-format configuration file;
- a FASTQ list containing one single-end FASTQ path per line;
- an output directory.

Every configuration key in `conf/ViruScanNIPT.conf.example` is mandatory and must occur exactly once. Unknown keys, duplicate keys, missing values and legacy aliases are rejected. All configured executable, reference, database-prefix and annotation-table paths must be absolute.

Before generating shells, the program checks:

- BWA, SAMtools and BLASTN are executable;
- the human reference FASTA and its five BWA index components (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) exist;
- the viral nucleotide BLAST database is available as a `.nal` alias or has `.nhr`, `.nin` and `.nsq` index files;
- both annotation tables exist and are non-empty;
- the Groupname table contains `ACCESSION` and `Groupname` columns;
- numeric configuration values are in the required ranges.

## 2. Sample-level shell generation

`ViruScanNIPT.pl` calls `bin/makeSampleShells.pl`. The output directory contains one folder per sample, with no extra workflow-stage directories:

```text
<outdir>/
├── <sampleA>/
│   └── <sampleA>.run.sh
└── <sampleB>/
    └── <sampleB>.run.sh
```

The generated shell is executable and is not run automatically. It uses Bash double-quoted arguments. Paths containing a single quote or newline cannot be represented in the required shell format and are rejected.

## 3. Per-sample analysis

Each `<sample>.run.sh` executes the following commands sequentially.

### 3.1 Align reads to the human reference

```bash
bwa mem -t BWA_THREADS HG38_REF sample.fq > sample.sam
```

The workflow currently accepts one FASTQ file per sample and therefore runs single-end BWA-MEM.

### 3.2 Convert, sort and index alignments

SAMtools converts SAM to BAM, sorts the BAM and creates its index:

```text
<sample>.bam
<sample>.sort.bam
<sample>.sort.bam.bai
```

### 3.3 Extract unmapped reads

Reads with the SAM unmapped flag are extracted from the sorted BAM:

```bash
samtools view -b -f 4 <sample>.sort.bam
samtools fasta <sample>.unmapped.bam
```

The resulting files are:

- `<sample>.unmapped.bam`
- `<sample>.unmapped.fa`

### 3.4 Search against the viral nucleotide database

BLASTN is run against `VIRAL_DB` using outfmt 7 and the configured `MAX_TARGET_SEQS`, `MAX_HSPS` and `BLAST_THREADS` values.

The raw BLAST output is:

- `<sample>.blast.out`

### 3.5 Parse and filter BLAST results

`bin/parseBlast.pl` applies the configured filters:

- subject accession begins with `NC_`, `AC_` or `NG_`;
- percent identity is greater than or equal to `IDENTITY_CUTOFF`;
- alignment length is greater than or equal to `ALIGNMENT_LENGTH_CUTOFF`;
- e-value is less than or equal to `EVALUE_CUTOFF`.

The parser writes:

- `<sample>_blast_raw_result.xls`: tab-delimited BLAST alignment fields for retained RefSeq-style accession prefixes, before the identity, length and e-value filters;
- `<sample>_blast_out_parsed.txt`: one representative qualifying subject hit per query, with the accession version suffix removed and the organism name looked up in `LOCUS_ORGANISM`.

When multiple hits pass the filters for a query, the current parser retains the first qualifying hit in the BLAST output order. BLAST output ordering is therefore part of the current best-hit behavior.

The parsed table has these columns:

```text
sample  query  subject  identity  alignment_length  evalue  viral_name
```

`subject` is versionless in the parsed output. `viral_name` is `NA` when the accession is not present in `locusOrganism.tsv`.

## 4. Intermediate-file retention

The `KEEP_INTERMEDIATE` setting controls only the selected alignment files:

- `KEEP_INTERMEDIATE=1`: retain SAM, BAM, sorted BAM and sorted BAM index files.
- `KEEP_INTERMEDIATE=0`: remove the SAM, BAM, sorted BAM and sorted BAM index after the sample shell reaches the cleanup step successfully.

The following files are retained in either mode:

- unmapped BAM and FASTA;
- BLAST output;
- parsed and raw parsed-result tables;
- BWA, SAMtools, BLAST and parser logs;
- the generated sample shell.

If a command fails, `set -euo pipefail` stops the shell before the final cleanup step, so intermediate files may remain for troubleshooting.

## 5. Cohort-level Summary

After all sample-level shells finish, run `bin/summarizeGroupnameResults.pl`. The main program prints a Summary command using the configured `GROUPNAME_TABLE` and `POSITIVE_HIT_CUTOFF`.

The Summary script:

1. Recursively finds files ending in `_blast_out_parsed.txt` under the input directory, or accepts one parsed file directly.
2. Reads the `sample`, `query` and `subject` columns.
3. Removes accession version suffixes before matching `subject` to `ACCESSION` in `viralGroupname.tsv`.
4. Excludes parsed rows whose `viral_name` contains `phage`, and does not assign rows with missing or `NA` Groupname values.
5. Counts each distinct query once per sample per Groupname. Multiple accessions assigned to the same Groupname do not inflate the count for the same query.
6. Calls a sample positive for a Groupname when the distinct query count is greater than or equal to `POSITIVE_HIT_CUTOFF`.
7. Writes sample-level, Groupname-level, co-detection and unmapped-accession tables.

If more than one parsed file represents the same sample ID, their query records are combined. A header-only parsed file still represents a sample; the sample ID is inferred from the filename. If a child sample directory contains a generated `.run.sh` but no parsed file, the Summary script stops with an error rather than silently omitting that sample from the cohort denominator.

See [outputFormat.md](outputFormat.md) for the exact output columns and statistical definitions.
