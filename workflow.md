# Workflow

ViruScanNIPT has two execution phases. The main program generates sample-level
shell scripts. Cohort-level Summary is run separately after the sample
analyses have finished.

## 1. Configuration and preflight checks

Run `ViruScanNIPT.pl` with three required inputs:

- A configuration file in the current format.
- A FASTQ list containing one single-end FASTQ path per line.
- An output directory.

Every key in `conf/ViruScanNIPT.conf.example` is mandatory and must occur
exactly once. Unknown keys, duplicate keys, missing values and legacy aliases
are rejected. All configured executable, reference, database-prefix and
annotation-table paths must be absolute.

Before generating sample shells, the main program checks that:

- BWA, SAMtools and BLASTN are executable.
- The human reference FASTA and its five BWA index components exist:
  `.amb`, `.ann`, `.bwt`, `.pac` and `.sa`.
- The viral nucleotide BLAST database is available as a `.nal` alias or has
  the `.nhr`, `.nin` and `.nsq` index files.
- Both annotation tables exist and are non-empty.
- The Groupname table contains `ACCESSION` and `Groupname` columns.
- Required numeric settings are within their permitted ranges.

## 2. Sample-level shell generation

`ViruScanNIPT.pl` calls `bin/makeSampleShells.pl`. The program creates one
sample directory directly under the specified output directory. It does not
create additional workflow-stage directories.

```text
<outdir>/
├── <sampleA>/
│   └── <sampleA>.run.sh
└── <sampleB>/
    └── <sampleB>.run.sh
```

Each generated shell is executable but is not run automatically. Arguments
are written using Bash double quotes. Paths containing a single quote or a
newline are rejected because they cannot be represented in the required
generated-shell format.

The sample name is derived from the FASTQ filename after removing `.fq`,
`.fq.gz`, `.fastq` or `.fastq.gz`. Duplicate sample names are not allowed.

## 3. Per-sample analysis

Each `<sample>.run.sh` executes the following steps in order.

### 3.1 Align reads to the human reference

BWA-MEM aligns the single-end FASTQ reads to the configured human reference:

```bash
bwa mem -t BWA_THREADS HG38_REF sample.fq > sample.sam
```

The current workflow accepts one FASTQ file per sample and therefore uses
single-end alignment.

### 3.2 Convert, sort and index alignments

SAMtools converts the SAM file to BAM, sorts the BAM file and creates its
index. The main alignment files are:

```text
<sample>.bam
<sample>.sort.bam
<sample>.sort.bam.bai
```

### 3.3 Extract unmapped reads

Reads carrying the SAM unmapped flag are extracted from the sorted BAM file.
SAMtools then converts the unmapped BAM records to FASTA.

```bash
samtools view -b -f 4 -o <sample>.unmapped.bam <sample>.sort.bam
samtools fasta <sample>.unmapped.bam > <sample>.unmapped.fa
```

The resulting files are:

- `<sample>.unmapped.bam`
- `<sample>.unmapped.fa`

### 3.4 Search against the viral nucleotide database

BLASTN searches the unmapped FASTA reads against `VIRAL_DB`. The workflow
uses outfmt 7 and the configured `MAX_TARGET_SEQS`, `MAX_HSPS` and
`BLAST_THREADS` values.

The raw BLAST output is:

- `<sample>.blast.out`

### 3.5 Parse and filter BLAST results

`bin/parseBlast.pl` applies the configured filters:

- The subject accession must begin with `NC_`, `AC_` or `NG_`.
- Percent identity must be greater than or equal to `IDENTITY_CUTOFF`.
- Alignment length must be greater than or equal to
  `ALIGNMENT_LENGTH_CUTOFF`.
- E-value must be less than or equal to `EVALUE_CUTOFF`.

The parser writes two files:

- `<sample>_blast_raw_result.xls` contains tab-delimited BLAST alignment
  fields for the retained RefSeq-style accession prefixes. Identity,
  alignment-length and e-value filters have not yet been applied to this
  output.
- `<sample>_blast_out_parsed.txt` contains one representative qualifying
  subject hit per query. The subject accession version suffix is removed,
  and the organism name is looked up in `LOCUS_ORGANISM`.

When multiple hits pass the filters for one query, the parser retains the
first qualifying hit in BLAST output order. The selected hit therefore
depends on the order of records in the BLAST output.

The parsed table contains these columns:

```text
sample  query  subject  identity  alignment_length  evalue  viral_name
```

The `subject` value is versionless. `viral_name` is `NA` when the accession
cannot be found in `locusOrganism.tsv`.

## 4. Intermediate-file retention

`KEEP_INTERMEDIATE` controls whether selected alignment files are retained:

- `KEEP_INTERMEDIATE=1` retains the SAM, BAM, sorted BAM and sorted BAM index.
- `KEEP_INTERMEDIATE=0` removes those files after the sample shell reaches
  its cleanup step successfully.

The following files are retained in both modes:

- Unmapped BAM and FASTA files.
- Raw BLAST output.
- Parsed and raw parsed-result tables.
- BWA, SAMtools, BLAST and parser logs.
- The generated sample shell.

Each generated shell uses `set -euo pipefail`. If a command fails, execution
stops before the final cleanup step. Intermediate files may therefore remain
for troubleshooting.

## 5. Cohort-level Summary

After all sample-level shells have finished, run
`bin/summarizeGroupnameResults.pl`. The main program prints a Summary command
using `GROUPNAME_TABLE` and `POSITIVE_HIT_CUTOFF` from the configuration.

The Summary program performs the following operations:

1. Recursively finds files ending in `_blast_out_parsed.txt` under the input
   directory. It can also accept one parsed file directly.
2. Reads the `sample`, `query` and `subject` columns.
3. Removes accession version suffixes before matching `subject` to
   `ACCESSION` in `viralGroupname.tsv`.
4. Excludes records whose `viral_name` contains `phage`, ignoring letter case.
   Records with missing or `NA` Groupname assignments are not counted.
5. Counts each distinct query once per sample per Groupname. Multiple
   accessions assigned to the same Groupname do not increase the count for
   the same sample-query pair.
6. Calls a sample positive for a Groupname when its distinct query count is
   greater than or equal to `POSITIVE_HIT_CUTOFF`.
7. Writes sample-level, Groupname-level, co-detection and unmapped-accession
   tables.

If multiple parsed files represent the same sample ID, their records are
combined before counting. A header-only parsed file still represents a
sample. In that case, the sample ID is inferred from the filename.

If a child sample directory contains a generated `.run.sh` but has no parsed
result file, Summary stops with an error. This prevents an incomplete sample
from being silently omitted from the cohort denominator.

See [outputFormat.md](outputFormat.md) for the output columns, counting units
and statistical definitions.

## 6. Interpretation

The Summary positive status is an operational threshold based on distinct
query reads assigned to a Groupname. It is not a clinical diagnosis and does
not establish active infection, maternal infection or fetal infection.
Interpretation requires read-level review, suitable controls, database
assessment and clinical validation.
