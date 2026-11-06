# ViruScanNIPT

ViruScanNIPT is a lightweight framework for viral sequence screening from low-depth
noninvasive prenatal testing (NIPT) sequencing data. It generates one executable
sample-level shell per FASTQ and provides a standalone cohort Summary script.

## 1. Workflow

The project has two independent phases:

1. **Sample-level analysis**: `ViruScanNIPT.pl` creates one directory per sample
   directly under `--outdir` and writes one `<sample>.run.sh` in that directory.
   It does not execute the generated shells.
2. **Cohort-level Summary**: after all sample shells finish, run
   `bin/summarizeGroupnameResults.pl`. It reads the parsed BLAST result files and
   maps their `subject` accessions to the curated Groupname table.

The sample workflow itself is:

```text
FASTQ
  |
  v
BWA mem against hg38
  |
  v
SAM/BAM conversion, sorting and indexing
  |
  v
Extract unmapped reads
  |
  v
BAM to FASTA
  |
  v
BLASTN against the viral reference database
  |
  v
parseBlast.pl
  |
  v
<sample>_blast_out_parsed.txt
```

The Summary workflow is independent of BWA, SAMtools and BLASTN:

```text
*_blast_out_parsed.txt + viralGroupname.tsv
                 |
                 v
       summarizeGroupnameResults.pl
                 |
                 v
       sample-level and Groupname-level tables
```

## 2. Repository structure

```text
ViruScanNIPT/
├── ViruScanNIPT.pl
├── README.md
├── .gitignore
│
├── conf/
│   └── ViruScanNIPT.conf.example
│
├── bin/
│   ├── makeSampleShells.pl
│   ├── parseBlast.pl
│   └── summarizeGroupnameResults.pl
│
├── database/
│   ├── annotation/
│   │   ├── README.md
│   │   ├── locusOrganism.tsv
│   │   ├── viralGroupname.tsv
│   │   └── source/
│   │       └── viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx
│   ├── metadata/
│   │   ├── README.md
│   │   ├── virushostdb.virus_host.tsv
│   │   └── RefSeq metadata tables
│   ├── scripts/
│   │   ├── extractViralGroupname.py
│   │   └── database annotation/building utilities
│   └── build/
│       └── accession/taxonomy intermediate files
│
├── tools/
│   └── viroStatVisualizer.py
│
├── docs/
│   ├── workflow.md
│   ├── outputFormat.md
│   └── deprecatedScripts.md
│
├── examples/
│   └── fq.list.example
│
├── tests/
│   ├── fixtures/
│   └── runSmokeTest.sh
│
└── deprecated/
    ├── summary/
    ├── legacyBinCommands.sh
    └── legacyDatabaseCommands.sh
```

### Main directories

- `bin/`: scripts used by the current sample-level workflow and Summary.
- `conf/`: software, reference, database and threshold configuration.
- `database/annotation/`: final accession-to-organism and accession-to-Groupname tables.
  See [`database/annotation/README.md`](database/annotation/README.md) for their
  sources, formats and rebuild instructions.
- `database/metadata/`: RefSeq metadata and the extracted Virus-Host DB reference table.
  See [`database/metadata/README.md`](database/metadata/README.md) for the
  Virus-Host DB source, extraction procedure and downstream use.
- `database/scripts/`: utilities for rebuilding or annotating reference metadata.
- `database/build/`: accession and taxonomy intermediate files.
- `tools/`: optional visualization utilities; not required by the core pipeline.
- `docs/`: workflow and output documentation.
- `tests/`: lightweight syntax and workflow smoke tests.
- `deprecated/`: old summary modules and command snippets retained temporarily for comparison.

Large hg38 and BLAST database files should remain outside the Git repository. Set their absolute paths in the configuration file.

## 3. Requirements

### Core workflow

- Linux
- Perl
- BWA
- SAMtools
- BLASTN
- hg38 reference FASTA
- Viral BLAST nucleotide database
- Accession-to-organism table (sample-level BLAST annotation)
- Accession-to-Groupname table (cohort Summary)

### Optional utilities

- `tools/viroStatVisualizer.py` requires Python 3, pandas, matplotlib and seaborn.
  It is not called by `ViruScanNIPT.pl` or the Summary script.
- Rebuilding `database/annotation/viralGroupname.tsv` from the bundled Excel
  source additionally requires Python 3 and `openpyxl`. This dependency is not
  required to run the core analysis or Summary workflow.

## 4. Configuration

Copy the example configuration:

```bash
cp conf/ViruScanNIPT.conf.example conf/ViruScanNIPT.conf
```

Edit the absolute paths in the configuration file. The key fields are:

| Key | Purpose |
|---|---|
| `BWA` | BWA executable |
| `SAMTOOLS` | SAMtools executable |
| `BLASTN` | BLASTN executable |
| `HG38_REF` | Human reference FASTA |
| `VIRAL_DB` | BLAST database prefix |
| `LOCUS_ORGANISM` | RefSeq accession-to-organism table |
| `GROUPNAME_TABLE` | RefSeq accession-to-Groupname table used by Summary |
| `BWA_THREADS` | BWA threads |
| `BLAST_THREADS` | BLASTN threads |
| `MAX_TARGET_SEQS` | BLASTN maximum target sequences |
| `MAX_HSPS` | BLASTN maximum HSPs |
| `EVALUE_CUTOFF` | BLAST e-value threshold |
| `IDENTITY_CUTOFF` | Minimum percent identity |
| `ALIGNMENT_LENGTH_CUTOFF` | Minimum alignment length |
| `POSITIVE_HIT_CUTOFF` | Distinct query-read threshold for a positive Groupname |
| `KEEP_INTERMEDIATE` | Keep or remove selected SAM/BAM intermediates |

For backward compatibility, `HUMAN_DNA_VIRUS_DB` is accepted as an alias for `GROUPNAME_TABLE`.
New configurations should use `GROUPNAME_TABLE`.

`database/annotation/locusOrganism.tsv` is extracted from the NCBI RefSeq
Viral GenBank flat file by reading each record's `LOCUS` and `ORGANISM` fields.
`database/annotation/viralGroupname.tsv` is generated from the `ACCESSION` and
`Groupname` columns in the `Accession_Grouping` worksheet of the supplied
`viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx` workbook.
The current mapping contains 298 unique accessions across 36 Groupname labels.
Groupname values are copied from the workbook and are not inferred from
Taxonomy, Genus or Virus-Host DB. The workbook and conversion script are retained
in the repository to make the mapping reproducible. See
[`database/annotation/README.md`](database/annotation/README.md) for the exact
source and rebuild commands.

The general virus-host association table `database/metadata/virushostdb.virus_host.tsv`
is extracted from the official Virus-Host DB download. It retains all host
associations; human-host filtering is performed only in downstream metadata
processing. See [`database/metadata/README.md`](database/metadata/README.md)
for the source, extraction procedure and downstream use.

## 5. Input

The FASTQ list contains one FASTQ file path per line. Both `.fq`, `.fq.gz`, `.fastq`
and `.fastq.gz` are supported.

Example:

```text
/path/to/sampleA.fq.gz
/path/to/sampleB.fastq.gz
```

Relative FASTQ paths are resolved first relative to the current working directory and
then relative to the directory containing `fq.list`.

Sample names are derived from the FASTQ filename after removing the FASTQ suffix.
Duplicate sample names are not allowed.

## 6. Generate sample shells

Run:

```bash
perl ViruScanNIPT.pl \
  --config conf/ViruScanNIPT.conf \
  --fq-list examples/fq.list.example \
  --outdir ViruScanNIPT_result
```

The output is intentionally flat:

```text
ViruScanNIPT_result/
├── CL100200740_L01_12/
│   ├── CL100200740_L01_12.run.sh
│   ├── CL100200740_L01_12.blast.out
│   ├── CL100200740_L01_12_blast_out_parsed.txt
│   └── sample log and intermediate files
├── Sample_B/
│   └── ...
└── ...
```

No `01_sample_shells`, `02_summary_shells`, `03_summary_analysis`, `04_logs` or
per-step subdirectories are created. Each sample's generated shell and its analysis
files are stored directly in that sample's directory.

The generated shell:
- uses Bash double-quoted arguments;
- contains no single-quote character;
- is executable;
- runs the sample-level steps sequentially when submitted;
- is not executed automatically by `ViruScanNIPT.pl`.

Paths containing a literal single-quote character are rejected by the shell generator,
because generated sample shells are required to contain no single-quote character.

To run one sample:

```bash
bash ViruScanNIPT_result/CL100200740_L01_12/CL100200740_L01_12.run.sh
```

To submit multiple samples, submit each `*.run.sh` using the local HPC scheduler.

## 7. Run the standalone Summary

After all sample shells have completed and the parsed result file exists in every
sample directory, run:

```bash
perl bin/summarizeGroupnameResults.pl \
  --input ViruScanNIPT_result \
  --group-table database/annotation/viralGroupname.tsv \
  --outdir ViruScanNIPT_result/Summary \
  --positive-hit-cutoff 2
```

A single parsed file can also be summarized:

```bash
perl bin/summarizeGroupnameResults.pl \
  --input ViruScanNIPT_result/CL100200740_L01_12/CL100200740_L01_12_blast_out_parsed.txt \
  --group-table database/annotation/viralGroupname.tsv \
  --outdir ViruScanNIPT_result/CL100200740_L01_12/Summary \
  --positive-hit-cutoff 2
```

The Groupname table must contain an accession column (`ACCESSION` or
`Reference_accession`) and a `Groupname` column. Other annotation columns are allowed.
Accession version suffixes such as `.1` or `.2` are removed before matching.

The Summary script searches recursively for `*_blast_out_parsed.txt`. If a sample
directory contains a generated `.run.sh` but its parsed result is missing, Summary
stops with an error rather than silently reducing the cohort denominator.

## 8. Summary outputs

| File | Description |
|---|---|
| `sampleGroupCounts.tsv` | One row per sample × Groupname, distinct query count and Positive/Negative status. Zero-count combinations are retained. |
| `sampleSummary.tsv` | Per-sample distinct assigned query count, number of groups with at least one read, number of positive groups and their names. |
| `groupSummary.tsv` | Group-level positive sample count/rate, total distinct queries and read-count distribution statistics among positive samples. |
| `groupReadCounts.tsv` | Positive sample-level read counts for each Groupname; also serves as input to the optional visualizer. |
| `groupCoDetection.tsv` | Pairwise Groupname co-detection counts and cohort percentages. |
| `sampleCoDetection.tsv` | Exact positive-Groupname combination for each sample with at least two positive groups. |
| `groupCombinationSummary.tsv` | Frequency of each exact multi-Groupname combination. |
| `unmappedAccessions.tsv` | Parsed accessions not found in the Groupname table, for annotation QC. |
| `summary.log` | Input files, sample count, mapped/unmapped records, threshold and output path. |

### Counting rules

1. A read is represented by its `query` identifier in the parsed BLAST table.
2. Each query is counted once per sample per Groupname, even if duplicate rows or
   multiple reference accessions within the same Groupname are present.
3. A sample is Positive for a Groupname when its distinct query count is greater than
   or equal to `POSITIVE_HIT_CUTOFF` (default 2).
4. The same sample can be positive for multiple Groupnames. Group-specific counts and
   percentages are therefore not mutually exclusive.
5. The cohort denominator is the number of unique samples represented by the parsed
   files supplied to Summary.
6. Rows whose `viral_name` contains `phage` are excluded, preserving the previous merge-stage filter.
7. Rows whose Groupname is `NA` or `N/A` are not treated as valid Groupname assignments.

## 9. Optional visualization

```bash
python tools/viroStatVisualizer.py \
  ViruScanNIPT_result/Summary/groupReadCounts.tsv
```

The plot is written to the same Summary directory.

## 10. Tests

Run the lightweight smoke test:

```bash
bash tests/runSmokeTest.sh
```

The test checks Perl syntax, verifies that sample shells are generated directly under
the output directory without single quotes or nested directories, and tests Summary
with small fixture data. It does not run BWA, SAMtools or BLASTN.

## 11. Interpretation

ViruScanNIPT is a research workflow for screening viral sequence signals in NIPT
sequencing data. A BLAST hit or Groupname-level positive call is not, by itself,
evidence of active infection, fetal infection or a clinical diagnosis. Interpretation
requires consideration of sequencing background, reference/database composition,
read-level evidence, sample context and orthogonal clinical validation.
