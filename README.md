# ViruScanNIPT

ViruScanNIPT is a research workflow for screening viral sequence signals in
single-end, low-depth noninvasive prenatal testing (NIPT) FASTQ data. It first
aligns reads to a human reference, extracts unmapped reads, searches those reads
against a viral nucleotide BLAST database, and summarizes parsed viral hits
using curated `Groupname` annotations.

The main program **generates sample-level shell scripts; it does not execute
them**. Cohort Summary is run separately after sample-level analyses finish.

## 1. Workflow overview

```text
Single-end FASTQ files
        |
        v
ViruScanNIPT.pl
  - validates the current configuration
  - creates one directory per sample
  - writes one executable <sample>.run.sh per sample
        |
        v
Run the generated sample shells
        |
        v
BWA-MEM against the human reference
        |
        v
SAMtools conversion, sorting and indexing
        |
        v
Extract unmapped reads and convert to FASTA
        |
        v
BLASTN against the viral nucleotide database
        |
        v
bin/parseBlast.pl
        |
        v
<sample>_blast_out_parsed.txt
        |
        v
bin/summarizeGroupnameResults.pl
        |
        v
Sample-level and Groupname-level Summary tables
```

See [workflow.md](workflow.md) for the detailed execution flow and
[outputFormat.md](outputFormat.md) for the output columns and counting rules.

## 2. Repository structure

```text
ViruScanNIPT/
├── .gitignore
├── README.md
├── workflow.md
├── outputFormat.md
├── ViruScanNIPT.pl
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
│   │   └── viralGroupname.tsv
│   ├── metadata/
│   │   ├── README.md
│   │   ├── virushostdb.virus_host.tsv
│   │   ├── viral.genomic.metadata.tsv
│   │   ├── viral.genomic.metadata.with_genus.tsv
│   │   ├── viral.genomic.metadata.with_genus.genomicDNA.tsv
│   │   ├── viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv
│   │   ├── viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv
│   │   └── viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv
│   └── scripts/
│       ├── README.md
│       ├── addGenusToMetadata.pl
│       ├── annotateVirusHostDbHost.pl
│       ├── extractAccessionFromGbff.pl
│       ├── extractAccessionTaxid.pl
│       ├── extractGbffMetadata.pl
│       ├── extractLocusOrganism.pl
│       ├── extractViralGroupname.py
│       ├── extractVirusHostInfo.pl
│       ├── filterViralGenomicDnaHost.pl
│       ├── filterVirusHostDbHumanHost.pl
│       └── parseBlastGenus.pl
│
├── examples/
│   └── fq.list.example
│
└── tests/
    ├── cmd.sh
    ├── simulated_FASTQ/
    │   ├── README.md
    │   ├── SimSample01.fq ... SimSample10.fq
    │   ├── fq.list
    │   ├── generate_simulated_fastq.py
    │   ├── reference_segments.fasta
    │   ├── reference_sources.tsv
    │   ├── sample_composition.tsv
    │   ├── simulation_manifest.tsv
    │   └── SHA256SUMS.txt
    └── ViruScanNIPT_result/
        ├── cmd.sh
        ├── all.sh
        └── SimSample01 ... SimSample10/
            └── *_blast_out_parsed.txt
```

### Main directories and files

- `ViruScanNIPT.pl`: validates the current configuration, checks required
  resources and invokes the sample-shell generator.
- `bin/`: scripts used by the active sample-level and cohort Summary workflows.
- `conf/`: the current configuration template. Old configuration aliases are
  not supported.
- `database/annotation/`: runtime accession-to-organism and
  accession-to-Groupname tables.
- `database/metadata/`: RefSeq-derived metadata snapshots and the Virus-Host DB
  association table. These files support annotation-resource preparation.
- `database/scripts/`: utilities for building or updating database metadata.
  See [database/scripts/README.md](database/scripts/README.md).
- `examples/`: example FASTQ-list file.
- `tests/simulated_FASTQ/`: ten synthetic positive-control FASTQ files.
- `tests/ViruScanNIPT_result/`: parsed BLAST result examples for the ten
  simulated samples. The existing `cmd.sh` and `all.sh` helper scripts are
  retained. `all.sh` contains machine-specific paths and should be reviewed
  before use on another system.
- `tests/cmd.sh`: existing helper shell retained in the repository.

Large human-reference FASTA files, BWA indexes and viral BLAST database files
should be stored outside the repository. Do not commit production sequencing
data or local runtime configuration files.

## 3. Requirements

### Runtime requirements

- Linux
- Bash
- Perl
- BWA with the `mem` command
- SAMtools
- NCBI BLAST+ (`blastn`)
- Human reference FASTA indexed for BWA
- Viral nucleotide BLAST database
- `database/annotation/locusOrganism.tsv`
- `database/annotation/viralGroupname.tsv`

The active Perl workflow uses standard Perl modules. Python is not required
for sample-level analysis or cohort Summary.

### Database-resource preparation

- Python 3 is required only when regenerating `viralGroupname.tsv`; the
  extractor uses the Python standard library.
- `gzip` is required when reading compressed GBFF input.
- The accession extraction utilities use `IO::Uncompress::Gunzip`.

## 4. Configuration

Create a working configuration file:

```bash
cp conf/ViruScanNIPT.conf.example conf/ViruScanNIPT.conf
```

Edit every value in `conf/ViruScanNIPT.conf`. All keys in the current template
are mandatory and must appear exactly once.

| Key | Required value |
|---|---|
| `BWA` | Absolute path to the BWA executable |
| `SAMTOOLS` | Absolute path to the SAMtools executable |
| `BLASTN` | Absolute path to the `blastn` executable |
| `HG38_REF` | Absolute path to the human reference FASTA |
| `VIRAL_DB` | Absolute prefix of the viral nucleotide BLAST database |
| `LOCUS_ORGANISM` | Absolute path to `database/annotation/locusOrganism.tsv` |
| `GROUPNAME_TABLE` | Absolute path to `database/annotation/viralGroupname.tsv` |
| `BWA_THREADS` | Positive integer; BWA threads per sample |
| `BLAST_THREADS` | Positive integer; BLASTN threads per sample |
| `MAX_TARGET_SEQS` | Positive integer passed to BLASTN |
| `MAX_HSPS` | Positive integer passed to BLASTN |
| `EVALUE_CUTOFF` | Positive numeric BLAST e-value threshold |
| `IDENTITY_CUTOFF` | Percent identity threshold, greater than 0 and at most 100 |
| `ALIGNMENT_LENGTH_CUTOFF` | Positive integer minimum alignment length |
| `POSITIVE_HIT_CUTOFF` | Positive integer number of distinct query reads required for a Groupname-positive call |
| `KEEP_INTERMEDIATE` | `1` keeps selected SAM/BAM files; `0` removes selected alignment intermediates after a successful run |

Configuration validation is strict:

- Only keys in `conf/ViruScanNIPT.conf.example` are accepted.
- Missing, duplicated, unknown or empty entries cause an error.
- Legacy aliases such as `HUMAN_DNA_VIRUS_DB` are not supported.
- Executable, reference, database-prefix and annotation-table paths must be
  absolute.
- Inline comments after values are not supported; put comments on separate
  lines.
- The configured BWA reference must have `.amb`, `.ann`, `.bwt`, `.pac` and
  `.sa` index files.
- The viral nucleotide BLAST database must have a `.nal` alias file or the
  `.nhr`, `.nin` and `.nsq` index files.
- `GROUPNAME_TABLE` must contain the `ACCESSION` and `Groupname` columns.
  Both values must be non-empty and each accession must be unique.

## 5. Input FASTQ list

The FASTQ list contains one single-end FASTQ path per line. Blank lines and
lines beginning with `#` are ignored.

Supported suffixes:

- `.fq`
- `.fq.gz`
- `.fastq`
- `.fastq.gz`

Example:

```text
/path/to/sampleA.fq.gz
/path/to/sampleB.fastq.gz
```

Relative FASTQ paths are resolved relative to the current working directory
or the directory containing the FASTQ list. Sample names are derived from the
FASTQ filename after removing its FASTQ suffix. Duplicate sample names are
rejected. Valid sample names contain only letters, digits, underscore,
period and hyphen.

## 6. Generate and run sample shells

Run from the project root:

```bash
perl ViruScanNIPT.pl \
  --config conf/ViruScanNIPT.conf \
  --fq-list examples/fq.list.example \
  --outdir ViruScanNIPT_result
```

The command creates one sample directory directly under `--outdir`:

```text
ViruScanNIPT_result/
├── sampleA/
│   └── sampleA.run.sh
└── sampleB/
    └── sampleB.run.sh
```

The generated shell is executable and is not run automatically. It uses
double-quoted arguments. Paths containing a single quote or newline are not
supported by the shell generator.

Run a sample manually:

```bash
bash ViruScanNIPT_result/sampleA/sampleA.run.sh
```

For multiple samples, submit each generated `*.run.sh` through the local HPC
scheduler. The sample shell runs commands sequentially and stops after a
failed command (`set -euo pipefail`).

After shell generation, `ViruScanNIPT.pl` prints the cohort Summary command
using `GROUPNAME_TABLE` and `POSITIVE_HIT_CUTOFF` from the configuration. Run
that command after all sample-level parsed result files are ready.

## 7. Run cohort Summary

The Summary program is independent of BWA, SAMtools and BLASTN. It reads the
parsed result files and maps each `subject` accession to a Groupname.

The equivalent command is:

```bash
perl bin/summarizeGroupnameResults.pl \
  --input ViruScanNIPT_result \
  --group-table database/annotation/viralGroupname.tsv \
  --outdir ViruScanNIPT_result/Summary \
  --positive-hit-cutoff 2
```

Replace `2` with the `POSITIVE_HIT_CUTOFF` value in your configuration if it
differs.

A single parsed file can also be summarized:

```bash
perl bin/summarizeGroupnameResults.pl \
  --input ViruScanNIPT_result/sampleA/sampleA_blast_out_parsed.txt \
  --group-table database/annotation/viralGroupname.tsv \
  --outdir ViruScanNIPT_result/sampleA/Summary \
  --positive-hit-cutoff 2
```

The Summary program recursively discovers files ending in
`_blast_out_parsed.txt`. It combines files with the same sample ID and counts
each query identifier once per sample per Groupname. A sample directory with
a generated `.run.sh` but no parsed result file is treated as an incomplete
sample and causes Summary to stop rather than silently reducing the cohort
denominator.

See [outputFormat.md](outputFormat.md) for the output columns and counting
rules.

## 8. Database resources

The sample-level parser uses `locusOrganism.tsv` to annotate BLAST subject
accessions with organism names. Cohort Summary uses `viralGroupname.tsv` to
map accessions to curated Groupname labels.

The two files have different sources:

- `locusOrganism.tsv` is extracted from the `LOCUS` and `ORGANISM` fields of
  NCBI RefSeq Viral GBFF records.
- `viralGroupname.tsv` is generated from
  `database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`.
  The extractor copies only the `ACCESSION` and `Groupname` columns, retaining
  the accession identifiers as supplied in the source file.

See [database/annotation/README.md](database/annotation/README.md) for
annotation-table sources, formats and rebuild commands.

`virushostdb.virus_host.tsv` is a general virus-host association table
extracted from Virus-Host DB. It retains recorded host associations; filtering
for `Homo sapiens` is a separate downstream operation. See
[database/metadata/README.md](database/metadata/README.md).

## 9. Simulated FASTQ data

`tests/simulated_FASTQ/` contains ten single-end FASTQ files, each with 1,000
reads of 35 bp. The reads are sampled from selected DNA-virus reference
regions for HSV-1, human cytomegalovirus, HPV16 and varicella-zoster virus.
This is an enriched positive-control simulation with no sequencing errors
and no human-background reads; it is not a model of the full complexity or
abundance distribution of clinical NIPT cfDNA.

The corresponding parsed BLAST result examples are stored under
`tests/ViruScanNIPT_result/SimSample*/`. They can be used as example inputs
for the cohort Summary program without rerunning BWA or BLASTN.

See [tests/simulated_FASTQ/README.md](tests/simulated_FASTQ/README.md) for
simulation design, reference traceability and limitations.

Regenerate the simulated FASTQ files into a separate directory:

```bash
python3 tests/simulated_FASTQ/generate_simulated_fastq.py \
  --outdir tests/simulated_FASTQ_regenerated
```

The generator uses a fixed random seed by default and does not overwrite the
committed simulation files when the separate output directory is used.

## 10. Interpretation

ViruScanNIPT is a research workflow for screening viral sequence signals in
NIPT sequencing data. A BLAST hit or Groupname-level positive call is not, by
itself, evidence of active infection, maternal infection, fetal infection or a
clinical diagnosis. Interpretation requires review of read-level evidence,
background and contamination controls, reference/database composition,
sample context and orthogonal clinical validation.
