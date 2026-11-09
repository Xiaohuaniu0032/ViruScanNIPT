# ViruScanNIPT

ViruScanNIPT is a research workflow for screening viral sequence signals in
single-end, low-depth noninvasive prenatal testing (NIPT) FASTQ data. The
workflow aligns reads to a human reference, extracts unmapped reads, searches
them against a viral nucleotide BLAST database, and summarizes the resulting
viral hits using curated `Groupname` annotations.

The main program **generates sample-level shell scripts but does not execute
them**. Cohort-level Summary is a separate step that runs after the sample
analyses are complete.

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
BWA-MEM alignment to the human reference
        |
        v
SAMtools BAM conversion, sorting and indexing
        |
        v
Extraction of unmapped reads and FASTA conversion
        |
        v
BLASTN search against the viral nucleotide database
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

For details, see [workflow.md](workflow.md). Output columns and counting rules
are documented in [outputFormat.md](outputFormat.md).

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

- `ViruScanNIPT.pl` validates the configuration and required resources, then
  invokes the sample-shell generator.
- `bin/` contains the scripts used by the active sample-level and cohort
  Summary workflows.
- `conf/` contains the current configuration template. Legacy configuration
  aliases are not supported.
- `database/annotation/` contains the runtime accession-to-organism and
  accession-to-Groupname lookup tables.
- `database/metadata/` contains RefSeq-derived metadata snapshots and the
  Virus-Host DB association table.
- `database/scripts/` contains utilities for preparing or updating database
  metadata. See [database/scripts/README.md](database/scripts/README.md).
- `examples/` contains an example FASTQ list.
- `tests/simulated_FASTQ/` contains ten synthetic positive-control FASTQ files.
- `tests/ViruScanNIPT_result/` contains parsed BLAST examples for the ten
  simulated samples. The existing `cmd.sh` and `all.sh` helper scripts are
  retained. The paths in `all.sh` are specific to the system where it was
  created and should be reviewed before use elsewhere.
- `tests/cmd.sh` is also retained from the original project.

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
- NCBI BLAST+ with `blastn`
- A human reference FASTA indexed for BWA
- A viral nucleotide BLAST database
- `database/annotation/locusOrganism.tsv`
- `database/annotation/viralGroupname.tsv`

The active Perl workflow uses standard Perl modules. Python is not required
for sample-level analysis or cohort Summary.

### Database preparation requirements

Python 3 is required only when regenerating `viralGroupname.tsv`. The extractor
uses the Python standard library and does not require third-party packages.
`gzip` is required when reading compressed GBFF input. The accession extraction
utilities use `IO::Uncompress::Gunzip`.

## 4. Database annotation resources

The two runtime annotation tables have different sources:

- `locusOrganism.tsv` is extracted from the `LOCUS` and `ORGANISM` fields of
  NCBI RefSeq Viral GenBank flat files.
- `viralGroupname.tsv` is generated from
  `database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`.
  Only the `ACCESSION` and `Groupname` columns are copied.

Regenerate `viralGroupname.tsv` from the repository root with:

```bash
python3 database/scripts/extractViralGroupname.py
```

The extractor preserves the source row order and accession identifiers. It
does not remove accession version suffixes or recalculate Groupname labels.
The Summary program normalizes accession version suffixes when matching parsed
BLAST subjects to the mapping table.

See [database/annotation/README.md](database/annotation/README.md) for the
source, format and update procedure for both annotation tables.

## 5. Configuration

Create a working configuration file:

```bash
cp conf/ViruScanNIPT.conf.example conf/ViruScanNIPT.conf
```

Edit every value in `conf/ViruScanNIPT.conf`. All keys in the current template
are mandatory and must appear exactly once.

| Key | Purpose |
|---|---|
| `BWA` | Absolute path to the BWA executable |
| `SAMTOOLS` | Absolute path to the SAMtools executable |
| `BLASTN` | Absolute path to the `blastn` executable |
| `HG38_REF` | Absolute path to the human reference FASTA |
| `VIRAL_DB` | Absolute prefix of the viral nucleotide BLAST database |
| `LOCUS_ORGANISM` | Absolute path to `database/annotation/locusOrganism.tsv` |
| `GROUPNAME_TABLE` | Absolute path to `database/annotation/viralGroupname.tsv` |
| `BWA_THREADS` | Number of BWA threads per sample |
| `BLAST_THREADS` | Number of BLASTN threads per sample |
| `MAX_TARGET_SEQS` | Maximum target sequences passed to BLASTN |
| `MAX_HSPS` | Maximum HSPs passed to BLASTN |
| `EVALUE_CUTOFF` | BLAST e-value threshold |
| `IDENTITY_CUTOFF` | Minimum percent identity |
| `ALIGNMENT_LENGTH_CUTOFF` | Minimum alignment length in bases |
| `POSITIVE_HIT_CUTOFF` | Minimum number of distinct query reads required for a Groupname-positive call |
| `KEEP_INTERMEDIATE` | `1` retains selected SAM/BAM files. `0` removes them after a successful sample run. |

Configuration validation is strict:

- Only keys listed in `conf/ViruScanNIPT.conf.example` are accepted.
- Missing, duplicated, unknown or empty entries cause an error.
- Legacy aliases such as `HUMAN_DNA_VIRUS_DB` are not supported.
- Executable, reference, database-prefix and annotation-table paths must be
  absolute.
- Inline comments after values are not supported. Put comments on separate
  lines.
- The configured BWA reference must have `.amb`, `.ann`, `.bwt`, `.pac` and
  `.sa` index files.
- The viral nucleotide BLAST database must have a `.nal` alias file or the
  `.nhr`, `.nin` and `.nsq` index files.
- `GROUPNAME_TABLE` must contain `ACCESSION` and `Groupname` columns. Both
  values must be non-empty, and each accession must be unique.

## 6. Input FASTQ list

The FASTQ list contains one single-end FASTQ path per line. Blank lines and
lines beginning with `#` are ignored.

Supported file suffixes:

- `.fq`
- `.fq.gz`
- `.fastq`
- `.fastq.gz`

Example:

```text
/path/to/sampleA.fq.gz
/path/to/sampleB.fastq.gz
```

Relative FASTQ paths are resolved first from the current working directory and
then from the directory containing the FASTQ list.

The sample name is derived from the FASTQ filename after removing its FASTQ
suffix. Duplicate sample names are rejected. Valid sample names contain only
letters, digits, underscores, periods and hyphens.

## 7. Generate and run sample shells

Run the following command from the repository root:

```bash
perl ViruScanNIPT.pl \
  --config conf/ViruScanNIPT.conf \
  --fq-list examples/fq.list.example \
  --outdir ViruScanNIPT_result
```

The program creates one sample directory directly under `--outdir`:

```text
ViruScanNIPT_result/
├── sampleA/
│   └── sampleA.run.sh
└── sampleB/
    └── sampleB.run.sh
```

The generated shell is executable and is not run automatically. It uses
double-quoted arguments. Paths containing a single quote or newline are
rejected by the shell generator.

Run a sample manually:

```bash
bash ViruScanNIPT_result/sampleA/sampleA.run.sh
```

For multiple samples, submit each generated `*.run.sh` through the local HPC
scheduler. Each shell runs the sample-level commands sequentially and stops
when a command fails.

After shell generation, `ViruScanNIPT.pl` prints the cohort Summary command
using `GROUPNAME_TABLE` and `POSITIVE_HIT_CUTOFF` from the configuration. Run
that command after all sample-level parsed result files are ready.

## 8. Run cohort Summary

The Summary program is independent of BWA, SAMtools and BLASTN. It reads the
parsed result files and maps each `subject` accession to a curated Groupname.

Run the Summary program with:

```bash
perl bin/summarizeGroupnameResults.pl \
  --input ViruScanNIPT_result \
  --group-table database/annotation/viralGroupname.tsv \
  --outdir ViruScanNIPT_result/Summary \
  --positive-hit-cutoff 2
```

Replace `2` with the `POSITIVE_HIT_CUTOFF` value in the configuration if it
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
`_blast_out_parsed.txt`. Files with the same sample ID are combined. Each query
is counted once per sample per Groupname. A sample directory containing a
generated `.run.sh` but no parsed result file is treated as incomplete and
causes Summary to stop rather than silently reducing the cohort denominator.

See [outputFormat.md](outputFormat.md) for all output columns and counting
rules.

## 9. Simulated FASTQ data

`tests/simulated_FASTQ/` contains ten single-end FASTQ files. Each file has
1,000 reads of 35 bp, for 10,000 reads in total. The reads are sampled from
selected DNA-virus reference regions for HSV-1, human cytomegalovirus, HPV16
and varicella-zoster virus.

This is an enriched positive-control simulation. It contains no sequencing
errors and no human-background reads, so it does not represent the full
complexity or abundance distribution of clinical NIPT cfDNA.

The corresponding parsed BLAST examples are stored under
`tests/ViruScanNIPT_result/SimSample*/`. They can be used to run the cohort
Summary without rerunning BWA or BLASTN.

See [tests/simulated_FASTQ/README.md](tests/simulated_FASTQ/README.md) for
the simulation design, reference traceability and limitations.

To regenerate the simulated FASTQ files into a separate directory, run:

```bash
python3 tests/simulated_FASTQ/generate_simulated_fastq.py \
  --outdir tests/simulated_FASTQ_regenerated
```

The generator uses a fixed random seed by default. The example command writes
to a separate directory and does not overwrite the committed simulation files.

## 10. Interpretation

ViruScanNIPT is a research workflow for screening viral sequence signals in
NIPT sequencing data. A BLAST hit or Groupname-level positive call is not, by
itself, evidence of active infection, maternal infection, fetal infection or
a clinical diagnosis.

Interpretation requires review of read-level evidence, background and
contamination controls, reference and database composition, sample context,
and orthogonal clinical validation.
