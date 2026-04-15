# ViruScanNIPT

ViruScanNIPT is a lightweight pipeline framework for viral sequence screening from low-depth noninvasive prenatal testing (NIPT) sequencing data. It is designed to organize sample-level processing and cohort-level summarization in a modular and reproducible way.

## Overview

### Main features
- Human read removal by alignment to the reference genome
- Viral sequence search against a BLAST database
- BLAST result parsing and filtering
- Viral accession annotation to organism names and predefined human DNA virus groups
- Sample-level group status calling
- Cohort-level ranking, co-detection analysis, and distribution statistics
- Downstream visualization support

### Intended use
This pipeline is intended for research use, especially for retrospective virome analysis of low-depth NIPT sequencing datasets.

---

## Workflow

The general workflow is:

1. Align raw FASTQ reads to the human reference genome
2. Extract non-human or unmapped reads
3. Search reads against a viral BLAST database
4. Parse and filter BLAST hits
5. Annotate viral hits to predefined virus groups
6. Generate sample-level and cohort-level summary results

`ViruScanNIPT.pl` does **not** directly run the whole analysis.  
Instead, it generates:
- per-sample shell scripts
- summary-level shell scripts
- output directory structure
- runtime copies of key input and configuration files

Users can then execute the generated shell scripts step by step in their own computing environment.

---

## Repository structure

```text
ViruScanNIPT/
├── ViruScanNIPT.pl
├── README.md
├── .gitignore
├── bin/
├── conf/
├── db/
└── test/
```

### Main components
- `ViruScanNIPT.pl`  
  Main framework generator of the pipeline

- `bin/`  
  Auxiliary scripts for downstream analysis, including:
  - `VSN_01_make_sample_shells.pl`
  - `VSN_05_parse_blast.pl`
  - `VSN_06_merge_sample_results.pl`
  - `VSN_07_annotate_human_dna_virus.pl`
  - `VSN_08_call_group_status.pl`
  - `VSN_09_group_positive_ranking.pl`
  - `VSN_10_codetection_analysis.pl`
  - `VSN_11_group_distribution_stats.pl`
  - `ViroStat_Visualizer.py`
  - `extract_locus_organism.pl`

- `conf/`  
  Configuration files, including `ViruScanNIPT.conf.example`

- `db/`  
  Annotation resource files, including:
  - `224_human_host_DNA_viruses.txt`
  - `locus_organism.tsv`

- `test/`  
  Test data and example inputs, including:
  - `fq.list`
  - `cmd.sh`

---

## Requirements

### Operating system
- Linux

### Required software
- Perl
- BWA
- SAMtools
- BLASTN
- Python 3

### Required resources
- Human reference genome
- Viral BLAST database
- Viral accession-to-organism annotation table
- Human DNA virus group annotation table

---

## Configuration

The pipeline uses a plain-text configuration file to store:
- software paths
- reference and database paths
- fixed resource parameters
- BLAST parameters
- filtering thresholds
- group-level calling threshold
- intermediate file retention policy

Example configuration file:

```text
conf/ViruScanNIPT.conf.example
```

For local use:

```bash
cp conf/ViruScanNIPT.conf.example conf/ViruScanNIPT.conf
```

Run-specific parameters such as the FASTQ list and output directory should be provided on the command line.

---

## Database and annotation resources

### VIRAL_DB
`VIRAL_DB` should be built from the NCBI RefSeq viral genome FASTA file.

Example:

```bash
wget https://ftp.ncbi.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gzip -dc viral.1.1.genomic.fna.gz > viral.genomic.fna

makeblastdb \
  -in viral.genomic.fna \
  -dbtype nucl \
  -parse_seqids \
  -title "NCBI viral genomic (nucl)" \
  -out viral
```

In the configuration file, `VIRAL_DB` should point to the BLAST database prefix, for example:

```text
/path/to/viral
```

### LOCUS_ORGANISM
`LOCUS_ORGANISM` is a tab-delimited mapping table from RefSeq viral accession (`LOCUS`) to organism name.

It can be generated from a RefSeq viral GenBank flat file such as `viral.1.genomic.gbff.gz` using the helper script:

```text
bin/extract_locus_organism.pl
```

Example:

```bash
perl bin/extract_locus_organism.pl viral.1.genomic.gbff.gz > locus_organism.tsv
```

### HUMAN_DNA_VIRUS_DB
`HUMAN_DNA_VIRUS_DB` is a manually curated annotation table used to map viral accessions to predefined human-host DNA virus groups.

In this project, `224_human_host_DNA_viruses.txt` was compiled based on a published human DNA virus resource:

- PMID: 35403226

Users may further revise or expand this table according to their own study design and grouping strategy.

---

## Input

The pipeline requires:
- a configuration file
- a FASTQ list file
- raw FASTQ files

### FASTQ list format
One FASTQ path per line, for example:

```text
/path/to/sample1.fq.gz
/path/to/sample2.fq.gz
/path/to/sample3.fq.gz
```

---

## Usage

### Basic command

```bash
perl ViruScanNIPT.pl \
  --config conf/ViruScanNIPT.conf \
  --fq-list test/fq.list \
  --outdir test/ViruScanNIPT_result
```

### Arguments
- `--config` : configuration file
- `--fq-list` : input FASTQ list
- `--outdir` : output directory
- `--help` : print help information

### Example
A minimal example command is already provided in:

```text
test/cmd.sh
```

### Before running
Please confirm that:
- all software paths are correct
- all reference and database files are available
- the FASTQ list is correctly formatted
- the output directory is writable

---

## Output

Running `ViruScanNIPT.pl` generates a working analysis framework rather than final biological results directly.

Typical outputs include:
- sample-level shell scripts
- summary-level shell scripts
- organized output directories
- copied configuration and input files
- log files for pipeline generation

After executing the generated shell scripts, downstream outputs may include:
- parsed BLAST result tables
- merged sample-level viral detection results
- annotated virus tables
- sample-level group status results
- group-level ranking tables
- co-detection analysis results
- group distribution statistics
- visualization-ready summary files

---

## Notes

- This repository is a pipeline **framework generator**, not a single-command end-to-end executor.
- Viral detection results should be interpreted with caution, especially for low-depth sequencing data.
- The reference genome, BLAST database, and annotation resources should be version-consistent.
- Keep `conf/ViruScanNIPT.conf.example` in the repository, and exclude production configuration files from version control.
- Test data and large result directories should also be excluded from version control when not required for public release.

---

## Contact

For questions, bug reports, or suggestions, please open an issue in this repository.
