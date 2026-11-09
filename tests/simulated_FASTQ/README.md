# Simulated FASTQ positive controls

This directory contains ten single-end FASTQ files. Each file contains
1,000 reads, and each read is 35 bp long. The complete dataset contains
10,000 reads.

## 1. Simulation design

- All reads are derived from DNA-virus reference gene or region fragments
  listed in `reference_segments.fasta`.
- Reads are exact 35-bp substrings sampled from the source fragments.
  Reverse-complement orientations are included.
- No sequencing errors or human-background reads are added. This is an
  intentionally enriched positive-control dataset designed to make viral
  BLAST hits easier to observe.
- Every base has FASTQ quality character `I`, corresponding to Phred+33 Q40.
- Each sample contains 1,000 reads. The virus composition varies by sample
  and is recorded in `simulation_manifest.tsv`.
- Read identifiers include the reference accession, source segment, start
  offset within the included segment and strand.
- The reference accessions are represented in the project's
  `database/annotation/viralGroupname.tsv`: NC_001806, NC_006273,
  NC_001526 and NC_001348.

## 2. Files

| File | Description |
|---|---|
| `SimSample01.fq` to `SimSample10.fq` | Ten simulated single-end FASTQ files. |
| `fq.list` | FASTQ list for `ViruScanNIPT.pl`. Paths are relative to this directory. |
| `simulation_manifest.tsv` | Read counts by sample and reference virus or segment. |
| `sample_composition.tsv` | Expected per-sample read counts by virus. |
| `reference_sources.tsv` | Accession, Groupname and source-page traceability. |
| `reference_segments.fasta` | Exact source fragments used to generate reads. |
| `generate_simulated_fastq.py` | Deterministic FASTQ generator with seed `20261106`. |
| `SHA256SUMS.txt` | SHA-256 checksums for the ten FASTQ files. |

Verify the FASTQ files from this directory:

```bash
sha256sum -c SHA256SUMS.txt
```

## 3. Reference sources and Groupname coverage

The four selected RefSeq accessions are included in
`database/annotation/viralGroupname.tsv`.

| Virus | Accession | Groupname |
|---|---|---|
| HSV-1 | NC_001806.2 | Herpesvirus 1/2 |
| Human cytomegalovirus | NC_006273.2 | Herpesvirus 5 |
| HPV16 | NC_001526.4 | Papillomavirus |
| Varicella-zoster virus | NC_001348.1 | Herpesvirus 3 |

The source pages are recorded in `reference_sources.tsv`. The HSV-1 sequences
are partial excerpts from three gene entries. The CMV gB, HPV16 L1 and VZV gE
sequences correspond to the gene intervals listed in
`reference_segments.fasta`.

The simulated reference accessions include version suffixes. The Groupname
mapping table retains the accession identifiers from the grouped metadata
source, while the Summary program removes version suffixes when matching
parsed BLAST subjects.

## 4. Run with ViruScanNIPT

Run the following command from the repository root. Replace the example
FASTQ-list path with the location of this directory if needed.

```bash
perl ViruScanNIPT.pl \
  --config conf/ViruScanNIPT.conf \
  --fq-list /path/to/simulated_FASTQ/fq.list \
  --outdir ViruScanNIPT_simulation_result
```

The pipeline first aligns reads to hg38 and then runs BLASTN on reads that
remain unmapped. Because this dataset contains only viral reads, most reads
are expected to remain unmapped to hg38 and proceed to viral BLASTN. A small
number may still align to human sequence by chance or because of local
sequence similarity.

## 5. Important limitations

These are synthetic positive-control reads, not real NIPT or cfDNA data.
They are sampled from selected reference gene or region fragments rather
than uniformly across complete viral genomes.

The dataset is not suitable for estimating clinical sensitivity, specificity,
prevalence, viral load or real-world background rates. A more comprehensive
validation dataset should use complete RefSeq FASTA sequences and include
human cfDNA background, realistic fragment-length distributions, sequencing
errors and contamination controls.
