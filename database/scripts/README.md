# Database scripts

This directory contains utilities for preparing and annotating viral reference
metadata. They are separate from the active sample-level workflow. Routine
analysis uses the scripts in `bin/`; database scripts are used when building or
updating annotation resources.

Run commands from the repository root unless a command explicitly changes
directory.

## Script inventory

| Script | Input | Output / purpose |
|---|---|---|
| `extractLocusOrganism.pl` | RefSeq Viral GBFF or GBFF.GZ | Writes a headerless two-column `LOCUS`–`ORGANISM` table to standard output; used to build `database/annotation/locusOrganism.tsv`. |
| `extractViralGroupname.py` | `database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv` | Extracts the `ACCESSION` and `Groupname` columns to `database/annotation/viralGroupname.tsv`. Uses the Python 3 standard library only. |
| `extractGbffMetadata.pl` | RefSeq Viral GBFF or GBFF.GZ | Writes `ACCESSION`, `ORGANISM`, `mol_type`, `host` and `taxon` metadata to standard output. |
| `extractAccessionFromGbff.pl` | GenBank/RefSeq GBFF or GBFF.GZ | Extracts unique accession values from GBFF records. |
| `extractAccessionTaxid.pl` | Accession list and NCBI `nucl_gb.accession2taxid` | Maps accessions to TaxIDs and reports accessions without a mapping. |
| `addGenusToMetadata.pl` | Metadata TSV, NCBI `nodes.dmp`, `names.dmp` | Appends a `Genus` column using NCBI Taxonomy lineage. |
| `filterViralGenomicDnaHost.pl` | Metadata TSV containing `mol_type` | Retains rows where `mol_type` is `genomic DNA`. |
| `extractVirusHostInfo.pl` | Virus-Host DB `virushostdb.tsv` | Extracts the four virus/host ID and name columns into `virushostdb.virus_host.tsv`. |
| `annotateVirusHostDbHost.pl` | Metadata TSV with `taxon`, and `virushostdb.virus_host.tsv` | Appends Virus-Host DB host names; one virus may expand to multiple rows. |
| `filterVirusHostDbHumanHost.pl` | Annotated metadata TSV | Retains rows whose `VirusHostDB_host_name` is exactly `Homo sapiens`. |
| `parseBlastGenus.pl` | BLAST results and NCBI accession/taxonomy files | Standalone genus-level BLAST annotation utility; not called by the current main workflow. |

## Rebuild annotation tables

### `locusOrganism.tsv`

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz

perl database/scripts/extractLocusOrganism.pl \
  viral.1.genomic.gbff.gz \
  > database/annotation/locusOrganism.tsv
```

The output is headerless. Column 1 is the GBFF `LOCUS` identifier and
column 2 is the `ORGANISM` name. Build it from the same RefSeq release used
for the viral BLAST database.

### `viralGroupname.tsv`

The authoritative source is the grouped metadata TSV:

`database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`

The extractor reads the `ACCESSION` and `Groupname` columns and retains the
source accession identifiers:

```bash
python3 database/scripts/extractViralGroupname.py
```

Equivalent explicit paths:

```bash
python3 database/scripts/extractViralGroupname.py \
  --input database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv \
  --output database/annotation/viralGroupname.tsv
```

The generated table has the header `ACCESSION` and `Groupname`. It preserves
the source row order and Groupname labels, and does not strip accession
version suffixes. No Excel file or third-party Python package is required.

### `virushostdb.virus_host.tsv`

See [database/metadata/README.md](../metadata/README.md) for the Virus-Host DB
source, extraction command, output schema and downstream host annotation
steps.

## RefSeq metadata preparation

The scripts support the following staged metadata preparation:

1. `extractGbffMetadata.pl` extracts accession-level metadata from RefSeq
   Viral GBFF records.
2. `addGenusToMetadata.pl` adds Genus using `nodes.dmp` and `names.dmp`.
3. `filterViralGenomicDnaHost.pl` retains genomic DNA records.
4. `annotateVirusHostDbHost.pl` adds Virus-Host DB host names using `taxon`.
5. `filterVirusHostDbHumanHost.pl` retains records whose Virus-Host DB host
   name is `Homo sapiens`.
6. Groupname curation is represented in
   `viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`.
   `extractViralGroupname.py` then produces the two-column runtime lookup.

Example commands for the initial metadata stages (use new output filenames
when rebuilding):

```bash
perl database/scripts/extractGbffMetadata.pl \
  viral.1.genomic.gbff.gz \
  > database/metadata/viral.genomic.metadata.tsv

perl database/scripts/addGenusToMetadata.pl \
  database/metadata/viral.genomic.metadata.tsv \
  /path/to/nodes.dmp \
  /path/to/names.dmp \
  database/metadata/viral.genomic.metadata.with_genus.tsv

perl database/scripts/filterViralGenomicDnaHost.pl \
  database/metadata/viral.genomic.metadata.with_genus.tsv \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.tsv

perl database/scripts/annotateVirusHostDbHost.pl \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.tsv \
  database/metadata/virushostdb.virus_host.tsv \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv

perl database/scripts/filterVirusHostDbHumanHost.pl \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv
```

The checked-in metadata TSVs are snapshots. Their row counts and contents
depend on the RefSeq, NCBI Taxonomy and Virus-Host DB releases used. Record
source release identifiers and retrieval dates when rebuilding them.

## Why some scripts have related functionality

`extractLocusOrganism.pl` and `extractGbffMetadata.pl` both read GenBank flat
files, but they produce different outputs: a lightweight runtime organism
lookup and a richer metadata table. `extractAccessionFromGbff.pl` extracts
accession identifiers only, whereas `extractGbffMetadata.pl` extracts the
primary accession and biological-source annotations.

`addGenusToMetadata.pl` annotates metadata using Taxonomy lineage.
`parseBlastGenus.pl` annotates BLAST results and is not part of the active
sample/Summary pipeline. These scripts have different inputs and outputs and
are retained as separate utilities.
