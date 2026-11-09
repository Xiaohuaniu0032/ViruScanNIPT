# Database scripts

This directory contains utilities for preparing and annotating viral
reference metadata. These scripts are separate from the active sample-level
analysis workflow. Routine sample analysis and cohort Summary use the scripts
in `bin/`. Database scripts are used when building or updating reference
resources.

Run commands from the repository root unless a command explicitly changes
the working directory.

## 1. Script inventory

| Script | Input | Output or purpose |
|---|---|---|
| `extractLocusOrganism.pl` | RefSeq Viral GBFF or GBFF.GZ | Writes a headerless `LOCUS` to `ORGANISM` table to standard output. Used to build `database/annotation/locusOrganism.tsv`. |
| `extractViralGroupname.py` | Grouped viral metadata TSV | Extracts `ACCESSION` and `Groupname` to `database/annotation/viralGroupname.tsv`. Uses the Python 3 standard library. |
| `extractGbffMetadata.pl` | RefSeq Viral GBFF or GBFF.GZ | Writes `ACCESSION`, `ORGANISM`, `mol_type`, `host` and `taxon` metadata to standard output. |
| `extractAccessionFromGbff.pl` | GenBank or RefSeq GBFF, including compressed input | Extracts accession identifiers from GBFF records. |
| `extractAccessionTaxid.pl` | Accession list and NCBI `nucl_gb.accession2taxid` | Maps accessions to TaxIDs and reports accessions without a mapping. |
| `addGenusToMetadata.pl` | Metadata TSV, NCBI `nodes.dmp` and `names.dmp` | Adds a `Genus` column using NCBI Taxonomy lineage. |
| `filterViralGenomicDnaHost.pl` | Metadata TSV containing `mol_type` | Retains rows where `mol_type` is `genomic DNA`. |
| `extractVirusHostInfo.pl` | Virus-Host DB `virushostdb.tsv` | Extracts the four virus and host ID and name columns into `virushostdb.virus_host.tsv`. |
| `annotateVirusHostDbHost.pl` | Metadata TSV with `taxon` and `virushostdb.virus_host.tsv` | Appends Virus-Host DB host names. A virus can expand to multiple rows. |
| `filterVirusHostDbHumanHost.pl` | Annotated metadata TSV | Retains rows whose `VirusHostDB_host_name` is exactly `Homo sapiens`. |
| `parseBlastGenus.pl` | BLAST results and NCBI accession or taxonomy files | Standalone genus-level BLAST annotation utility. It is not called by the current main workflow. |

## 2. Rebuild annotation tables

### 2.1. `locusOrganism.tsv`

Download the RefSeq Viral GenBank flat file and extract the `LOCUS` and
`ORGANISM` fields:

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz

perl database/scripts/extractLocusOrganism.pl \
  viral.1.genomic.gbff.gz \
  > database/annotation/locusOrganism.tsv
```

The output is headerless. Column 1 is the GBFF `LOCUS` identifier and
column 2 is the `ORGANISM` name. Build this table from the same RefSeq release
used for the viral BLAST database.

### 2.2. `viralGroupname.tsv`

The authoritative source is:

`database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`

The extractor reads the `ACCESSION` and `Groupname` columns and retains the
accession identifiers exactly as they appear in the source table.

Run from the repository root:

```bash
python3 database/scripts/extractViralGroupname.py
```

To specify the paths explicitly:

```bash
python3 database/scripts/extractViralGroupname.py \
  --input database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv \
  --output database/annotation/viralGroupname.tsv
```

The generated table has the `ACCESSION` and `Groupname` header. The extractor
preserves source row order and Groupname labels. It does not remove accession
version suffixes. No Excel file or third-party Python package is required.

### 2.3. `virushostdb.virus_host.tsv`

See [database/metadata/README.md](../metadata/README.md) for the Virus-Host DB
source, extraction command, output schema and downstream host annotation
steps.

## 3. RefSeq metadata preparation

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
   `extractViralGroupname.py` then generates the two-column runtime lookup.

Example commands for the initial metadata stages are shown below. Use new
output filenames when rebuilding so that existing snapshots are not
overwritten unintentionally.

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

The checked-in metadata files are snapshots. Their row counts and contents
depend on the RefSeq, NCBI Taxonomy and Virus-Host DB releases used. Record
source release identifiers and retrieval dates when rebuilding them.

## 4. Related script functions

Some scripts read related source data but produce different outputs:

- `extractLocusOrganism.pl` creates a lightweight runtime organism lookup.
- `extractGbffMetadata.pl` creates a richer accession-level metadata table.
- `extractAccessionFromGbff.pl` extracts accession identifiers without
  producing the full metadata table.
- `addGenusToMetadata.pl` annotates metadata using Taxonomy lineage.
- `parseBlastGenus.pl` annotates BLAST results at the genus level and is not
  part of the active sample-level or cohort Summary workflow.

These utilities are retained because their inputs, outputs and purposes
differ.
