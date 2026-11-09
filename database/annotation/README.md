# Annotation resources

This directory contains the two annotation tables required by the active
ViruScanNIPT workflow:

- `locusOrganism.tsv` maps RefSeq Viral `LOCUS` identifiers to organism names.
- `viralGroupname.tsv` maps reference accessions to curated `Groupname` labels.

The two tables have different sources and are rebuilt independently.

## 1. `locusOrganism.tsv`

### Purpose

`locusOrganism.tsv` provides the organism name for a viral reference
accession. `bin/parseBlast.pl` uses this table to populate the `viral_name`
column in each `<sample>_blast_out_parsed.txt` file.

### Source

The source is the NCBI RefSeq Viral GenBank flat file:

- File: `viral.1.genomic.gbff.gz`
- Download directory:
  <https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/>

`database/scripts/extractLocusOrganism.pl` extracts the `LOCUS` identifier and
`ORGANISM` name from each GenBank record. It does not infer organism names
from NCBI Taxonomy or Virus-Host DB.

### Rebuild

Run the following commands from the repository root:

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz

perl database/scripts/extractLocusOrganism.pl \
  viral.1.genomic.gbff.gz \
  > database/annotation/locusOrganism.tsv
```

The script accepts compressed `.gbff.gz` files and uncompressed `.gbff`
files. `gzip` must be available when the input is compressed.

### Format

The file is headerless and contains two tab-delimited columns.

| Column | Meaning |
|---|---|
| 1 | `LOCUS` identifier from the GBFF record. |
| 2 | Organism name from the GBFF `ORGANISM` field. |

Example:

```text
NC_029549    Emaravirus tritici
NC_029550    Emaravirus tritici
```

The current repository snapshot contains 19,019 records. The identifiers in
this table are versionless. `bin/parseBlast.pl` removes the version suffix
from each BLAST subject accession before looking up the organism name.

### Update notes

- Rebuild this table from the same RefSeq Viral release used to build the
  viral BLAST database.
- Update the FASTA database and annotation table together to avoid mixing
  release snapshots.
- Record the source release and retrieval date when updating the resource.
- This table supplies organism names. It does not assign Groupname labels
  or host species.

## 2. `viralGroupname.tsv`

### Purpose

`viralGroupname.tsv` maps viral reference accessions to the curated
`Groupname` labels used by `bin/summarizeGroupnameResults.pl`. The Summary
program uses each parsed BLAST `subject` accession to assign query reads to
Groupname categories.

### Authoritative source

The sole source for this mapping is:

`database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`

This tab-delimited metadata table contains the following columns:

```text
ACCESSION  ORGANISM  mol_type  host  taxon  Genus  VirusHostDB_host_name  Groupname
```

Only `ACCESSION` and `Groupname` are extracted. The remaining columns stay
in the source metadata table and are not copied to `viralGroupname.tsv`.

Groupname labels are taken directly from the source `Groupname` column.
They are not recalculated from `ORGANISM`, `Genus`, NCBI Taxonomy files or
Virus-Host DB. No Excel workbook is required.

### Generate `viralGroupname.tsv`

Run the following command from the repository root:

```bash
python3 database/scripts/extractViralGroupname.py
```

The default paths are:

- Input:
  `database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`
- Output: `database/annotation/viralGroupname.tsv`

The paths can also be specified explicitly:

```bash
python3 database/scripts/extractViralGroupname.py \
  --input database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv \
  --output database/annotation/viralGroupname.tsv
```

The script uses only the Python 3 standard library. It:

1. Reads the source TSV using the `ACCESSION` and `Groupname` header names.
2. Extracts these two columns in the original row order.
3. Writes a TSV with the `ACCESSION` and `Groupname` header.
4. Preserves accession identifiers exactly as supplied by the source. It
   does not remove accession version suffixes or rewrite identifiers.
5. Stops with an error if a required column is missing, a value is blank or
   an accession is duplicated.

### Current snapshot and format

The current source table contains 298 data rows, 298 unique accessions and
36 distinct Groupname labels. The generated `viralGroupname.tsv` contains
the same 298 accession-to-Groupname records.

Example:

```text
ACCESSION    Groupname
NC_017993    Papillomavirus
NC_017994    Papillomavirus
```

The `ACCESSION` field is retained as a named column. When matching parsed
BLAST subjects, the Summary program removes accession version suffixes for
comparison. This normalization does not modify the source mapping table.

### Update procedure

When the grouped metadata table is updated:

1. Update
   `viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`.
2. Confirm that `ACCESSION` and `Groupname` are present and populated.
3. Run `extractViralGroupname.py` to regenerate the runtime mapping table.
4. Review the row count, unique accession count and number of Groupname labels.
5. Keep the grouped metadata source and generated mapping table from the same
   annotation snapshot.

Do not edit `viralGroupname.tsv` independently of its source table unless
the corresponding change is also made in the grouped metadata TSV.
