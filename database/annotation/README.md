# Annotation resources

This directory contains the two annotation tables used by the active
ViruScanNIPT workflow:

- `locusOrganism.tsv`: RefSeq Viral `LOCUS` to `ORGANISM` lookup.
- `viralGroupname.tsv`: curated `ACCESSION` to `Groupname` lookup.

The two files have different sources and are rebuilt independently.

## 1. `locusOrganism.tsv`

### Purpose

`locusOrganism.tsv` maps a viral reference sequence identifier to its organism
name. `bin/parseBlast.pl` uses this table to populate the `viral_name` column
in each `<sample>_blast_out_parsed.txt` file.

### Source

The source is the NCBI RefSeq Viral GenBank flat file:

- File: `viral.1.genomic.gbff.gz`
- Download directory:
  <https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/>

The project script `database/scripts/extractLocusOrganism.pl` extracts the
`LOCUS` identifier and `ORGANISM` name from each GenBank record. It does not
infer organism names from NCBI Taxonomy or Virus-Host DB.

### Rebuild

Run from the repository root:

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz

perl database/scripts/extractLocusOrganism.pl \
  viral.1.genomic.gbff.gz \
  > database/annotation/locusOrganism.tsv
```

The script accepts compressed `.gbff.gz` and uncompressed `.gbff` input.
`gzip` must be available for compressed input.

### Format

The file is headerless and contains two tab-delimited columns:

| Column | Meaning |
|---|---|
| 1 | `LOCUS` identifier from the GBFF record |
| 2 | Organism name from the GBFF `ORGANISM` field |

Example:

```text
NC_029549    Emaravirus tritici
NC_029550    Emaravirus tritici
```

The checked-in snapshot contains 19,019 records. The identifiers are
versionless. `bin/parseBlast.pl` removes the version suffix from each BLAST
subject accession before looking up the organism name.

### Update notes

- Rebuild this table from the same RefSeq Viral release used to build the
  viral BLAST database.
- Update the FASTA database and annotation table together to avoid mixing
  release snapshots.
- Record the source release and retrieval date when updating the table.
- This table supplies organism names only; it does not assign Groupname or
  host species.

## 2. `viralGroupname.tsv`

### Purpose

`viralGroupname.tsv` maps viral reference accessions to the curated `Groupname`
labels used by `bin/summarizeGroupnameResults.pl`. The Summary script uses the
`subject` accession in parsed BLAST results to assign each query read to a
Groupname.

### Authoritative source

The sole source for this mapping is:

`database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`

This is the curated, tab-delimited metadata table in `database/metadata/`.
It contains the following columns:

```text
ACCESSION  ORGANISM  mol_type  host  taxon  Genus  VirusHostDB_host_name  Groupname
```

Only the `ACCESSION` and `Groupname` columns are extracted. The other columns
remain in the source metadata table and are not copied into
`viralGroupname.tsv`.

The Groupname values are taken directly from the `Groupname` column. They are
not recalculated from `ORGANISM`, `Genus`, NCBI Taxonomy files or
Virus-Host DB. No Excel workbook is required or used.

### Generate `viralGroupname.tsv`

Run from the repository root:

```bash
python3 database/scripts/extractViralGroupname.py
```

The default input and output paths are:

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

1. Reads the source TSV using the header names `ACCESSION` and `Groupname`.
2. Extracts these two columns in the original row order.
3. Writes a two-column TSV with the header `ACCESSION` and `Groupname`.
4. Preserves the accession identifier as supplied in the source. It does not
   remove accession version suffixes or rewrite accession identifiers.
5. Stops with an error if either required column is absent, a value is blank,
   or an `ACCESSION` occurs more than once.

### Current snapshot and format

The current source table contains 298 data rows, 298 unique `ACCESSION`
values and 36 distinct Groupname labels. The generated `viralGroupname.tsv`
contains the same 298 accession-to-Groupname records.

The output format is:

```text
ACCESSION    Groupname
NC_017993    Papillomavirus
NC_017994    Papillomavirus
```

The `ACCESSION` field is retained as a named column. Accession version
normalization, when needed for matching parsed BLAST subjects, is performed
by the Summary program and does not modify this source annotation table.

### Update notes

When the grouped metadata table is revised:

1. Update `viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv`.
2. Confirm that `ACCESSION` and `Groupname` are present and populated.
3. Run `extractViralGroupname.py` to regenerate the mapping table.
4. Review the number of rows, unique accessions and Groupname labels.
5. Keep the grouped metadata source and generated mapping table from the same
   annotation snapshot.

Do not edit `viralGroupname.tsv` independently of its source table unless the
change is also reflected in the grouped metadata TSV.
