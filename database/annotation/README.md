# Annotation resources

This directory contains two accession-level annotation tables used by
ViruScanNIPT. They have different sources and are built independently:

- `locusOrganism.tsv`: RefSeq Viral GenBank `LOCUS` to `ORGANISM` mapping.
- `viralGroupname.tsv`: curated `ACCESSION` to `Groupname` mapping extracted
  from the project's supplied Excel workbook.

## 1. `locusOrganism.tsv`

### Purpose

`locusOrganism.tsv` maps a viral reference sequence identifier to its organism
name. The sample-level parser (`bin/parseBlast.pl`) uses this table to populate
the `viral_name` column in `<sample>_blast_out_parsed.txt`.

### Source

The table is extracted from the NCBI RefSeq Viral GenBank flat file:

- File: `viral.1.genomic.gbff.gz`
- NCBI RefSeq release directory:
  <https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/>

Each GBFF record contains a `LOCUS` line and an `ORGANISM` annotation. The
project script `database/scripts/extractLocusOrganism.pl` reads these two
fields from every record and writes them as a tab-delimited pair. It does not
use NCBI Taxonomy or Virus-Host DB to infer the organism name.

### Rebuild

Run from the repository root:

```bash
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.genomic.gbff.gz

perl database/scripts/extractLocusOrganism.pl \
  viral.1.genomic.gbff.gz \
  > database/annotation/locusOrganism.tsv
```

The script accepts either a compressed `.gbff.gz` file or an uncompressed
`.gbff` file. `gzip` must be available when the input is compressed.

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

The bundled snapshot contains 19,019 records. The identifiers in this table
are versionless. `bin/parseBlast.pl` removes a version suffix from each BLAST
subject accession before looking it up, allowing versioned BLAST accessions to
match this table.

### Update and consistency notes

- Rebuild this table from the same RefSeq Viral release used to build the
  viral BLAST database. Update the FASTA database and annotation table
  together to avoid mixing release snapshots.
- Record the download date or RefSeq release snapshot when updating the
  resource.
- This table supplies organism names only. It does not assign `Groupname`,
  host species or higher taxonomy.

## 2. `viralGroupname.tsv`

### Purpose

`viralGroupname.tsv` maps viral reference accessions to the project's curated
`Groupname` labels. The cohort Summary script
(`bin/summarizeGroupnameResults.pl`) uses this table to convert BLAST
`subject` accessions into Groupname-level counts and detection statuses.

### Source workbook

The source is the Excel workbook supplied for this project:

`viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx`

A copy is retained in this repository at:

`database/annotation/source/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx`

The workbook contains two worksheets:

- `Accession_Grouping`: accession-level annotation records. This is the
  worksheet used to build `viralGroupname.tsv`.
- `Groupname_Summary`: a grouped summary for review. It is not used to create
  the accession-level mapping file.

The `Accession_Grouping` worksheet includes columns such as `ACCESSION`,
`ORGANISM`, `mol_type`, `host`, `taxon`, `Genus`, `VirusHostDB_host_name` and
`Groupname`. Only `ACCESSION` and `Groupname` are extracted. The Groupname
labels are taken directly from this workbook; they are not inferred from
`ORGANISM`, `Genus`, NCBI Taxonomy (`nodes.dmp` / `names.dmp`) or Virus-Host DB.

### How `viralGroupname.tsv` is generated

Run from the repository root:

```bash
python3 -m pip install openpyxl

python3 database/scripts/extractViralGroupname.py \
  --input database/annotation/source/viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx \
  --output database/annotation/viralGroupname.tsv
```

The script performs the following steps:

1. Reads the `Accession_Grouping` worksheet.
2. Locates the `ACCESSION` and `Groupname` columns by their header names.
3. Copies these two values for each data row, preserving the source row order
   and Groupname spelling.
4. Writes a tab-delimited file with the header `ACCESSION` and `Groupname`.
5. Stops with an error if either required value is blank or an accession is
   duplicated, so that an ambiguous mapping is not silently introduced.

The script does not read or recalculate the `Groupname_Summary` worksheet.
It also does not merge, split, rename or otherwise reinterpret Groupname
labels.

### Current table snapshot and format

The supplied workbook's `Accession_Grouping` worksheet contains 298 data rows,
298 unique accessions and 36 distinct Groupname labels. There are no blank
`ACCESSION` or `Groupname` values and no duplicated accessions in this
snapshot. The generated TSV therefore contains 299 lines including its header.

The output is tab-delimited and has a header:

```text
ACCESSION    Groupname
NC_017993    Papillomavirus
NC_017994    Papillomavirus
NC_030297    Torque-Teno virus
```

The source workbook's accession values are retained as supplied. The current
snapshot uses versionless accessions. The Summary script removes a version
suffix from BLAST subject accessions before matching, so a versioned subject
can still match its versionless table entry.

### Updating the table

When the Groupname assignments change:

1. Update or replace the source workbook in
   `database/annotation/source/`.
2. Confirm that the accession-level data remain in the
   `Accession_Grouping` worksheet and that the required columns are named
   `ACCESSION` and `Groupname`.
3. Run `extractViralGroupname.py` to regenerate the TSV.
4. Review the number of records, unique accessions and Groupname labels in the
   script output; compare them with the source workbook.
5. Run the project smoke test and review `unmappedAccessions.tsv` after a new
   cohort Summary run.

Do not use the `locusOrganism.tsv` extraction script or the Virus-Host DB
scripts to generate Groupname values. The Excel workbook is the authoritative
source for the Groupname assignments bundled with this project.
