# Database metadata resources

This directory contains a general Virus-Host DB association table and staged RefSeq-derived viral metadata snapshots. These resources support database construction and annotation. The active sample-level pipeline directly requires only the annotation tables in `database/annotation/`.

## 1. `virushostdb.virus_host.tsv`

### Purpose

`virushostdb.virus_host.tsv` is a four-column subset of the Virus-Host DB virus-host association table. It is used to associate a virus Taxonomy ID with its recorded host names. It is a general virus-host reference table, **not a human-host-only table**.

### Source

The source database is Virus-Host DB, maintained by the Bioinformatics Center, Institute for Chemical Research, Kyoto University:

- Website: <https://www.genome.jp/virushostdb/>
- Download directory: <https://www.genome.jp/ftp/db/virushostdb/>
- Source file: <https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv>
- Release information: `dbrel.txt` in the same download directory.

Virus-Host DB records associations between viruses and hosts using NCBI Taxonomy identifiers. Consult the official [Virus-Host DB README](https://www.genome.jp/ftp/db/virushostdb/README) and release notes for source and curation details.

### How the file is generated

The original `virushostdb.tsv` contains more columns than this project needs. The script `database/scripts/extractVirusHostInfo.pl` extracts the following four columns, in this order:

| Original source column | Output column | Description |
|---|---|---|
| `virus tax id` | `virus_tax_id` | NCBI Taxonomy ID of the virus |
| `virus name` | `virus_name` | Virus name in the source database |
| `host tax id` | `host_tax_id` | NCBI Taxonomy ID of the host |
| `host name` | `host_name` | Host name in the source database |

The output header is created by replacing spaces in the source column names with underscores. Values are copied from the corresponding source fields; the extraction script does not reclassify viruses or hosts.

Run:

```bash
cd database/metadata

wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv

perl ../scripts/extractVirusHostInfo.pl
```

The script reads `virushostdb.tsv` and writes `virushostdb.virus_host.tsv` in the current working directory. It checks that the four required source columns are present. If the upstream header or file format changes, inspect the source header and update the script before rebuilding the table.

### Output format and retention rules

The output is tab-delimited and has this header:

```text
virus_tax_id    virus_name    host_tax_id    host_name
```

- All host associations in the source file are retained, including non-human hosts.
- A virus associated with multiple hosts can occupy multiple rows.
- Blank host TaxID or host-name values are retained as blank fields.
- The extraction script does not deduplicate rows or filter by host species.
- Other source fields, such as evidence, PMID, lineage, disease and sample type, are not included in this four-column subset.

The checked-in snapshot contains 42,742 data rows, 148 exact four-column duplicate rows and 65 rows with blank host TaxID and host-name fields. These are properties of this snapshot, not fixed properties of the upstream database. The row count is not the number of unique virus-host pairs.

## 2. Downstream host annotation

`database/scripts/annotateVirusHostDbHost.pl` uses the `taxon` column in a viral metadata table to look up `virus_tax_id` in `virushostdb.virus_host.tsv`. For each virus TaxID, it appends each unique associated `host_name` as a separate row. If no matching TaxID is found, it appends `NA`. Duplicate host names for the same virus TaxID are collapsed by this annotation step, while the source table itself remains unchanged.

Example:

```bash
perl database/scripts/annotateVirusHostDbHost.pl \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.tsv \
  database/metadata/virushostdb.virus_host.tsv \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv
```

The separate script `database/scripts/filterVirusHostDbHumanHost.pl` can then retain rows whose `VirusHostDB_host_name` is exactly `Homo sapiens`:

```bash
perl database/scripts/filterVirusHostDbHumanHost.pl \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv \
  database/metadata/viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv
```

This human-host restriction is a downstream metadata filter. It must not be applied when constructing the general `virushostdb.virus_host.tsv` resource.

## 3. RefSeq-derived metadata snapshots

The following tab-delimited files are staged outputs from the RefSeq viral metadata preparation workflow. They are reference snapshots, not runtime inputs to the sample-level workflow.

| File | Data rows | Columns / role |
|---|---:|---|
| `viral.genomic.metadata.tsv` | 19,019 | Base RefSeq metadata: `ACCESSION`, `ORGANISM`, `mol_type`, `host`, `taxon`. |
| `viral.genomic.metadata.with_genus.tsv` | 19,019 | Base metadata with `Genus` added from NCBI Taxonomy. |
| `viral.genomic.metadata.with_genus.genomicDNA.tsv` | 9,887 | Genomic DNA records only. |
| `viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv` | 12,013 | Genomic DNA records expanded by Virus-Host DB host associations. |
| `viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv` | 298 | Human-host subset (`VirusHostDB_host_name` = `Homo sapiens`). |
| `viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv` | 298 | Human-host subset with the curated `Groupname` column. |

The checked-in row counts describe the current repository snapshot. They may change when RefSeq, NCBI Taxonomy, Virus-Host DB or the curated Groupname assignments are updated.

`viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv` is the
authoritative tab-delimited source for the runtime Groupname lookup.
`database/scripts/extractViralGroupname.py` extracts its `ACCESSION` and
`Groupname` columns, in source row order, into
`database/annotation/viralGroupname.tsv`. The accession identifiers are
retained as supplied in the source table. No Excel workbook is needed.

For the scripts and commands used to rebuild these stages, see [database/scripts/README.md](../scripts/README.md).

## 4. Updating and citation

When updating Virus-Host DB:

1. Download the latest `virushostdb.tsv` from the official download directory.
2. Record the retrieval date and release identifier from `dbrel.txt`.
3. Rebuild the four-column subset with `extractVirusHostInfo.pl`.
4. Keep the source download and generated table from the same release.
5. Re-run downstream host annotation and review unmatched TaxIDs.

The Virus-Host DB publication is:

Mihara T, Nishimura Y, Shimizu Y, et al. Linking virus genomes with host taxonomy. *Viruses*. 2016;8(3):66. <https://doi.org/10.3390/v8030066>

Check the official database website for the current release, terms of use and citation requirements before redistributing or publishing results based on an updated snapshot.
