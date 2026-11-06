# Virus-Host DB metadata resource

## 1. File purpose

`virushostdb.virus_host.tsv` is a four-column subset of the Virus-Host DB
virus-host association table. It is used by the database annotation utilities
to associate a virus Taxonomy ID with its recorded host names. It is a general
virus-host reference table, **not a human-host-only table**.

## 2. Source

The source database is Virus-Host DB, maintained by the Bioinformatics Center,
Institute for Chemical Research, Kyoto University:

- Website: <https://www.genome.jp/virushostdb/>
- Download directory: <https://www.genome.jp/ftp/db/virushostdb/>
- Source file: [`virushostdb.tsv`](https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv)
- Release information: `dbrel.txt` in the same download directory.

Virus-Host DB records associations between viruses and hosts using NCBI
Taxonomy identifiers. The database documentation describes the underlying
sources and curation. Consult the official
[Virus-Host DB README](https://www.genome.jp/ftp/db/virushostdb/README) and
release notes for details about a particular snapshot.

## 3. How `virushostdb.virus_host.tsv` is generated

The original `virushostdb.tsv` contains more columns than this project needs.
The script `database/scripts/extractVirusHostInfo.pl` extracts the following
four columns in this order:

| Original source column | Output column | Description |
|---|---|---|
| `virus tax id` | `virus_tax_id` | NCBI Taxonomy ID of the virus |
| `virus name` | `virus_name` | Virus name in the source database |
| `host tax id` | `host_tax_id` | NCBI Taxonomy ID of the host |
| `host name` | `host_name` | Host name in the source database |

The output header is created by replacing spaces in the source column names
with underscores. The values are copied from the corresponding source fields;
the extraction script does not reclassify viruses or hosts.

### Rebuild

The extraction script expects the downloaded source file to be named
`virushostdb.tsv` in the current working directory. It writes
`virushostdb.virus_host.tsv` to that same directory. Run:

```bash
cd database/metadata

wget -c https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv

perl ../scripts/extractVirusHostInfo.pl
```

The output will be:

```text
database/metadata/virushostdb.virus_host.tsv
```

The script checks that all four source columns are present. If the upstream
header or file format changes, inspect the source header and update the script
before rebuilding the table.

## 4. Output format and retention rules

The output is a tab-delimited table with this header:

```text
virus_tax_id    virus_name    host_tax_id    host_name
```

Important behavior:

- All host associations in the source file are retained, including non-human
  hosts.
- A virus associated with multiple hosts can occupy multiple rows.
- Blank host TaxID or host-name values are retained as blank fields.
- The extraction script does not deduplicate rows or filter by host species.
- Other source fields, such as evidence, PMID, lineage, disease and sample
  type, are not included in this four-column subset.

The bundled snapshot contains 42,742 data rows. It includes 148 exact
four-column duplicate rows and 65 rows with blank host TaxID and host-name
fields. These are properties of this extracted snapshot; the row count should
not be interpreted as the number of unique virus-host pairs.

## 5. Downstream use in this repository

`database/scripts/annotateVirusHostDbHost.pl` uses the `taxon` column in a
viral metadata table to look up `virus_tax_id` in this resource. For each
virus TaxID, it appends each unique associated `host_name` as a separate row.
If no matching TaxID is found, it appends `NA`. Duplicate host names for the
same virus TaxID are collapsed by this annotation step, while the source table
itself remains unchanged.

Example:

```bash
perl database/scripts/annotateVirusHostDbHost.pl \
  viral.genomic.metadata.with_genus.genomicDNA.tsv \
  database/metadata/virushostdb.virus_host.tsv \
  viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv
```

The separate script `database/scripts/filterVirusHostDbHumanHost.pl` can then
retain rows whose `VirusHostDB_host_name` is exactly `Homo sapiens`:

```bash
perl database/scripts/filterVirusHostDbHumanHost.pl \
  viral.genomic.metadata.with_genus.genomicDNA.host_annotated.tsv \
  viral.genomic.metadata.with_genus.genomicDNA.host_human.tsv
```

This human-host restriction is a downstream analysis filter. It must not be
applied when constructing the general `virushostdb.virus_host.tsv` resource.

## 6. Updating and citation

When updating this resource:

1. Download the latest `virushostdb.tsv` from the official download directory.
2. Record the retrieval date and release identifier from `dbrel.txt`.
3. Rebuild the four-column subset with `extractVirusHostInfo.pl`.
4. Keep the source download and generated table from the same release.
5. Re-run downstream host annotation and review unmatched TaxIDs.

The Virus-Host DB publication is:

Mihara T, Nishimura Y, Shimizu Y, et al. Linking virus genomes with host
taxonomy. *Viruses*. 2016;8(3):66.
<https://doi.org/10.3390/v8030066>

Check the official database website for the current release, terms of use and
citation requirements before redistributing or publishing results based on a
new snapshot.
