# Deprecated and removable scripts

The following six scripts were used by the previous multi-stage Summary workflow.
Their responsibilities have been consolidated into `bin/summarizeGroupnameResults.pl`.
They are temporarily retained under `deprecated/summary/` for result comparison and
should be removed after the new Summary outputs have been validated against the
previous release.

| Previous script | Function now handled by |
|---|---|
| `VSN_06_merge_sample_results.pl` | `summarizeGroupnameResults.pl` reads parsed files directly |
| `VSN_07_annotate_human_dna_virus.pl` | Accession-to-Groupname lookup inside Summary |
| `VSN_08_call_group_status.pl` | Sample × Groupname counts and status inside Summary |
| `VSN_09_group_positive_ranking.pl` | `groupSummary.tsv` |
| `VSN_10_codetection_analysis.pl` | `groupCoDetection.tsv` |
| `VSN_11_group_distribution_stats.pl` | `groupSummary.tsv` and `groupReadCounts.tsv` |

Other housekeeping items:

- `bin/ViroStat_Visualizer.py` has been moved to `tools/viroStatVisualizer.py`; it is
  optional and is no longer launched by the main program.
- The old `bin/cmd.sh` and `db/cmd.sh` command snippets are retained under
  `deprecated/` for reference only. They are not part of the runtime workflow.
- The duplicate `bin/extract_gbff_metadata.pl` was identical to the database copy.
  A single canonical copy is retained as `database/scripts/extractGbffMetadata.pl`.
- `extractLocusOrganism.pl` is retained as a database utility, while the richer
  `extractGbffMetadata.pl` is the recommended source for accession-level metadata.
- The old `bin/` locations for database-construction utilities and intermediate
  accession/taxonomy files have been moved to `database/scripts/` and
  `database/build/`, respectively.

After validating the new Summary output, the six scripts under
`deprecated/summary/` can be deleted.


The deprecated `VSN_07_annotate_human_dna_virus.pl` belongs to the previous
Summary implementation and may refer to a legacy accession-to-Groupname input
file. Do not use that legacy input as the current annotation source. The active
workflow uses `database/annotation/viralGroupname.tsv`, which is generated from
the `Accession_Grouping` worksheet of the bundled Excel workbook; see
`database/annotation/README.md`.
