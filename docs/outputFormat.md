# Summary output formats

All files are tab-delimited text (TSV), regardless of the `.tsv` extension.

## sampleGroupCounts.tsv

One row for each sample × valid Groupname.

- `Sample`: sample identifier.
- `Groupname`: Groupname from the annotation table.
- `DistinctQueryCount`: number of unique query identifiers assigned to the Groupname.
- `Status`: Positive when the count is greater than or equal to the configured cutoff.

## sampleSummary.tsv

One row per sample.

- `TotalDistinctAssignedQueries`: union of query identifiers assigned to all valid Groupnames.
- `GroupnamesWithAtLeastOneRead`: number of Groupnames with at least one assigned query.
- `PositiveGroupCount`: number of Groupnames meeting the positive cutoff.
- `PositiveGroupnames`: semicolon-delimited list, or `NA`.

## groupSummary.tsv

One row per Groupname.

- `PositiveSampleCount`: number of samples meeting the positive cutoff.
- `TotalSampleCount`: number of samples represented by parsed result files.
- `PositiveRatePercent`: positive sample count divided by total sample count.
- `TotalDistinctQueries`: sum of distinct query counts across samples for this Groupname.
- `MeanReadsPerPositiveSample`, `MedianReadsPerPositiveSample`, `Q1`, `Q3`, `Min`, `Max`: distribution among positive samples only.

## groupReadCounts.tsv

One row per positive sample × Groupname. It contains the Groupname, sample identifier and
distinct query count. It is intended for analysis and visualization, not public
supplementary tables containing de-identified aggregate data.

## groupCoDetection.tsv

One row per observed pair of positive Groupnames. A sample with three or more positive
Groupnames contributes to more than one pairwise combination.

## unmappedAccessions.tsv

Accession identifiers present in parsed BLAST results but absent from the Groupname
table. This is a reference annotation quality-control output, not a list of positive
clinical detections.


## sampleCoDetection.tsv

One row per sample with at least two positive Groupnames. The combination is the
alphabetically sorted list of all positive Groupnames in that sample.

## groupCombinationSummary.tsv

One row per exact Groupname combination, ranked by the number of samples containing
that exact combination. This preserves the previous pipeline's exact-combination
summary. `groupCoDetection.tsv` separately reports pairwise co-detection, so a sample
with three or more positive Groupnames contributes to multiple pairwise rows.
