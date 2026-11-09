#!/usr/bin/env python3
"""Build viralGroupname.tsv from the curated grouped RefSeq metadata TSV."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

DEFAULT_INPUT = Path(
    "database/metadata/"
    "viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.tsv"
)
DEFAULT_OUTPUT = Path("database/annotation/viralGroupname.tsv")
REQUIRED_COLUMNS = ("ACCESSION", "Groupname")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract ACCESSION and Groupname from the curated grouped viral "
            "metadata TSV and write database/annotation/viralGroupname.tsv."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help=f"Source grouped metadata TSV (default: {DEFAULT_INPUT})",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Output mapping TSV (default: {DEFAULT_OUTPUT})",
    )
    args = parser.parse_args()

    if not args.input.is_file():
        parser.error(f"Input TSV does not exist: {args.input}")

    records: list[tuple[str, str]] = []
    seen_accessions: set[str] = set()

    try:
        with args.input.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                parser.error(f"Input TSV is empty: {args.input}")

            missing = [
                column for column in REQUIRED_COLUMNS
                if column not in reader.fieldnames
            ]
            if missing:
                parser.error(
                    "Required column(s) missing from input TSV: "
                    + ", ".join(missing)
                )

            for line_number, row in enumerate(reader, start=2):
                accession_value = row.get("ACCESSION")
                groupname_value = row.get("Groupname")

                accession = (
                    "" if accession_value is None else accession_value.strip()
                )
                groupname = (
                    "" if groupname_value is None else groupname_value.strip()
                )

                if not accession or not groupname:
                    parser.error(
                        f"Blank ACCESSION or Groupname at input line "
                        f"{line_number}; please check {args.input}."
                    )

                if accession in seen_accessions:
                    parser.error(
                        f"Duplicated ACCESSION '{accession}' at input line "
                        f"{line_number}."
                    )

                # Keep the ACCESSION identifier as supplied by the source TSV.
                # Do not remove RefSeq version suffixes or rewrite accession IDs.
                seen_accessions.add(accession)
                records.append((accession, groupname))

    except (OSError, UnicodeError, csv.Error) as exc:
        parser.error(f"Cannot read input TSV {args.input}: {exc}")

    if not records:
        parser.error(f"No data rows found in input TSV: {args.input}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    try:
        with args.output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(
                handle, delimiter="\t", lineterminator="\n"
            )
            writer.writerow(REQUIRED_COLUMNS)
            writer.writerows(records)
    except OSError as exc:
        parser.error(f"Cannot write output TSV {args.output}: {exc}")

    group_count = len({groupname for _, groupname in records})
    print(f"Source TSV       : {args.input}")
    print(f"Output TSV       : {args.output}")
    print(f"Data rows        : {len(records)}")
    print(f"Unique ACCESSION : {len(seen_accessions)}")
    print(f"Groupname labels : {group_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
