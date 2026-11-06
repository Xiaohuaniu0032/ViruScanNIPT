#!/usr/bin/env python3
"""Extract ACCESSION and Groupname from the curated Excel annotation workbook."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

try:
    from openpyxl import load_workbook
except ImportError as exc:
    print(
        "[ERROR] Missing Python dependency. Install openpyxl with "
        "'python3 -m pip install openpyxl'.",
        file=sys.stderr,
    )
    raise SystemExit(1) from exc


DEFAULT_INPUT = Path(
    "database/annotation/source/"
    "viral.genomic.metadata.with_genus.genomicDNA.host_human.grouped.xlsx"
)
DEFAULT_OUTPUT = Path("database/annotation/viralGroupname.tsv")
SHEET_NAME = "Accession_Grouping"
REQUIRED_COLUMNS = ("ACCESSION", "Groupname")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build viralGroupname.tsv from the Accession_Grouping worksheet "
            "of the curated ViruScanNIPT annotation workbook."
        )
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT,
                        help=f"Source Excel workbook (default: {DEFAULT_INPUT})")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                        help=f"Output TSV (default: {DEFAULT_OUTPUT})")
    args = parser.parse_args()

    if not args.input.is_file():
        parser.error(f"Input workbook does not exist: {args.input}")

    workbook = load_workbook(args.input, read_only=True, data_only=True)
    if SHEET_NAME not in workbook.sheetnames:
        parser.error(
            f"Worksheet '{SHEET_NAME}' not found. Available sheets: "
            + ", ".join(workbook.sheetnames)
        )

    worksheet = workbook[SHEET_NAME]
    rows = worksheet.iter_rows(values_only=True)
    try:
        header = next(rows)
    except StopIteration:
        parser.error(f"Worksheet '{SHEET_NAME}' is empty")

    header_index = {
        str(value).strip(): index
        for index, value in enumerate(header)
        if value is not None and str(value).strip()
    }
    missing = [name for name in REQUIRED_COLUMNS if name not in header_index]
    if missing:
        parser.error("Required column(s) missing: " + ", ".join(missing))

    accession_index = header_index["ACCESSION"]
    groupname_index = header_index["Groupname"]

    records: list[tuple[str, str]] = []
    seen_accessions: set[str] = set()

    for excel_row, row in enumerate(rows, start=2):
        accession_value = row[accession_index] if accession_index < len(row) else None
        groupname_value = row[groupname_index] if groupname_index < len(row) else None

        accession = "" if accession_value is None else str(accession_value).strip()
        groupname = "" if groupname_value is None else str(groupname_value).strip()

        if not accession or not groupname:
            parser.error(
                f"Blank ACCESSION or Groupname at worksheet row {excel_row}; "
                "please check the source workbook."
            )
        if accession in seen_accessions:
            parser.error(
                f"Duplicated ACCESSION '{accession}' at worksheet row {excel_row}."
            )

        seen_accessions.add(accession)
        records.append((accession, groupname))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(REQUIRED_COLUMNS)
        writer.writerows(records)

    workbook.close()

    group_count = len({groupname for _, groupname in records})
    print(f"Source workbook : {args.input}")
    print(f"Worksheet       : {SHEET_NAME}")
    print(f"Output          : {args.output}")
    print(f"Data rows       : {len(records)}")
    print(f"Unique accessions: {len(seen_accessions)}")
    print(f"Groupname labels: {group_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
