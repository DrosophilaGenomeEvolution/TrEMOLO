#!/usr/bin/env python3
"""Convert legacy whole-assembly TE calls into the public BED artifact."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path


EXPECTED_HEADER = (
    "sseqid",
    "qseqid",
    "pident",
    "size_per",
    "size_el",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


def parse_query_identifier(value: str, line_number: int) -> tuple[str, int, int, str]:
    cluster, separator, location = value.partition("::")
    if not separator or not cluster or not location:
        raise ValueError(
            f"line {line_number}: malformed whole-assembly query identifier: {value!r}"
        )
    coordinate_text, strand_separator, strand = location.rpartition(":")
    chrom, coordinate_separator, interval = coordinate_text.rpartition(":")
    if (
        not strand_separator
        or strand not in {"+", "-"}
        or not coordinate_separator
        or not chrom
        or "-" not in interval
    ):
        raise ValueError(
            f"line {line_number}: malformed whole-assembly query identifier: {value!r}"
        )
    start_text, end_text = interval.split("-", 1)
    try:
        start = int(start_text)
        end = int(end_text)
    except ValueError as error:
        raise ValueError(
            f"line {line_number}: non-integer whole-assembly coordinates: {value!r}"
        ) from error
    if start < 0 or end < start:
        raise ValueError(
            f"line {line_number}: invalid whole-assembly interval: {value!r}"
        )
    return chrom, start, end, cluster


def read_calls(path: Path) -> list[tuple[str, int, int, str]]:
    records = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as error:
            raise ValueError(f"empty whole-assembly TE table: {path}") from error
        if header != EXPECTED_HEADER:
            raise ValueError(
                f"{path}: unexpected whole-assembly TE schema: {header!r}"
            )
        for line_number, fields in enumerate(reader, 2):
            if not fields or not any(fields):
                continue
            if len(fields) != len(EXPECTED_HEADER):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(EXPECTED_HEADER)} columns, "
                    f"got {len(fields)}"
                )
            family = fields[0]
            if not family or "|" in family:
                raise ValueError(
                    f"{path}:{line_number}: invalid TE family for BED output: {family!r}"
                )
            chrom, start, end, cluster = parse_query_identifier(
                fields[1], line_number
            )
            records.append((chrom, start, end, f"{family}|{cluster}"))
    return records


def atomic_write_bed(path: Path, records: list[tuple[str, int, int, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent), text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as handle:
            for record in records:
                handle.write("\t".join(map(str, record)) + "\n")
        os.replace(str(temporary), str(path))
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("calls", type=Path)
    parser.add_argument("output", type=Path)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    atomic_write_bed(args.output, read_calls(args.calls))


if __name__ == "__main__":
    main()
