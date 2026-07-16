#!/usr/bin/env python3
"""Semantic regression check for the historical TrEMOLO test run."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_EXPECTED = HERE / "work_test_expected.json"
DEFAULT_ACTUAL = HERE.parents[2] / "work_test"


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_te_infos(path: Path) -> dict:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))

    if not rows:
        raise AssertionError(f"Empty result file: {path}")

    header, records = rows[0], rows[1:]
    widths = Counter(map(len, records))
    if len(widths) != 1:
        raise AssertionError(f"Inconsistent record widths in {path}: {dict(widths)}")

    chromosomes = Counter(row[0] for row in records)
    families = Counter(row[3].split("|", 1)[0] for row in records)
    event_types = Counter(row[-1] for row in records)
    sources = Counter()
    for row in records:
        if "_INSIDER" in row[13]:
            sources["INSIDER"] += 1
        elif "_OUTSIDER" in row[13]:
            sources["OUTSIDER"] += 1
        else:
            sources["UNKNOWN"] += 1
    calls_with_tsd = sum(row[5] != "NONE" for row in records)

    return {
        "records": len(records),
        "columns_in_header": len(header),
        "columns_in_records": next(iter(widths), 0),
        "chromosomes": dict(sorted(chromosomes.items())),
        "sources": dict(sorted(sources.items())),
        "event_types": dict(sorted(event_types.items())),
        "families": dict(sorted(families.items())),
        "calls_with_tsd": calls_with_tsd,
        "sha256": sha256(path),
    }


def compare(expected: dict, actual: dict, strict: bool) -> list[str]:
    errors = []
    semantic_keys = (
        "records",
        "columns_in_header",
        "columns_in_records",
        "chromosomes",
        "sources",
        "event_types",
        "families",
        "calls_with_tsd",
    )
    for key in semantic_keys:
        if actual[key] != expected[key]:
            errors.append(f"te_infos.{key}: expected {expected[key]!r}, got {actual[key]!r}")
    if strict and actual["sha256"] != expected["sha256"]:
        errors.append(
            "te_infos.sha256: exact file changed "
            f"({expected['sha256']} -> {actual['sha256']})"
        )
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("actual", nargs="?", type=Path, default=DEFAULT_ACTUAL)
    parser.add_argument("--expected", type=Path, default=DEFAULT_EXPECTED)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="also require the byte-for-byte TE_INFOS.bed checksum",
    )
    args = parser.parse_args()

    expected = json.loads(args.expected.read_text())
    errors = []
    te_infos = args.actual / "TE_INFOS.bed"
    if not te_infos.is_file():
        errors.append(f"missing output: {te_infos}")
    else:
        errors.extend(compare(expected["te_infos"], read_te_infos(te_infos), args.strict))

    for relative, expected_lines in expected["key_outputs"].items():
        path = args.actual / relative
        if not path.is_file():
            errors.append(f"missing output: {path}")
            continue
        actual_lines = line_count(path)
        if actual_lines != expected_lines:
            errors.append(
                f"{relative}: expected {expected_lines} lines, got {actual_lines}"
            )

    if errors:
        print("TrEMOLO regression check: FAILED")
        for error in errors:
            print(f"- {error}")
        return 1

    mode = "strict" if args.strict else "semantic"
    print(f"TrEMOLO regression check: PASSED ({mode})")
    print(f"Reference: {expected['name']}")
    for limitation in expected["known_limitations"]:
        print(f"WARNING: {limitation}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
