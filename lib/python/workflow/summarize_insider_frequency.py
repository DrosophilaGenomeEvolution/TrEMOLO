#!/usr/bin/env python3
"""Summarize INSIDER breakpoint depths into the legacy frequency table.

The historical ``FREQ_INSIDERv2`` shell block sampled one base on either
side of every candidate insertion and combined those depths with deletion
evidence.  This module preserves its output and arithmetic while making the
inputs explicit and validating them before replacing the destination file.

In particular, an absent flank has depth zero, the mean uses integer floor
division, and the percentage ratio is truncated to six decimal places before
being multiplied by 100 and printed with four decimal places.  The last detail
also applies when deletion evidence makes read support negative.
"""

from __future__ import print_function

import argparse
import csv
import os
import tempfile
from pathlib import Path


HEADER = (
    "chrom",
    "position",
    "total_depth",
    "depth_empty_site",
    "read_support",
    "read_support_percent",
    "TE",
    "info_TE",
)


def non_negative_int(value):
    """Return an integer suitable for a non-negative CLI option."""

    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be greater than or equal to zero")
    return parsed


def build_parser():
    """Build the command-line parser without parsing at import time."""

    parser = argparse.ArgumentParser(
        description="Summarize INSIDER depth and deletion evidence",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "insertions", type=Path, help="four-column INSERTION_TE.bed"
    )
    parser.add_argument(
        "deletions", type=Path, help="five-column DEL_NB.bed"
    )
    parser.add_argument(
        "depths", type=Path, help="chromosome:position depth table"
    )
    parser.add_argument("output", type=Path, help="DEPTH_TE_INSIDER.csv output")
    parser.add_argument(
        "--depth-margin",
        type=non_negative_int,
        default=30,
        help="distance outside each insertion boundary sampled for depth",
    )
    return parser


def parse_non_negative(value, path, line_number, label):
    """Parse a non-negative integer with source-aware diagnostics."""

    try:
        parsed = int(value)
    except ValueError:
        raise ValueError(
            "{}:{}: {} must be an integer: {}".format(
                path, line_number, label, value
            )
        )
    if parsed < 0:
        raise ValueError(
            "{}:{}: {} must be non-negative: {}".format(
                path, line_number, label, value
            )
        )
    return parsed


def read_insertions(path):
    """Read candidate BED records in order, retaining duplicate candidates."""

    candidates = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if not fields or all(field == "" for field in fields):
                continue
            if len(fields) < 4:
                raise ValueError(
                    "{}:{}: expected at least four tab-delimited columns".format(
                        path, line_number
                    )
                )

            chromosome, start_text, end_text, name = fields[:4]
            if not chromosome:
                raise ValueError(
                    "{}:{}: chromosome must not be empty".format(path, line_number)
                )
            start = parse_non_negative(
                start_text, path, line_number, "insertion start"
            )
            end = parse_non_negative(end_text, path, line_number, "insertion end")
            if start >= end:
                raise ValueError(
                    "{}:{}: insertion start must be smaller than end".format(
                        path, line_number
                    )
                )
            if not name or "|" not in name or not name.split("|", 1)[0]:
                raise ValueError(
                    "{}:{}: candidate name must have the form TE|event".format(
                        path, line_number
                    )
                )

            candidates.append((chromosome, start, end, name))
    return candidates


def read_deletion_counts(path):
    """Read deletion support keyed by the complete ``TE|event`` name."""

    counts = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if not fields or all(field == "" for field in fields):
                continue
            if len(fields) < 5:
                raise ValueError(
                    "{}:{}: expected at least five tab-delimited columns".format(
                        path, line_number
                    )
                )
            name = fields[3]
            if not name:
                raise ValueError(
                    "{}:{}: deletion candidate name must not be empty".format(
                        path, line_number
                    )
                )
            count = parse_non_negative(
                fields[4], path, line_number, "deletion count"
            )
            if name in counts and counts[name] != count:
                raise ValueError(
                    "{}:{}: conflicting deletion counts for {}".format(
                        path, line_number, name
                    )
                )
            counts[name] = count
    return counts


def read_depths(path):
    """Read one-based samtools depth values keyed by ``chromosome:position``."""

    depths = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if not fields or all(field == "" for field in fields):
                continue
            if len(fields) < 2:
                raise ValueError(
                    "{}:{}: expected a position key and depth".format(
                        path, line_number
                    )
                )
            key = fields[0]
            chromosome, separator, position_text = key.rpartition(":")
            if not separator or not chromosome:
                raise ValueError(
                    "{}:{}: malformed depth position key: {}".format(
                        path, line_number, key
                    )
                )
            position = parse_non_negative(
                position_text, path, line_number, "depth position"
            )
            if position == 0:
                raise ValueError(
                    "{}:{}: depth position must use one-based coordinates".format(
                        path, line_number
                    )
                )
            depth = parse_non_negative(fields[1], path, line_number, "depth")
            normalized_key = "{}:{}".format(chromosome, position)
            if normalized_key in depths and depths[normalized_key] != depth:
                raise ValueError(
                    "{}:{}: conflicting depths for {}".format(
                        path, line_number, normalized_key
                    )
                )
            depths[normalized_key] = depth
    return depths


def truncate_division(numerator, denominator):
    """Divide integers with Bash's truncation-toward-zero semantics."""

    if denominator == 0:
        raise ZeroDivisionError("integer division by zero")
    quotient = abs(numerator) // abs(denominator)
    if (numerator < 0) != (denominator < 0):
        return -quotient
    return quotient


def legacy_percentage(support, total_depth):
    """Format the percentage exactly like the historical shell expression."""

    if total_depth == 0:
        return "0.0000"

    # Historical Bash arithmetic evaluated this expression before printf:
    #   ((10**6 * support / total_depth) * 100)e-6
    # Integer division therefore truncated the ratio at 1e-6.  After the
    # multiplication and scientific exponent, the value has exactly four
    # decimal places, so constructing it as an integer avoids float drift.
    scaled_ratio = truncate_division(10 ** 6 * support, total_depth)
    sign = "-" if scaled_ratio < 0 else ""
    magnitude = abs(scaled_ratio)
    return "{}{}.{:04d}".format(sign, magnitude // 10000, magnitude % 10000)


def summarize(candidates, deletion_counts, depths, depth_margin=30):
    """Return output rows while preserving candidate order and duplicates."""

    rows = []
    for chromosome, start, end, name in candidates:
        left_key = "{}:{}".format(chromosome, start - depth_margin + 1)
        right_key = "{}:{}".format(chromosome, end + depth_margin)
        left_depth = depths.get(left_key, 0)
        right_depth = depths.get(right_key, 0)
        total_depth = (left_depth + right_depth) // 2
        family = name.split("|", 1)[0]

        if total_depth == 0:
            empty_depth = 0
            support = 0
        else:
            empty_depth = deletion_counts.get(name, 0)
            support = total_depth - empty_depth

        rows.append(
            (
                chromosome,
                str(start),
                str(total_depth),
                str(empty_depth),
                str(support),
                legacy_percentage(support, total_depth),
                family,
                name,
            )
        )
    return rows


def render(rows):
    """Render the complete LF-terminated TSV output."""

    lines = ["\t".join(HEADER)]
    lines.extend("\t".join(row) for row in rows)
    return "\n".join(lines) + "\n"


def atomic_write(path, content):
    """Replace ``path`` atomically after writing complete UTF-8 content."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=str(path.parent), prefix=".{}.".format(path.name), suffix=".tmp"
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, str(path))
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def main(argv=None):
    """Run the command-line program."""

    args = build_parser().parse_args(argv)
    candidates = read_insertions(args.insertions)
    deletion_counts = read_deletion_counts(args.deletions)
    depths = read_depths(args.depths)
    rows = summarize(
        candidates,
        deletion_counts,
        depths,
        depth_margin=args.depth_margin,
    )
    atomic_write(args.output, render(rows))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
