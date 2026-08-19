#!/usr/bin/env python3
"""Prepare the BED inputs used by legacy INSIDER frequency estimation.

The historical ``FREQ_INSIDERv2`` rule derived four interval tables from the
``INSERTION.csv`` candidate table inside one large shell block.  This module
keeps the same candidate order, coordinates, and duplicate handling while
validating the input before writing any output.  Unlike the historical awk
commands, the CSV header is never interpreted as a candidate.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


REQUIRED_COLUMNS = ("sseqid", "qseqid")


def non_negative_int(value):
    """Return an integer suitable for an argparse non-negative option."""

    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be greater than or equal to zero")
    return parsed


def build_parser():
    """Build the command-line parser without parsing arguments at import."""

    parser = argparse.ArgumentParser(
        description="Prepare intervals for INSIDER allele-frequency estimation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "candidates", type=Path, help="tab-delimited INSERTION.csv candidate table"
    )
    parser.add_argument(
        "insertion_bed", type=Path, help="output candidate insertion BED"
    )
    parser.add_argument(
        "outer_flanks_bed", type=Path, help="output outer-flank BED"
    )
    parser.add_argument(
        "inner_flanks_bed", type=Path, help="output inner-flank BED"
    )
    parser.add_argument(
        "depth_flanks_bed", type=Path, help="output depth-sampling BED"
    )
    parser.add_argument(
        "--outer-margin",
        type=non_negative_int,
        default=100,
        help="distance from each insertion boundary to the outer flank",
    )
    parser.add_argument(
        "--inner-margin",
        type=non_negative_int,
        default=5,
        help="distance inside each insertion boundary to the inner flank",
    )
    parser.add_argument(
        "--depth-margin",
        type=non_negative_int,
        default=30,
        help="distance from each insertion boundary to sample read depth",
    )
    return parser


def parse_qseqid(qseqid, path, line_number):
    """Return ``event_id, chromosome, start, end`` from one query name."""

    fields = qseqid.split(":")
    if len(fields) < 6:
        raise ValueError(
            "{}:{}: malformed qseqid (expected at least six ':'-delimited "
            "fields): {}".format(path, line_number, qseqid)
        )

    event_id = fields[0]
    chromosome = fields[4]
    coordinates = fields[5].split("-")
    if not event_id or not chromosome or len(coordinates) != 2 or not all(coordinates):
        raise ValueError(
            "{}:{}: malformed qseqid event, chromosome, or coordinate field: {}".format(
                path, line_number, qseqid
            )
        )

    try:
        start, end = (int(coordinate) for coordinate in coordinates)
    except ValueError:
        raise ValueError(
            "{}:{}: qseqid coordinates must be integers: {}".format(
                path, line_number, fields[5]
            )
        )
    if start < 0 or end < 0:
        raise ValueError(
            "{}:{}: qseqid coordinates must be non-negative: {}".format(
                path, line_number, fields[5]
            )
        )
    if start >= end:
        raise ValueError(
            "{}:{}: qseqid start must be smaller than end: {}".format(
                path, line_number, fields[5]
            )
        )
    return event_id, chromosome, start, end


def read_candidates(path):
    """Read and validate INSIDER insertion candidates in their source order."""

    candidates = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("{}: missing candidate header".format(path))
        missing = [name for name in REQUIRED_COLUMNS if name not in reader.fieldnames]
        if missing:
            raise ValueError(
                "{}: missing candidate column(s): {}".format(path, ", ".join(missing))
            )

        for row in reader:
            if not row or all(value in (None, "") for value in row.values()):
                continue
            sseqid = row.get("sseqid")
            qseqid = row.get("qseqid")
            if not sseqid or not qseqid:
                raise ValueError(
                    "{}:{}: empty sseqid or qseqid".format(path, reader.line_num)
                )
            event_id, chromosome, start, end = parse_qseqid(
                qseqid, path, reader.line_num
            )
            candidates.append(
                {
                    "family": sseqid,
                    "event_id": event_id,
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "line_number": reader.line_num,
                }
            )
    return candidates


def bed_line(chromosome, start, end, name=None):
    """Format one validated BED interval."""

    if start < 0 or end < 0 or start >= end:
        raise ValueError(
            "invalid BED interval {}:{}-{}".format(chromosome, start, end)
        )
    fields = [chromosome, str(start), str(end)]
    if name is not None:
        fields.append(name)
    return "\t".join(fields) + "\n"


def append_probe(lines, chromosome, start, end):
    """Append an in-bounds probe, omitting coordinates left of a contig."""

    if start < 0 or end <= 0:
        return
    lines.append(bed_line(chromosome, start, end))


def prepare_intervals(candidates, outer_margin=100, inner_margin=5, depth_margin=30):
    """Return the four legacy interval tables as deterministic strings."""

    insertion_lines = []
    outer_lines = []
    inner_lines = []
    depth_lines = []

    for candidate in candidates:
        chromosome = candidate["chromosome"]
        start = candidate["start"]
        end = candidate["end"]
        name = "{}|{}".format(candidate["family"], candidate["event_id"])

        try:
            insertion_lines.append(bed_line(chromosome, start, end, name))
            append_probe(
                outer_lines,
                chromosome,
                start - outer_margin,
                start - outer_margin + 1,
            )
            append_probe(
                outer_lines,
                chromosome,
                end + outer_margin,
                end + outer_margin + 1,
            )
            append_probe(
                inner_lines,
                chromosome,
                start + inner_margin,
                start + inner_margin + 1,
            )
            append_probe(
                inner_lines,
                chromosome,
                end - inner_margin,
                end - inner_margin + 1,
            )
            append_probe(
                depth_lines,
                chromosome,
                start - depth_margin,
                start - depth_margin + 1,
            )
            append_probe(
                depth_lines,
                chromosome,
                end + depth_margin - 1,
                end + depth_margin,
            )
        except ValueError as error:
            raise ValueError(
                "candidate line {}: {}".format(candidate["line_number"], error)
            )

    return (
        "".join(insertion_lines),
        "".join(outer_lines),
        "".join(inner_lines),
        "".join(depth_lines),
    )


def write_text(path, content):
    """Create a parent directory and write LF-terminated UTF-8 text."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(content)


def main(argv=None):
    """Run the command-line program."""

    args = build_parser().parse_args(argv)
    candidates = read_candidates(args.candidates)
    outputs = prepare_intervals(
        candidates,
        outer_margin=args.outer_margin,
        inner_margin=args.inner_margin,
        depth_margin=args.depth_margin,
    )
    for path, content in zip(
        (
            args.insertion_bed,
            args.outer_flanks_bed,
            args.inner_flanks_bed,
            args.depth_flanks_bed,
        ),
        outputs,
    ):
        write_text(path, content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
