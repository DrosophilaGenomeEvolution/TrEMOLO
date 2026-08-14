#!/usr/bin/env python3
"""Summarize legacy OUTSIDER read evidence into frequency tables.

This is a Python 3 replacement for the two AWK programs in the historical
``FREQUENCEv2`` rule.  The compatibility behaviour is intentional:

* records are ordered by the composite ``event:read:evidence`` key before
  they are counted;
* only the first record for an event/read pair contributes, so an ``E``
  record can mask insertion or clipped-read evidence that sorts after it;
* candidates are grouped by event identifier alone and only ``INS`` and
  ``DEL`` event types are emitted; and
* numbers use AWK's default six-significant-digit output representation.

These quirks must remain isolated here until a corrected frequency model is
introduced deliberately.  In particular, changing the evidence precedence
would change existing TrEMOLO results.
"""

from __future__ import print_function

import argparse
import re
import shutil
from pathlib import Path


NUMBER_PREFIX = re.compile(
    r"^[\t ]*([+-]?(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[eE][+-]?[0-9]+)?)"
)


def awk_number(value):
    """Return AWK-like numeric conversion for a whitespace field."""

    match = NUMBER_PREFIX.match(value)
    if match is None:
        return 0.0
    return float(match.group(1))


def awk_format(value):
    """Format a number like the default numeric output of mawk/gawk."""

    return format(float(value), ".6g")


def read_count_records(path, pre_sorted=False):
    """Yield validated evidence records in legacy aggregation order."""

    records = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split()
            if len(fields) < 17:
                raise ValueError(
                    "{}:{}: expected at least 17 fields, found {}".format(
                        path, line_number, len(fields)
                    )
                )
            if pre_sorted:
                yield fields
            else:
                # The shell pipeline sorted the entire prefixed record, not
                # only the prefix. Keeping the original line as a tie-breaker
                # retains that detail for duplicated evidence rows.
                sort_record = "{}:{}:{} {}".format(
                    fields[2], fields[1], fields[0], line
                )
                records.append((sort_record, fields))
    if pre_sorted:
        return
    records.sort(key=lambda record: record[0])
    for _, fields in records:
        yield fields


def new_event(fields):
    """Initialize the state used by the historical frequency state machine."""

    return {
        "total": 0,
        "deleted": 0,
        "inserted": 0,
        "insertion_only": 0,
        "clipped": 0,
        "sv_type": fields[16],
        "read_support": awk_number(fields[4]),
        "event_id": fields[2],
        "te": fields[3],
        "position": fields[9],
        "current_read": "",
    }


def consume_record(state, fields):
    """Apply one legacy evidence row to an event state."""

    evidence_type = fields[0]
    read_id = fields[1]
    if evidence_type != "E" and state["current_read"] != read_id:
        if state["sv_type"] == "INS":
            state["inserted"] += 1
        if state["sv_type"] == "DEL":
            state["deleted"] += 1

        state["current_read"] = read_id
        state["total"] += 1
        if evidence_type in ("S", "H"):
            state["clipped"] += 1
        elif evidence_type == "I":
            state["insertion_only"] += 1
    elif evidence_type == "E" and state["current_read"] != read_id:
        state["total"] += 1
        state["current_read"] = read_id


def finalize_event(state):
    """Return one legacy output row, or ``None`` for ignored event types."""

    sv_type = state["sv_type"]
    support = state["read_support"]
    total = state["total"]

    if sv_type == "INS":
        inserted = max(state["inserted"], support)
        insertion_only = state["insertion_only"]
        if state["event_id"].startswith("sniffles") and insertion_only < support:
            insertion_only = support
        if insertion_only == 0:
            insertion_only = support
        if total < inserted or total == 0:
            total = inserted

        without_clipped = total - state["clipped"]
        if without_clipped < insertion_only or without_clipped == 0:
            without_clipped = total

        return [
            state["position"],
            state["event_id"],
            state["te"],
            awk_format(insertion_only),
            awk_format(state["clipped"]),
            awk_format(without_clipped),
            awk_format(inserted),
            awk_format(total),
            awk_format((insertion_only / without_clipped) * 100.0),
            awk_format((inserted / total) * 100.0),
            sv_type,
        ]

    if sv_type == "DEL":
        deleted = max(state["deleted"], support)
        if total < deleted or total == 0:
            total = deleted
        reference_reads = total - deleted
        frequency = (reference_reads / total) * 100.0
        return [
            state["position"],
            state["event_id"],
            state["te"],
            awk_format(deleted),
            ".",
            awk_format(total),
            awk_format(reference_reads),
            awk_format(total),
            awk_format(frequency),
            awk_format(frequency),
            sv_type,
        ]

    # HARD and SOFT identifiers do not encode INS/DEL in field 17 in the
    # legacy counter and were silently absent from its frequency output.
    return None


def summarize_counts(read_counts, pre_sorted=False):
    """Return raw legacy frequency rows from ``COUNT_READS.txt``."""

    rows = []
    state = None
    for fields in read_count_records(read_counts, pre_sorted=pre_sorted):
        if state is None or state["event_id"] != fields[2]:
            if state is not None:
                row = finalize_event(state)
                if row is not None:
                    rows.append(row)
            state = new_event(fields)
        consume_record(state, fields)

    if state is not None:
        row = finalize_event(state)
        if row is not None:
            rows.append(row)
    return rows


def read_support_counts(paths):
    """Build the AWK-style support dictionary in concatenation order."""

    support = {}
    for path in paths:
        with path.open() as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.strip()
                if not line:
                    continue
                fields = line.split()
                if len(fields) < 2:
                    raise ValueError(
                        "{}:{}: expected a key and a count".format(path, line_number)
                    )
                # Later records overwrite earlier ones, just as dic[$1]=$2
                # did after concatenating Sniffles then direct-call support.
                support[fields[0]] = awk_number(fields[1])
    return support


def precise_rows(frequency_rows, support):
    """Apply the second legacy AWK correction to insertion rows."""

    rows = []
    for fields in frequency_rows:
        read_support = support.get("{}:{}".format(fields[1], fields[2]), 0.0)
        if read_support <= 0:
            continue

        clipped = awk_number(fields[4])
        no_clipped_total = awk_number(fields[5])
        total = awk_number(fields[7])
        if read_support == no_clipped_total + 1:
            no_clipped_total += 1
            total += 1

        rows.append(
            [
                fields[0],
                fields[1],
                fields[2],
                awk_format(read_support),
                fields[4],
                awk_format(no_clipped_total),
                awk_format(read_support + clipped),
                awk_format(total),
                awk_format((read_support / no_clipped_total) * 100.0),
                awk_format(((read_support + clipped) / total) * 100.0),
                "INS",
            ]
        )
    return rows


def write_rows(path, rows):
    """Write tab-separated rows with the same final newline as AWK."""

    with path.open("w") as handle:
        for fields in rows:
            handle.write("\t".join(fields))
            handle.write("\n")


def concatenate_support(paths, destination):
    """Optionally materialize the exact Sniffles/direct byte concatenation."""

    with destination.open("wb") as output_handle:
        for path in paths:
            with path.open("rb") as input_handle:
                shutil.copyfileobj(input_handle, output_handle)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("read_counts", type=Path, help="sorted COUNT_READS.txt")
    parser.add_argument("sniffles_support", type=Path, help="Sniffles COUNT_TE_IN_RS.txt")
    parser.add_argument("direct_support", type=Path, help="direct-call COUNT_TE_IN_RS.txt")
    parser.add_argument("frequency_output", type=Path, help="FREQUENCY_TE_INS.tsv output")
    parser.add_argument(
        "precise_output", type=Path, help="FREQUENCY_TE_INS_PRECISE.tsv output"
    )
    parser.add_argument(
        "--combined-support-output",
        type=Path,
        help="optional byte concatenation of Sniffles then direct support",
    )
    parser.add_argument(
        "--pre-sorted",
        action="store_true",
        help="input is already ordered by the legacy event:read:type key",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    frequency = summarize_counts(args.read_counts, pre_sorted=args.pre_sorted)
    support_paths = (args.sniffles_support, args.direct_support)
    support = read_support_counts(support_paths)

    write_rows(args.frequency_output, frequency)
    write_rows(args.precise_output, precise_rows(frequency, support))
    if args.combined_support_output is not None:
        concatenate_support(support_paths, args.combined_support_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
