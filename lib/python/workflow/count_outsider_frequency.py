#!/usr/bin/env python3
"""Count read evidence around OUTSIDER TE candidates.

This is a dependency-light, deterministic replacement for
``lib/python/parsing/getFrequency.py``.  The 19-column output deliberately
retains the historical coordinate and CIGAR semantics because the downstream
frequency aggregation currently relies on them.
"""

from __future__ import print_function

import argparse
import csv
import multiprocessing
from pathlib import Path

try:
    import pysam
except ImportError:  # Allow parser and unit tests to import this module.
    pysam = None


REQUIRED_CANDIDATE_COLUMNS = ("sseqid", "qseqid")

_BAM = None
_TE_SIZES = None
_SV_SIZES = None
_WINDOW = None
_MIN_SIZE = None
_SIZE_PERCENT = None


def positive_int(value):
    """Argparse converter accepting strictly positive integers."""

    try:
        parsed = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def build_parser():
    """Build the command-line parser without parsing at import time."""

    parser = argparse.ArgumentParser(
        description="Count BAM read evidence around OUTSIDER TE candidates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", type=Path, help="indexed input BAM")
    parser.add_argument(
        "candidates", type=Path, help="tab-delimited candidate CSV with a header"
    )
    parser.add_argument("te_size", type=Path, help="tab-delimited TE_SIZE.tsv")
    parser.add_argument("sv_size", type=Path, help="tab-delimited SV_SIZE.tsv")
    parser.add_argument("output", type=Path, help="19-column read-evidence output")
    parser.add_argument(
        "-p",
        "--per-size",
        dest="per_size",
        type=int,
        default=80,
        help="minimum insertion length as a percentage of the TE length",
    )
    parser.add_argument(
        "-w",
        "--window",
        type=int,
        default=50,
        help="maximum distance around each candidate breakpoint",
    )
    parser.add_argument(
        "-s",
        "--min-size",
        dest="min_size",
        type=int,
        default=200,
        help="minimum soft/hard-clipped sequence length",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=positive_int,
        default=1,
        help="number of BAM-reading worker processes",
    )
    return parser


def read_candidates(path):
    """Read the legacy tab-delimited candidate table in file order."""

    candidates = []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("{}: missing candidate header".format(path))
        missing = [
            column
            for column in REQUIRED_CANDIDATE_COLUMNS
            if column not in reader.fieldnames
        ]
        if missing:
            raise ValueError(
                "{}: missing candidate column(s): {}".format(path, ", ".join(missing))
            )
        for line_number, row in enumerate(reader, 2):
            if not row or all(value in (None, "") for value in row.values()):
                continue
            sseqid = row.get("sseqid")
            qseqid = row.get("qseqid")
            if not sseqid or not qseqid:
                raise ValueError(
                    "{}:{}: empty sseqid or qseqid".format(path, line_number)
                )
            candidates.append((sseqid, qseqid, line_number))
    return candidates


def read_te_sizes(path):
    """Read ``TE_SIZE.tsv`` while retaining its historical string values."""

    sizes = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, 1):
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(
                    "{}:{}: expected two tab-delimited columns".format(
                        path, line_number
                    )
                )
            sizes[row[0]] = row[1]
    return sizes


def read_sv_sizes(path):
    """Read ``SV_SIZE.tsv``, keyed by the event ID embedded in column one."""

    sizes = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, 1):
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(
                    "{}:{}: expected two tab-delimited columns".format(
                        path, line_number
                    )
                )
            fields = row[0].split(":")
            if len(fields) <= 4:
                raise ValueError(
                    "{}:{}: sequence name has no event ID: {}".format(
                        path, line_number, row[0]
                    )
                )
            sizes[fields[4]] = row[1]
    return sizes


def parse_candidate(candidate):
    """Extract the fields used by the legacy counter from one candidate."""

    sseqid, qseqid, line_number = candidate
    fields = qseqid.split(":")
    if len(fields) <= 5:
        raise ValueError(
            "candidate line {} has a malformed qseqid: {}".format(line_number, qseqid)
        )
    event_id = fields[4]
    event_fields = event_id.split(".")
    if len(event_fields) <= 1:
        raise ValueError(
            "candidate line {} has a malformed event ID: {}".format(
                line_number, event_id
            )
        )
    try:
        start = int(fields[2])
        end = int(fields[3])
    except ValueError:
        raise ValueError(
            "candidate line {} has non-integer coordinates: {}:{}".format(
                line_number, fields[2], fields[3]
            )
        )
    return {
        "sseqid": sseqid,
        "qseqid": qseqid,
        "event_id": event_id,
        "read_support": fields[5],
        "sv_type": event_fields[1],
        "chromosome": fields[0],
        "start": start,
        "end": end,
    }


def initialize_worker(bam_path, te_sizes, sv_sizes, window, min_size, per_size):
    """Open one independent BAM handle in each worker process."""

    global _BAM, _TE_SIZES, _SV_SIZES, _WINDOW, _MIN_SIZE, _SIZE_PERCENT
    if pysam is None:
        raise RuntimeError("pysam is required to count BAM evidence")
    _BAM = pysam.AlignmentFile(str(bam_path), "rb")
    _TE_SIZES = te_sizes
    _SV_SIZES = sv_sizes
    _WINDOW = window
    _MIN_SIZE = min_size
    _SIZE_PERCENT = per_size


def format_record(record_type, infos, trailing_size, sequence):
    """Format and guard the historical 19-space-delimited-column layout."""

    fields = [record_type] + [str(value) for value in infos] + [
        str(trailing_size),
        str(sequence),
    ]
    if len(fields) != 19:
        raise AssertionError("internal error: expected 19 output columns")
    return " ".join(fields) + "\n"


def count_candidate(candidate):
    """Return every legacy evidence row for one candidate, in BAM order."""

    parsed = parse_candidate(candidate)
    sseqid = parsed["sseqid"]
    qseqid = parsed["qseqid"]
    event_id = parsed["event_id"]
    read_support = parsed["read_support"]
    sv_type = parsed["sv_type"]
    chromosome = parsed["chromosome"]
    start = parsed["start"]
    end = parsed["end"]

    # Let missing keys fail loudly, just as the legacy implementation did.
    sv_size = _SV_SIZES[event_id]
    te_size = _TE_SIZES[sseqid]
    locus = "{}:{}-{}".format(chromosome, start, end)
    records = []
    hard_clip_identities = set()

    reads = _BAM.fetch(
        str(chromosome), max(start - _WINDOW, 1), end + _WINDOW
    )
    for read in reads:
        reference_start = read.reference_start
        sequence = read.seq
        read_name = read.query_name
        reference_name = read.reference_name

        count_ref = 0
        count_read = 0
        count_read_real = 0

        if read_name:
            found = False
            for operation, length in read.cigartuples:
                # These operation sets intentionally omit N and X.  Changing
                # them would change the historical breakpoint coordinates.
                if operation in (0, 2, 7):
                    count_ref += length
                if operation in (0, 1, 7, 4):
                    count_read += length
                if operation in (0, 1, 7, 4, 5):
                    count_read_real += length

                breakpoint = reference_start + count_ref
                distance = min(abs(breakpoint - end), abs(breakpoint - start))
                infos = [
                    read_name,
                    event_id,
                    sseqid,
                    read_support,
                    sv_size,
                    te_size,
                    ".",
                    breakpoint,
                    locus,
                    qseqid,
                    distance,
                    reference_start,
                    count_ref,
                    count_read,
                    count_read_real,
                    sv_type,
                ]

                if sv_type == "INS":
                    if max(start - _WINDOW, 1) <= breakpoint <= end + _WINDOW:
                        if operation == 1 and length >= (
                            (_SIZE_PERCENT / 100) * int(te_size)
                        ):
                            infos[6] = length
                            if sequence:
                                evidence_sequence = sequence[
                                    count_read - length : count_read
                                ]
                            else:
                                evidence_sequence = "."
                            records.append(
                                format_record(
                                    "I", infos, length, evidence_sequence
                                )
                            )
                            found = True

                        # The outer >= and inner > checks are both retained.
                        if operation == 5 and length >= _MIN_SIZE:
                            if sequence and length > _MIN_SIZE:
                                identity = ":".join(
                                    [
                                        str(reference_name),
                                        str(breakpoint),
                                        str(read_name),
                                        str(count_read_real),
                                        "L",
                                        str(len(sequence)),
                                    ]
                                )
                            else:
                                identity = None
                            if (
                                identity is not None
                                and identity not in hard_clip_identities
                            ):
                                infos[6] = length
                                hard_clip_identities.add(identity)
                                records.append(format_record("H", infos, ".", "."))
                                found = True

                        if operation == 4 and length >= _MIN_SIZE:
                            if sequence and length > _MIN_SIZE:
                                infos[6] = length
                                evidence_sequence = sequence[
                                    count_read - length : count_read
                                ]
                                records.append(
                                    format_record(
                                        "S", infos, length, evidence_sequence
                                    )
                                )
                                found = True
                else:
                    deletion_start = breakpoint - length
                    if (
                        max(start - _WINDOW, 1)
                        <= deletion_start
                        <= end + _WINDOW
                    ):
                        if operation == 2 and length >= int(sv_size):
                            infos[7] = deletion_start
                            infos[6] = length
                            if sequence:
                                if (
                                    0 < count_read - 30
                                    and count_read + 30 < len(sequence)
                                ):
                                    evidence_sequence = (
                                        sequence[count_read - 30 : count_read]
                                        + ":"
                                        + sequence[count_read : count_read + 30]
                                    )
                                else:
                                    evidence_sequence = "."
                            else:
                                evidence_sequence = "."
                            records.append(
                                format_record("D", infos, 20, evidence_sequence)
                            )
                            found = True

            distance = min(
                abs((reference_start + count_ref) - end),
                abs((reference_start + count_ref) - start),
            )
            if not found:
                infos = [
                    read_name,
                    event_id,
                    sseqid,
                    read_support,
                    sv_size,
                    te_size,
                    ".",
                    reference_start + count_ref,
                    locus,
                    qseqid,
                    distance,
                    reference_start,
                    count_ref,
                    count_read,
                    count_read_real,
                    sv_type,
                ]
                records.append(format_record("E", infos, ".", "."))

    return "".join(records)


def run(args):
    """Count all candidates in parallel and write them in input order."""

    if pysam is None:
        raise RuntimeError("pysam is required to count BAM evidence")

    candidates = read_candidates(args.candidates)
    te_sizes = read_te_sizes(args.te_size)
    sv_sizes = read_sv_sizes(args.sv_size)
    initializer_args = (
        args.bam,
        te_sizes,
        sv_sizes,
        args.window,
        args.min_size,
        args.per_size,
    )

    with args.output.open("w", newline="") as handle:
        # imap preserves candidate order and propagates worker failures while
        # streaming completed candidates instead of retaining the whole
        # cohort's read evidence in memory.
        with multiprocessing.Pool(
            processes=args.threads,
            initializer=initialize_worker,
            initargs=initializer_args,
        ) as pool:
            for result in pool.imap(count_candidate, candidates, chunksize=1):
                handle.write(result)


def main(argv=None):
    """CLI entry point."""

    args = build_parser().parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
