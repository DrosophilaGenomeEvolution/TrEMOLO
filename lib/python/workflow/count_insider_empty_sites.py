#!/usr/bin/env python3
"""Count empty-site reads for legacy-compatible INSIDER frequencies.

This is a deterministic replacement for
``lib/python/parsing/find_deletions.py``.  Its CIGAR-coordinate semantics are
intentionally historical: only ``M``, ``D``, and ``=`` consume reference
coordinates.  In particular, ``N`` and ``X`` are omitted so that existing
TrEMOLO results remain reproducible while the legacy frequency model is being
migrated.
"""

from __future__ import print_function

import argparse
import atexit
import multiprocessing
import os
from pathlib import Path
import tempfile

try:
    import pysam
except ImportError:  # Allow parser/unit-test imports outside the container.
    pysam = None


_BAM = None


def positive_int(value):
    """Argparse converter accepting strictly positive integers."""

    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("must be an integer")
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return parsed


def nonnegative_int(value):
    """Argparse converter accepting zero or a positive integer."""

    try:
        parsed = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("must be an integer")
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be greater than or equal to zero")
    return parsed


def build_parser():
    """Build the command-line parser without parsing arguments at import."""

    parser = argparse.ArgumentParser(
        description="Count empty-site BAM reads for INSIDER TE candidates",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("bam", type=Path, help="indexed input BAM")
    parser.add_argument(
        "candidates",
        type=Path,
        help="four-column INSIDER insertion BED, without a header",
    )
    parser.add_argument("output", type=Path, help="five-column count output")
    parser.add_argument(
        "--size-window",
        type=nonnegative_int,
        default=30,
        help="allowed difference between deletion and candidate lengths",
    )
    parser.add_argument(
        "--breakpoint-distance",
        type=positive_int,
        default=30,
        help="strict maximum distance from a deletion endpoint to a breakpoint",
    )
    parser.add_argument(
        "--threads",
        type=positive_int,
        default=1,
        help="number of independent BAM-reading worker processes",
    )
    return parser


def read_candidates(path):
    """Read candidate BED rows in input order.

    Returned tuples contain all fields required by a worker plus the source
    line number used in validation messages.
    """

    candidates = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    "{}:{}: expected at least four tab-delimited columns".format(
                        path, line_number
                    )
                )
            chromosome = fields[0]
            if not chromosome:
                raise ValueError("{}:{}: empty chromosome".format(path, line_number))
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise ValueError(
                    "{}:{}: start and end must be integers".format(path, line_number)
                )
            if start < 0:
                raise ValueError("{}:{}: start must not be negative".format(path, line_number))
            if end < start:
                raise ValueError(
                    "{}:{}: end must be greater than or equal to start".format(
                        path, line_number
                    )
                )
            name = fields[3].strip()
            if not name:
                raise ValueError("{}:{}: empty candidate name".format(path, line_number))
            candidates.append((chromosome, start, end, name, line_number))
    return candidates


def close_worker_bam():
    """Close the process-local BAM handle, if one was opened."""

    global _BAM
    if _BAM is not None:
        _BAM.close()
        _BAM = None


def initialize_worker(bam_path):
    """Open one independent indexed BAM handle in a worker process."""

    global _BAM
    if pysam is None:
        raise RuntimeError("pysam is required to count INSIDER empty sites")
    _BAM = pysam.AlignmentFile(str(bam_path), "rb")
    atexit.register(close_worker_bam)


def count_candidate_with_bam(candidate, bam, size_window, breakpoint_distance):
    """Count legacy-compatible empty-site reads for one candidate."""

    chromosome, start, end, name, _line_number = candidate
    candidate_size = end - start
    minimum_size = candidate_size - size_window
    maximum_size = candidate_size + size_window
    supporting_reads = set()

    for read in bam.fetch(str(chromosome), int(start), int(end)):
        read_name = read.query_name
        if not read_name or read_name in supporting_reads:
            continue

        reference_offset = 0
        for operation, length in read.cigartuples:
            # Compatibility contract: the legacy implementation omitted N
            # and X even though both consume reference coordinates.
            if operation in (0, 2, 7):  # M, D, =
                reference_offset += length

            deletion_endpoint = reference_offset + read.reference_start
            if (
                operation == 2
                and minimum_size <= length <= maximum_size
                and (
                    abs(deletion_endpoint - start) < breakpoint_distance
                    or abs(deletion_endpoint - end) < breakpoint_distance
                )
            ):
                supporting_reads.add(read_name)
                break

    return "\t".join(
        (chromosome, str(start), str(end), name, str(len(supporting_reads)))
    ) + "\n"


def count_candidate(task):
    """Multiprocessing adapter using the process-local BAM handle."""

    candidate, size_window, breakpoint_distance = task
    if _BAM is None:
        raise RuntimeError("INSIDER BAM worker was not initialized")
    return count_candidate_with_bam(
        candidate, _BAM, size_window, breakpoint_distance
    )


def count_all(bam_path, candidates, size_window, breakpoint_distance, threads):
    """Return formatted result rows in deterministic candidate order."""

    if pysam is None:
        raise RuntimeError("pysam is required to count INSIDER empty sites")
    if not candidates:
        return []

    if threads == 1:
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            return [
                count_candidate_with_bam(
                    candidate, bam, size_window, breakpoint_distance
                )
                for candidate in candidates
            ]

    tasks = [
        (candidate, size_window, breakpoint_distance) for candidate in candidates
    ]
    worker_count = min(threads, len(tasks))
    chunk_size = max(1, len(tasks) // (worker_count * 4))
    with multiprocessing.Pool(
        processes=worker_count,
        initializer=initialize_worker,
        initargs=(str(bam_path),),
    ) as pool:
        return pool.map(count_candidate, tasks, chunksize=chunk_size)


def atomic_write(path, rows):
    """Replace ``path`` only after every row has been written successfully."""

    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=path.name + ".", suffix=".tmp", dir=str(path.parent)
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.writelines(rows)
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
    candidates = read_candidates(args.candidates)
    rows = count_all(
        args.bam,
        candidates,
        args.size_window,
        args.breakpoint_distance,
        args.threads,
    )
    atomic_write(args.output, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
