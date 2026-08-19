#!/usr/bin/env python3

import argparse
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
import count_insider_empty_sites as empty_site_counter  # noqa: E402


class FakeRead(object):
    def __init__(self, name, reference_start, cigartuples):
        self.query_name = name
        self.reference_start = reference_start
        self.cigartuples = cigartuples


class FakeBam(object):
    def __init__(self, reads=None, reads_by_region=None):
        self.reads = list(reads or [])
        self.reads_by_region = dict(reads_by_region or {})
        self.fetch_calls = []

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception, traceback):
        return False

    def fetch(self, chromosome, start, end):
        region = (chromosome, start, end)
        self.fetch_calls.append(region)
        return iter(self.reads_by_region.get(region, self.reads))


class FakePysam(object):
    def __init__(self, bam):
        self.bam = bam
        self.open_calls = []

    def AlignmentFile(self, path, mode):
        self.open_calls.append((path, mode))
        return self.bam


def candidate(start=100, end=200, name="roo|Assemblytics_w_1", chromosome="chr1"):
    return (chromosome, start, end, name, 1)


def count(reads, value=None, size_window=30, breakpoint_distance=30):
    return empty_site_counter.count_candidate_with_bam(
        value or candidate(),
        FakeBam(reads),
        size_window,
        breakpoint_distance,
    )


class EmptySiteCountingTests(unittest.TestCase):
    def test_deletion_length_window_is_inclusive_at_both_boundaries(self):
        rows = count(
            [
                FakeRead("lower-bound", 130, [(2, 70)]),
                FakeRead("upper-bound", 70, [(2, 130)]),
                FakeRead("below-window", 131, [(2, 69)]),
                FakeRead("above-window", 69, [(2, 131)]),
            ],
            size_window=30,
            breakpoint_distance=1,
        )

        self.assertEqual(
            rows,
            "chr1\t100\t200\troo|Assemblytics_w_1\t2\n",
        )

    def test_breakpoint_distance_is_strict(self):
        rows = count(
            [
                FakeRead("inside-distance", 71, [(2, 100)]),
                FakeRead("at-distance", 70, [(2, 100)]),
            ],
            size_window=0,
            breakpoint_distance=30,
        )

        self.assertTrue(rows.endswith("\t1\n"))

    def test_deletion_endpoint_is_evaluated_after_consuming_the_deletion(self):
        # Before consuming this 70-base D, the reference position is 130 and
        # is not within one base of either breakpoint.  Its endpoint is 200.
        rows = count(
            [FakeRead("endpoint-after-D", 130, [(2, 70)])],
            size_window=30,
            breakpoint_distance=1,
        )

        self.assertTrue(rows.endswith("\t1\n"))

    def test_legacy_reference_offset_includes_m_d_equal_but_omits_x_and_n(self):
        value = candidate(start=100, end=170)
        rows = count(
            [
                FakeRead("M-consumes", 0, [(0, 100), (2, 70)]),
                FakeRead("D-consumes", 99, [(2, 1), (2, 70)]),
                FakeRead("equal-consumes", 0, [(7, 100), (2, 70)]),
                FakeRead("X-is-omitted", 100, [(8, 10), (2, 70)]),
                FakeRead("N-is-omitted", 100, [(3, 10), (2, 70)]),
            ],
            value=value,
            size_window=0,
            breakpoint_distance=1,
        )

        self.assertEqual(rows, "chr1\t100\t170\troo|Assemblytics_w_1\t5\n")

    def test_query_names_are_counted_at_most_once_per_candidate(self):
        rows = count(
            [
                FakeRead("same-read", 100, [(2, 100)]),
                FakeRead("same-read", 100, [(2, 100)]),
                FakeRead("", 100, [(2, 100)]),
            ],
            size_window=0,
            breakpoint_distance=1,
        )

        self.assertTrue(rows.endswith("\t1\n"))

    def test_serial_counting_preserves_candidate_order_deterministically(self):
        candidates = [
            candidate(100, 200, "roo|event-2", "chr2"),
            candidate(300, 400, "copia|event-1", "chr1"),
        ]
        bam = FakeBam()
        fake_pysam = FakePysam(bam)

        with patch.object(empty_site_counter, "pysam", fake_pysam):
            first = empty_site_counter.count_all(
                Path("unused.bam"), candidates, 30, 30, 1
            )
            second = empty_site_counter.count_all(
                Path("unused.bam"), candidates, 30, 30, 1
            )

        expected = [
            "chr2\t100\t200\troo|event-2\t0\n",
            "chr1\t300\t400\tcopia|event-1\t0\n",
        ]
        self.assertEqual(first, expected)
        self.assertEqual(second, expected)
        self.assertEqual(
            bam.fetch_calls,
            [
                ("chr2", 100, 200),
                ("chr1", 300, 400),
                ("chr2", 100, 200),
                ("chr1", 300, 400),
            ],
        )


class EmptySiteInputAndCliTests(unittest.TestCase):
    def test_cli_accepts_only_positive_worker_counts(self):
        parser = empty_site_counter.build_parser()
        positional = ["input.bam", "insertions.bed", "counts.bed"]

        self.assertEqual(
            parser.parse_args(positional + ["--threads", "3"]).threads,
            3,
        )
        self.assertEqual(empty_site_counter.positive_int("2"), 2)
        with self.assertRaises(argparse.ArgumentTypeError):
            empty_site_counter.positive_int("not-an-integer")
        for value in ("0", "-1"):
            with self.subTest(value=value):
                with redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        parser.parse_args(positional + ["--threads", value])

    def test_zero_candidates_replace_stale_output_with_an_empty_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "INSERTION_TE.bed"
            candidates.write_bytes(b"")
            output = root / "DEL_NB.bed"
            output.write_text("stale result\n")

            with patch.object(empty_site_counter, "pysam", object()):
                result = empty_site_counter.main(
                    [
                        str(root / "unused.bam"),
                        str(candidates),
                        str(output),
                        "--threads",
                        "8",
                    ]
                )

            self.assertEqual(result, 0)
            self.assertEqual(output.read_bytes(), b"")


if __name__ == "__main__":
    unittest.main()
