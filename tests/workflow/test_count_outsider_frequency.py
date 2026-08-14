#!/usr/bin/env python3

import argparse
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
import count_outsider_frequency as frequency_counter  # noqa: E402


class FakeRead(object):
    def __init__(
        self,
        name,
        reference_start,
        cigartuples,
        sequence="A" * 80,
        reference_name="chr1",
    ):
        self.query_name = name
        self.reference_start = reference_start
        self.cigartuples = cigartuples
        self.seq = sequence
        self.reference_name = reference_name


class FakeBam(object):
    def __init__(self, reads):
        self.reads = list(reads)
        self.fetch_calls = []

    def fetch(self, chromosome, start, end):
        self.fetch_calls.append((chromosome, start, end))
        return iter(self.reads)


def candidate(event_id, support="1", start=100, end=100, te="roo"):
    return (
        te,
        "chr1:sample:{}:{}:{}:{}".format(start, end, event_id, support),
        2,
    )


def configure_counter(bam, event_id, sv_size="5", te_size="10"):
    frequency_counter._BAM = bam
    frequency_counter._TE_SIZES = {"roo": te_size}
    frequency_counter._SV_SIZES = {event_id: sv_size}
    frequency_counter._WINDOW = 50
    frequency_counter._MIN_SIZE = 3
    frequency_counter._SIZE_PERCENT = 80


class FrequencyInputParsingTests(unittest.TestCase):
    def test_reads_and_parses_candidate_and_size_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "candidates.tsv"
            candidates.write_text(
                "sseqid\tqseqid\textra\n"
                "roo\tchr2:sample:101:103:sniffles.INS.7:4\tignored\n"
                "\t\t\n"
            )
            te_sizes = root / "TE_SIZE.tsv"
            te_sizes.write_text("roo\t123\ncopia\t456\nroo\t124\n")
            sv_sizes = root / "SV_SIZE.tsv"
            sv_sizes.write_text(
                "chr2:sample:101:103:sniffles.INS.7:first\t80\n"
                "chr2:sample:101:103:sniffles.INS.7:last\t81\n"
            )

            parsed_candidates = frequency_counter.read_candidates(candidates)

            self.assertEqual(
                parsed_candidates,
                [("roo", "chr2:sample:101:103:sniffles.INS.7:4", 2)],
            )
            self.assertEqual(
                frequency_counter.parse_candidate(parsed_candidates[0]),
                {
                    "sseqid": "roo",
                    "qseqid": "chr2:sample:101:103:sniffles.INS.7:4",
                    "event_id": "sniffles.INS.7",
                    "read_support": "4",
                    "sv_type": "INS",
                    "chromosome": "chr2",
                    "start": 101,
                    "end": 103,
                },
            )
            self.assertEqual(
                frequency_counter.read_te_sizes(te_sizes),
                {"roo": "124", "copia": "456"},
            )
            self.assertEqual(
                frequency_counter.read_sv_sizes(sv_sizes),
                {"sniffles.INS.7": "81"},
            )

    def test_candidate_table_requires_the_named_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "candidates.tsv"
            path.write_text("family\tidentifier\nroo\tevent\n")

            with self.assertRaisesRegex(ValueError, "sseqid, qseqid"):
                frequency_counter.read_candidates(path)

    def test_cli_accepts_only_positive_worker_counts(self):
        parser = frequency_counter.build_parser()
        arguments = ["input.bam", "calls.tsv", "te.tsv", "sv.tsv", "out.txt"]

        self.assertEqual(parser.parse_args(arguments + ["-t", "3"]).threads, 3)
        self.assertEqual(frequency_counter.positive_int("2"), 2)
        with self.assertRaises(argparse.ArgumentTypeError):
            frequency_counter.positive_int("not-an-integer")
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                parser.parse_args(arguments + ["--threads", "0"])


class CountCandidateTests(unittest.TestCase):
    def test_insertion_and_empty_evidence_are_deterministic(self):
        event_id = "sniffles.INS.1"
        bam = FakeBam(
            [
                FakeRead("insertion-read", 95, [(0, 5), (1, 8), (0, 5)]),
                FakeRead("reference-read", 90, [(0, 20)]),
            ]
        )
        configure_counter(bam, event_id)
        value = candidate(event_id, support="2")

        first = frequency_counter.count_candidate(value)
        second = frequency_counter.count_candidate(value)

        self.assertEqual(first, second)
        rows = [line.split() for line in first.splitlines()]
        self.assertEqual([row[0] for row in rows], ["I", "E"])
        self.assertEqual([row[1] for row in rows], ["insertion-read", "reference-read"])
        self.assertEqual(rows[0][17], "8")
        self.assertEqual(rows[0][18], "A" * 8)
        self.assertTrue(all(len(row) == 19 for row in rows))
        self.assertEqual(bam.fetch_calls, [("chr1", 50, 150), ("chr1", 50, 150)])

    def test_deletion_evidence_uses_the_deletion_breakpoint(self):
        event_id = "sniffles.DEL.2"
        bam = FakeBam(
            [
                FakeRead("deletion-read", 95, [(0, 5), (2, 5), (0, 5)]),
                FakeRead("reference-read", 90, [(0, 20)]),
            ]
        )
        configure_counter(bam, event_id, sv_size="5")

        rows = [
            line.split()
            for line in frequency_counter.count_candidate(candidate(event_id)).splitlines()
        ]

        self.assertEqual([row[0] for row in rows], ["D", "E"])
        self.assertEqual(rows[0][8], "100")
        self.assertEqual(rows[0][16], "DEL")
        self.assertEqual(rows[0][17:], ["20", "."])
        self.assertTrue(all(len(row) == 19 for row in rows))

    def test_hard_clip_without_bam_sequence_becomes_empty_evidence(self):
        event_id = "sniffles.INS.3"
        bam = FakeBam(
            [FakeRead("sequence-omitted", 95, [(5, 10), (0, 10)], sequence=None)]
        )
        configure_counter(bam, event_id)

        rows = frequency_counter.count_candidate(candidate(event_id)).splitlines()

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0].split()[0], "E")
        self.assertEqual(len(rows[0].split()), 19)


if __name__ == "__main__":
    unittest.main()
