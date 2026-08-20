#!/usr/bin/env python3

import argparse
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/TSD"))
from getTSDgenome import (  # noqa: E402
    call_tsds,
    find_tsd,
    parse_candidate,
    positive_integer,
)


ORACLE_RECORD = (
    "3406:X_RaGOO_RaGOO:8083209-8083219:1360|Assemblytics_w_888\t"
    "1360\tAssemblytics_w_888\tGAGTCAAACC\tGTCAAACTAT\t.\t."
)


class CandidateParsingTests(unittest.TestCase):
    def test_parses_insider_flank_record(self):
        parsed, error = parse_candidate(ORACLE_RECORD, line_number=2)

        self.assertIsNone(error)
        self.assertEqual(parsed.te, "1360")
        self.assertEqual(parsed.identifier, "Assemblytics_w_888")
        self.assertEqual(parsed.te_size, 3406)
        self.assertEqual(parsed.genome_position, 8083209)

    def test_rejects_wrong_width_and_ambiguous_bases(self):
        parsed, error = parse_candidate("one\ttwo", line_number=3)
        self.assertIsNone(parsed)
        self.assertIn("expected 7", error)

        parsed, error = parse_candidate(ORACLE_RECORD.replace("GAGT", "GART"))
        self.assertIsNone(parsed)
        self.assertIn("unsupported ambiguous base(s): R", error)


class ThreadValidationTests(unittest.TestCase):
    def test_positive_threads_only(self):
        self.assertEqual(positive_integer("8"), 8)
        for value in ("0", "-1", "bad"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    positive_integer(value)


class TsdCallingTests(unittest.TestCase):
    def test_reproduces_first_insider_oracle_call(self):
        parsed, error = parse_candidate(ORACLE_RECORD)
        self.assertIsNone(error)
        self.assertEqual(
            find_tsd(parsed),
            [
                "1360|Assemblytics_w_888",
                "GTCAAAC",
                "8083208",
                "3407",
                "3406:X_RaGOO_RaGOO:8083209-8083219:1360|Assemblytics_w_888",
            ],
        )

    def test_zero_candidates_do_not_start_workers(self):
        self.assertEqual(call_tsds([], 4), [])


if __name__ == "__main__":
    unittest.main()
