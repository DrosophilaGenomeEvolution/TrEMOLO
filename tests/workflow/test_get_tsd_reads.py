#!/usr/bin/env python3

import argparse
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/TSD"))
from getTSDreads import find_tsd, parse_candidate, positive_integer  # noqa: E402


ORACLE_RECORD = (
    "8114:X_RaGOO_RaGOO:<INS>:22714472:22714473:"
    "TrEMOLO.INS.84:36:IMPRECISE:30\t"
    "Max-element\tTrEMOLO.INS.84\tTGATATTGTCC\tGTCCTATCAT\t"
    "CCAAACACCT\tCGTCATCAAT"
)


class CandidateParsingTests(unittest.TestCase):
    def test_parses_seven_column_record(self):
        parsed, error = parse_candidate(ORACLE_RECORD, line_number=7)

        self.assertIsNone(error)
        self.assertEqual(parsed.te, "Max-element")
        self.assertEqual(parsed.identifier, "TrEMOLO.INS.84")
        self.assertEqual(parsed.te_size, 8114)
        self.assertEqual(parsed.genome_position, 22714472)
        self.assertEqual(parsed.flank_left, "TGATATTGTCC")
        self.assertEqual(parsed.genome_flank_right, "CGTCATCAAT")

    def test_rejects_record_with_wrong_column_count(self):
        parsed, error = parse_candidate("one\ttwo\tthree", line_number=3)

        self.assertIsNone(parsed)
        self.assertEqual(
            error, "line 3: expected 7 tab-separated columns, got 3"
        )

    def test_rejects_unsupported_ambiguous_base(self):
        record = ORACLE_RECORD.replace("TGATATTGTCC", "TGATARTGTCC")

        parsed, error = parse_candidate(record, line_number=11)

        self.assertIsNone(parsed)
        self.assertEqual(error, "line 11: unsupported ambiguous base(s): R")


class ThreadValidationTests(unittest.TestCase):
    def test_accepts_positive_integer(self):
        self.assertEqual(positive_integer("1"), 1)
        self.assertEqual(positive_integer("12"), 12)

    def test_rejects_zero_negative_and_non_integer_values(self):
        for value in ("0", "-1", "not-an-integer"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    positive_integer(value)


class TsdCallingTests(unittest.TestCase):
    def test_reproduces_max_element_oracle_call(self):
        parsed, error = parse_candidate(ORACLE_RECORD)
        self.assertIsNone(error)

        self.assertEqual(
            find_tsd(parsed),
            [
                "Max-element|TrEMOLO.INS.84",
                "GTCC",
                "22714472",
                "8114",
                (
                    "8114:X_RaGOO_RaGOO:<INS>:22714472:22714473:"
                    "TrEMOLO.INS.84:36:IMPRECISE:30"
                ),
            ],
        )


if __name__ == "__main__":
    unittest.main()
