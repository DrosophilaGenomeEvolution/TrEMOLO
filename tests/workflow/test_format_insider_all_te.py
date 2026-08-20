#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "lib/python/workflow/format_insider_all_te.py"
sys.path.insert(0, str(SCRIPT.parent))
from format_insider_all_te import EXPECTED_HEADER, read_calls  # noqa: E402


def row(family, query):
    return [
        family,
        query,
        "99.0",
        "95.0",
        "100",
        "0",
        "0",
        "1",
        "101",
        "1",
        "101",
        "1e-20",
        "200",
    ]


def write_calls(path, records, header=EXPECTED_HEADER):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(records)


class FormatInsiderAllTeTests(unittest.TestCase):
    def test_formats_legacy_identifiers_in_input_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = root / "ALL_TE.csv"
            output = root / "POSITION_ALL_TE.bed"
            write_calls(
                calls,
                [
                    row("roo", "13::chr1:100-250:+"),
                    row("copia", "7::assembly:hap1:300-500:-"),
                ],
            )

            subprocess.run(
                [sys.executable, str(SCRIPT), str(calls), str(output)],
                check=True,
            )

            self.assertEqual(
                output.read_text(),
                "chr1\t100\t250\troo|13\n"
                "assembly:hap1\t300\t500\tcopia|7\n",
            )

    def test_header_only_input_replaces_stale_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = root / "ALL_TE.csv"
            output = root / "POSITION_ALL_TE.bed"
            write_calls(calls, [])
            output.write_text("STALE\n")

            subprocess.run(
                [sys.executable, str(SCRIPT), str(calls), str(output)],
                check=True,
            )

            self.assertEqual(output.read_bytes(), b"")

    def test_rejects_an_unexpected_schema(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "ALL_TE.csv"
            write_calls(path, [], header=("wrong", "schema"))
            with self.assertRaisesRegex(ValueError, "unexpected.*schema"):
                read_calls(path)

    def test_rejects_malformed_identifiers_and_families(self):
        invalid = (
            row("roo", "missing-structure"),
            row("roo", "1::chr1:300-200:+"),
            row("bad|family", "1::chr1:100-200:+"),
        )
        for record in invalid:
            with self.subTest(record=record[:2]), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / "ALL_TE.csv"
                write_calls(path, [record])
                with self.assertRaises(ValueError):
                    read_calls(path)


if __name__ == "__main__":
    unittest.main()
