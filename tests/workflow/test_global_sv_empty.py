#!/usr/bin/env python3

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "lib/python/parsing/global_sv.py"

RAW_HEADER = (
    "sseqid\tqseqid\tpident\tsize_per\tsize_el\tmismatch\tgapopen\t"
    "qstart\tqend\tsstart\tsend\tevalue\tbitscore\n"
)
COMBINE_HEADER = (
    "sseqid\tqseqid\tgrain_pident\tsize_per\tsize_el\tqstart\tqend\t"
    "sstart\tsend\n"
)


def run_global_sv(blast, database, raw, combined, extra_options=(), no_site=False):
    command = [sys.executable]
    if no_site:
        command.append("-S")
    command.extend(
        [
            str(SCRIPT),
            *extra_options,
            "--combine_name",
            str(combined),
            str(blast),
            str(database),
            str(raw),
        ]
    )
    return subprocess.run(command, capture_output=True, text=True, check=False)


class GlobalSvEmptyInputTests(unittest.TestCase):
    def test_empty_blast_replaces_stale_outputs_with_exact_headers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast = root / "empty.bln"
            database = root / "te.fasta"
            raw = root / "calls.tsv"
            combined = root / "combined.tsv"

            blast.write_bytes(b"")
            database.write_text(">roo\nACGT\n")
            raw.write_text("stale raw result\n")
            combined.write_text("stale combined result\n")

            result = run_global_sv(
                blast,
                database,
                raw,
                combined,
                extra_options=(
                    "--min-pident",
                    "100",
                    "--min-size-percent",
                    "100",
                    "--combine",
                ),
                no_site=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(raw.read_text(), RAW_HEADER)
            self.assertEqual(combined.read_text(), COMBINE_HEADER)

    def test_non_empty_malformed_blast_does_not_succeed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast = root / "malformed.bln"
            database = root / "te.fasta"
            raw = root / "calls.tsv"
            combined = root / "combined.tsv"

            blast.write_text("query\troo\tnot-enough-columns\n")
            database.write_text(">roo\nACGT\n")

            result = run_global_sv(blast, database, raw, combined)

            self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()
