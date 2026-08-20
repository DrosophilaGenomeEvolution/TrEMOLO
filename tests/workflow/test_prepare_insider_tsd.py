#!/usr/bin/env python3

import argparse
import sys
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from prepare_insider_tsd import prepare, positive_integer  # noqa: E402


class FakeIndexedFasta(object):
    def __init__(self, records):
        self.records = OrderedDict(records)
        self.references = tuple(self.records)
        self.lengths = tuple(len(sequence) for sequence in self.records.values())
        self.closed = False

    def fetch(self, name, start, end):
        return self.records[name][start:end]

    def close(self):
        self.closed = True


class PrepareInsiderTsdTests(unittest.TestCase):
    def test_reproduces_legacy_tables_and_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fasta"
            genome.write_text(">chr1\n" + "ACGT" * 20 + "\n")
            index = root / "genome.fasta.fai"
            index.write_text("chr1\t80\t6\t80\t81\n")
            insertion = root / "insertions.bed"
            insertion.write_text("chr1\t20\t25\troo|ins1\t+\n")
            deletion = root / "deletions.bed"
            deletion.write_text("chr1\t50\t50\tcopia|del1\t-\n")
            args = SimpleNamespace(
                genome=genome,
                genome_index=index,
                insertion_bed=insertion,
                deletion_bed=deletion,
                flank_size=4,
                merged_bed=root / "merged.bed",
                flank_bed=root / "flanks.bed",
                formatted_flanks=root / "formatted.bed",
            )
            indexed = FakeIndexedFasta({"chr1": "ACGT" * 20})
            with patch(
                "prepare_insider_tsd.open_indexed_genome", return_value=indexed
            ):
                self.assertEqual(prepare(args), 2)

            self.assertTrue(indexed.closed)
            self.assertEqual(
                args.merged_bed.read_text(),
                "chr1\t20\t25\troo|ins1\nchr1\t50\t50\tcopia|del1\n",
            )
            self.assertEqual(
                args.flank_bed.read_text(),
                "chr1\t16\t20\t5:roo|ins1:FK_L\n"
                "chr1\t25\t29\t5:roo|ins1:FK_R\n"
                "chr1\t46\t50\t0:copia|del1:FK_L\n"
                "chr1\t50\t54\t0:copia|del1:FK_R\n",
            )
            self.assertEqual(
                args.formatted_flanks.read_text(),
                "5:chr1:16-20:roo|ins1\troo\tins1\tACGT\tCGTA\t.\t.\n"
                "0:chr1:46-50:copia|del1\tcopia\tdel1\tGTAC\tGTAC\t.\t.\n",
            )

    def test_zero_candidates_replace_stale_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fasta"
            genome.write_text(">chr1\nACGT\n")
            index = root / "genome.fasta.fai"
            index.write_text("chr1\t4\t6\t4\t5\n")
            insertion = root / "insertions.bed"
            deletion = root / "deletions.bed"
            insertion.write_text("")
            deletion.write_text("")
            outputs = [root / name for name in ("merged", "flank", "formatted")]
            for output in outputs:
                output.write_text("stale\n")
            args = SimpleNamespace(
                genome=genome,
                genome_index=index,
                insertion_bed=insertion,
                deletion_bed=deletion,
                flank_size=1,
                merged_bed=outputs[0],
                flank_bed=outputs[1],
                formatted_flanks=outputs[2],
            )
            indexed = FakeIndexedFasta({"chr1": "ACGT"})
            with patch(
                "prepare_insider_tsd.open_indexed_genome", return_value=indexed
            ):
                self.assertEqual(prepare(args), 0)
            self.assertEqual([path.read_bytes() for path in outputs], [b"", b"", b""])

    def test_rejects_invalid_flank_size(self):
        for value in ("0", "-1", "bad"):
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    positive_integer(value)


if __name__ == "__main__":
    unittest.main()
