#!/usr/bin/env python3

import os
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from prepare_inputs import prepare  # noqa: E402


class PrepareInputsTests(unittest.TestCase):
    def test_normalizes_fasta_and_creates_relative_links(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fa"
            genome.write_text(">chr1\nACGT\n")
            reads = root / "reads.fastq"
            reads.write_text("@r1\nAC\n+\nII\n")
            te_database = root / "te.fa"
            te_database.write_text(">roo\nac\ngt\n>copia\nttaa\n")
            output = root / "output"
            manifest = prepare(genome, te_database, output, sample=reads)
            self.assertFalse(manifest["pseudonyms_used"])
            self.assertEqual(manifest["te_records"], 2)
            self.assertEqual(
                (output / "INPUT/te_database.fasta").read_text(),
                ">roo\nACGT\n>copia\nTTAA\n",
            )
            self.assertTrue((output / "INPUT/genome.fasta").is_symlink())
            link_target = os.readlink(str(output / "INPUT/genome.fasta"))
            self.assertFalse(Path(link_target).is_absolute())

    def test_renames_all_headers_when_one_is_unsafe(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fa"
            genome.write_text(">chr1\nACGT\n")
            te_database = root / "te.fa"
            te_database.write_text(">safe\nAC\n>unsafe|name\nGT\n")
            output = root / "output"
            manifest = prepare(genome, te_database, output)
            self.assertTrue(manifest["pseudonyms_used"])
            self.assertEqual(
                (output / "INPUT/te_headers.tsv").read_text(),
                "original\tprepared\nsafe\tTrEMOLOTE1\nunsafe|name\tTrEMOLOTE2\n",
            )

    def test_rejects_invalid_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fa"
            genome.write_text(">chr1\nACGT\n")
            te_database = root / "te.fa"
            te_database.write_text("ACGT\n")
            with self.assertRaisesRegex(ValueError, "sequence before first FASTA header"):
                prepare(genome, te_database, root / "output")


if __name__ == "__main__":
    unittest.main()
