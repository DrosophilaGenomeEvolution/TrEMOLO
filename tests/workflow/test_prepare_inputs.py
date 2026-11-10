#!/usr/bin/env python3

import os
import gzip
import hashlib
import shutil
import subprocess
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

    def test_gzip_is_decompressed_only_when_requested_and_can_be_indexed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome, database, reads = root / "genome.fa", root / "te.fa", root / "reads.data"
            genome.write_text(">chr1\nACGT\n")
            database.write_text(">roo\nACGT\n")
            # No .gz suffix: compression must be detected from the stream.
            fastq = b"@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nJJJJ\n"
            reads.write_bytes(gzip.compress(fastq))
            original = reads.read_bytes()
            output = root / "output"
            linked = prepare(genome, database, output, sample=reads)
            prepared = output / "INPUT/reads.fastq"
            self.assertTrue(prepared.is_symlink())
            self.assertEqual(linked["sample_preparation"]["method"], "symlink")
            manifest = prepare(genome, database, output, sample=reads, decompress_sample=True)
            self.assertFalse(prepared.is_symlink())
            self.assertEqual(prepared.read_bytes(), fastq)
            self.assertEqual(reads.read_bytes(), original)
            self.assertEqual(manifest["inputs"]["sample"]["sha256"], hashlib.sha256(original).hexdigest())
            self.assertEqual(manifest["sample_preparation"],
                             {"source_compression": "gzip", "method": "decompressed"})
            if shutil.which("samtools"):
                subprocess.run(["samtools", "fqidx", str(prepared)], check=True, capture_output=True)
                result = subprocess.run(["samtools", "fqidx", str(prepared), "r2"],
                                        check=True, capture_output=True)
                self.assertEqual(result.stdout, b"@r2\nTGCA\n+\nJJJJ\n")

    def test_plain_fastq_stays_linked_even_when_decompression_requested(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome, database, reads = root / "genome.fa", root / "te.fa", root / "reads.fastq"
            genome.write_text(">chr1\nACGT\n")
            database.write_text(">roo\nACGT\n")
            reads.write_text("@r1\nAC\n+\nII\n")
            manifest = prepare(genome, database, root / "output", sample=reads, decompress_sample=True)
            self.assertTrue((root / "output/INPUT/reads.fastq").is_symlink())
            self.assertEqual(manifest["sample_preparation"],
                             {"source_compression": "none", "method": "symlink"})

    def test_truncated_gzip_does_not_publish_partial_fastq(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome, database, reads = root / "genome.fa", root / "te.fa", root / "reads.fastq.gz"
            genome.write_text(">chr1\nACGT\n")
            database.write_text(">roo\nACGT\n")
            original = gzip.compress(b"@r1\nAC\n+\nII\n")[:-4]
            reads.write_bytes(original)
            with self.assertRaises(EOFError):
                prepare(genome, database, root / "output", sample=reads, decompress_sample=True)
            self.assertFalse((root / "output/INPUT/reads.fastq").exists())
            self.assertEqual(list((root / "output/INPUT").glob(".reads-*")), [])
            self.assertEqual(reads.read_bytes(), original)
            self.assertFalse((root / "output/.inputs.prepared").exists())

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
