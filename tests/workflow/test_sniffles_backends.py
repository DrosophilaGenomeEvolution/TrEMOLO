import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
sys.path.insert(0, str(ROOT / "lib/python/parsing"))
from run_sniffles import command
from sniffles_vcf import records
from integrate_outsider_te import build as integrate


class SnifflesTests(unittest.TestCase):
    def test_generation_specific_commands_and_version_mismatch(self):
        old = command("sniffles1", "sniffles", "1.0.12b", "a.bam", "g.fa", "a.vcf", 8, 1)
        self.assertEqual(old, ["sniffles", "-t", "8", "--report-seq", "-s", "1",
                               "-m", "a.bam", "-v", "a.vcf", "-n", "-1"])
        self.assertIn("--report_seq", command("sniffles1", "s", "1.0.10", "b", "g", "v", 1, 1))
        new = command("sniffles2", "sniffles2", "2.8.1", "a.bam", "g.fa", "a.vcf", 8, 3)
        self.assertIn("--output-rnames", new)
        self.assertIn("--reference", new)
        self.assertIn("--minsupport", new)
        self.assertNotIn("-m", new)
        for caller, version in [("sniffles1", "2.8.1"), ("sniffles2", "1.0.12b")]:
            with self.assertRaises(ValueError):
                command(caller, "s", version, "b", "g", "v", 1, 1)

    def fixture(self, root, rows):
        path = root / "input.vcf"
        path.write_text("##fileformat=VCFv4.2\n##source=Sniffles2_2.8.1\n"
                        "##contig=<ID=chr1,length=1000>\n"
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + rows)
        return path

    def row(self, identifier="Sniffles2.INS.1S0", ref="A", alt="AACGN", info=None):
        return "chr1\t10\t{}\t{}\t{}\t60\tPASS\t{}\n".format(identifier, ref, alt, info or
            "SVTYPE=INS;SVLEN=4;END=10;SUPPORT=3;IMPRECISE;RNAMES=r1,r2,r3")

    def test_anchor_removed_but_breakpoint_and_unknown_bases_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = self.fixture(root, self.row())
            record = list(records(path))[0]
            self.assertEqual((record["start"], record["end"], record["sequence"]), (10, 10, "ACGN"))
            self.assertEqual(record["precision"], "IMPRECISE")
            output = root / "seq.fa"
            subprocess.run([sys.executable, str(ROOT / "lib/python/parsing/get_seq_vcf.py"),
                            "-c", ".", "-m", "0", str(path), str(output)], check=True, capture_output=True)
            self.assertEqual(output.read_text(),
                ">chr1:<INS>:10:10:sniffles.INS.Sniffles2.INS.1S0:3:IMPRECISE\nACGN\n")
            subprocess.run([sys.executable, str(ROOT / "lib/python/parsing/extract_region_reads_vcf.py"),
                            "-c", ".", "-d", str(root / "reads"), str(path)], check=True, capture_output=True)
            self.assertEqual(next((root / "reads").iterdir()).read_text(), "r1\nr2\nr3\n")

    def test_unanchored_alt_does_not_lose_a_real_base(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(Path(directory), self.row(ref="N", alt="NACG"))
            self.assertEqual(list(records(path))[0]["sequence"], "NACG")

    def test_normalized_native_id_survives_genome_integration(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            record = list(records(self.fixture(root, self.row())))[0]
            identifier = record["id"]
            header = "chr1:<INS>:10:10:{}:3:IMPRECISE:0".format(identifier)
            def write(name, text):
                path = root / name
                path.write_text(text)
                return path
            args = SimpleNamespace(
                genome=write("genome.fa", ">chr1\n" + "A" * 30 + "\n"),
                te_database=write("te.fa", ">roo\nACGN\n"),
                merged_bed=write("merged.bed", "chr1\t10\t11\troo|" + identifier + "\n"),
                sniffles_calls=write("calls.tsv", "sseqid\tqseqid\nroo\t" + header + ":+\n"),
                direct_calls=write("direct.tsv", "sseqid\tqseqid\n"),
                sniffles_fasta=write("seq.fa", ">" + header + "\n" + record["sequence"] + "\n"),
                direct_fasta=write("empty.fa", ""),
                **{key: root / key for key in ("canonical_genome", "observed_genome",
                    "canonical_bed", "observed_bed", "canonical_public_bed", "observed_public_bed", "audit")})
            self.assertEqual(integrate(args), 1)
            self.assertEqual(args.observed_genome.read_text(), ">chr1\n" + "A" * 10 + "ACGN" + "A" * 20 + "\n")
            self.assertIn("chr1\t10\t14\troo|" + identifier, args.observed_public_bed.read_text())

    def test_deletion_and_symbolic_alleles(self):
        with tempfile.TemporaryDirectory() as directory:
            path = self.fixture(Path(directory), self.row(identifier="DEL.1", ref="AACGN", alt="A",
                info="PRECISE;SVTYPE=DEL;SVLEN=-4;END=14;SUPPORT=2") +
                self.row(identifier="DEL.2", ref="N", alt="<DEL>", info="SVTYPE=DEL;SVLEN=-9;SUPPORT=2") +
                self.row(identifier="INS.3", alt="<INS>"))
            result = list(records(path))
            self.assertEqual((result[0]["sequence"], result[0]["end"]), ("ACGN", 14))
            self.assertIsNone(result[1]["sequence"])
            self.assertEqual(result[1]["end"], 19)
            self.assertIsNone(result[2]["sequence"])

    def test_empty_vcf_is_a_valid_empty_fasta(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = self.fixture(root, "")
            output = root / "seq.fa"
            subprocess.run([sys.executable, str(ROOT / "lib/python/parsing/get_seq_vcf.py"),
                            "-c", ".", str(path), str(output)], check=True, capture_output=True)
            self.assertEqual(output.read_bytes(), b"")

    def test_malformed_lengths_and_duplicate_ids_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for rows in (self.row(alt="AC"), self.row() + self.row()):
                with self.subTest(rows=rows), self.assertRaises(ValueError):
                    list(records(self.fixture(root, rows)))

    def test_runner_propagates_failure_and_does_not_write_success_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            executable = root / "fake-sniffles"
            executable.write_text("#!/bin/sh\nif [ \"$1\" = --version ]; then echo 'Sniffles2, Version 2.8.1'; exit 0; fi\nexit 7\n")
            executable.chmod(0o755)
            result = subprocess.run([sys.executable, str(ROOT / "lib/python/workflow/run_sniffles.py"),
                "--caller", "sniffles2", "--executable", str(executable), "--bam", "a.bam",
                "--reference", "a.fa", "--vcf", str(root / "out.vcf"), "--metadata",
                str(root / "run.json"), "--threads", "1"], capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((root / "run.json").exists())


if __name__ == "__main__":
    unittest.main()
