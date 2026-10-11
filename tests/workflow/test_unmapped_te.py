import csv
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from unmapped_te import extract, annotate, resolve_parameters, MATCH_COLUMNS, HSP_COLUMNS


class UnmappedTETests(unittest.TestCase):
    def test_parameter_validation(self):
        for values in ({"MIN_PIDENT": True}, {"MIN_PIDENT": 101},
                       {"MIN_ALIGNED_BP": 0}, {"MIN_CHAIN_GAP": 1.5},
                       {"MAX_EVALUE": float("nan")}, {"AMBIGUITY_OVERLAP": 0},
                       {"TYPO": 5}, []):
            with self.subTest(values=values), self.assertRaises(ValueError):
                resolve_parameters(values)
        self.assertEqual(resolve_parameters()["MIN_PIDENT"], 80)

    def test_extraction_uses_original_sequences_and_excludes_any_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sam = root / "reads.sam"
            sam.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n"
                           "unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII\n"
                           "mixed\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n"
                           "mapped\t0\tchr1\t1\t0\t4M\t*\t0\t0\tAAAA\tIIII\n"
                           "mixed\t2048\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\n"
                           "secondary\t256\tchr1\t1\t60\t4M\t*\t0\t0\tAAAA\tIIII\n"
                           "secondary\t4\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII\n")
            fastq = root / "reads.fastq.gz"
            with gzip.open(fastq, "wt") as handle:
                for name, sequence in (("unmapped", "ACGT"), ("mapped", "AAAA"),
                                       ("mixed", "AAAA"), ("secondary", "AAAA")):
                    handle.write("@{} description\n{}\n+\nIIII\n".format(name, sequence))
            stats = extract(sam, fastq, root / "out.fa", root / "reads.tsv", root / "stats.json")
            self.assertEqual(stats, dict(total_reads=4, unmapped_candidates=3,
                                         excluded_with_mapped_alignment=2, unmapped_reads=1))
            self.assertEqual((root / "out.fa").read_text(), ">unmapped_000000000001\nACGT\n")
            self.assertIn("\tunmapped\t4\n", (root / "reads.tsv").read_text())

    def test_extraction_rejects_missing_or_duplicate_selected_reads(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sam = root / "reads.sam"
            sam.write_text("@HD\tVN:1.6\nq\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n")
            fastq = root / "reads.fastq"
            for content, reason in (("@other\nACGT\n+\nIIII\n", "absent"),
                                    ("@q\nACGT\n+\nIIII\n" * 2, "Duplicate")):
                fastq.write_text(content)
                with self.subTest(reason=reason), self.assertRaisesRegex(ValueError, reason):
                    extract(sam, fastq, root / "out.fa", root / "reads.tsv", root / "stats.json")

    def test_annotation_keeps_ambiguity_fragments_strands_and_repeated_copies(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "names.tsv").write_text('original\tprepared\n"TE A|variant: test"\tTE1\nTE B\tTE2\n')
            (root / "reads.tsv").write_text("query_id\tread_id\tread_length\nu1\tread-original\t2000\n")
            (root / "extract.json").write_text(json.dumps(dict(unmapped_reads=1, total_reads=1)))
            (root / "hits.bln").write_text(
                "TE1\tu1\t90\t200\t20\t0\t1\t200\t101\t300\t1e-30\t200\t500\t2000\n"
                "TE1\tu1\t90\t200\t20\t0\t301\t500\t401\t600\t1e-30\t200\t500\t2000\n"
                "TE2\tu1\t95\t500\t25\t0\t1\t500\t101\t600\t1e-30\t500\t500\t2000\n"
                "TE1\tu1\t90\t500\t50\t0\t1\t500\t1001\t1500\t1e-30\t450\t500\t2000\n"
                "TE1\tu1\t90\t100\t10\t0\t1\t100\t1900\t1801\t1e-20\t100\t500\t2000\n"
                "TE2\tu1\t40\t100\t60\t0\t1\t100\t1601\t1700\t1e-20\t50\t500\t2000\n")
            summary = annotate(root / "hits.bln", root / "names.tsv", root / "reads.tsv",
                               root / "extract.json", root / "matches.tsv", root / "hsps.tsv", root / "summary.json")
            with (root / "matches.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(summary["te_matches"], 4)
            self.assertEqual(summary["components"], 3)
            self.assertEqual(summary["ambiguous_components"], 1)
            fragmented = next(row for row in rows if row["fragment_count"] == "2")
            self.assertEqual(fragmented["te_family"], "TE A|variant: test")
            self.assertEqual(fragmented["read_segments"], "100-300;400-600")
            self.assertEqual(float(fragmented["te_coverage"]), 80)
            self.assertEqual(fragmented["read_aligned_bp"], "400")
            self.assertEqual(fragmented["candidate_count"], "2")
            self.assertEqual({row["assignment"] for row in rows[:2]}, {"primary", "alternative"})
            reverse = next(row for row in rows if row["strand"] == "-")
            self.assertEqual((reverse["read_start"], reverse["read_end"]), ("1800", "1900"))
            self.assertEqual(reverse["coverage_status"], "partial")
            with (root / "hsps.tsv").open() as handle:
                hsps = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(hsps), 6)
            self.assertEqual(sum(row["filter_status"] == "rejected" for row in hsps), 1)

    def test_empty_catalogue_has_headers_and_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "names.tsv").write_text("original\tprepared\nTE\tTE\n")
            (root / "reads.tsv").write_text("query_id\tread_id\tread_length\n")
            (root / "extract.json").write_text('{"unmapped_reads": 0, "total_reads": 0}')
            (root / "hits.bln").touch()
            summary = annotate(root / "hits.bln", root / "names.tsv", root / "reads.tsv",
                               root / "extract.json", root / "matches.tsv", root / "hsps.tsv", root / "summary.json")
            self.assertEqual(summary["reads_with_te"], 0)
            self.assertEqual((root / "matches.tsv").read_text(), "\t".join(MATCH_COLUMNS) + "\n")
            self.assertEqual((root / "hsps.tsv").read_text(), "\t".join(HSP_COLUMNS) + "\n")

    def test_short_repeated_te_fragments_are_separate_matches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "names.tsv").write_text("original\tprepared\nTE\tTE\n")
            (root / "reads.tsv").write_text("query_id\tread_id\tread_length\nu1\tr1\t300\n")
            (root / "extract.json").write_text('{"unmapped_reads": 1, "total_reads": 1}')
            (root / "hits.bln").write_text(
                "TE\tu1\t100\t80\t0\t0\t1\t80\t1\t80\t1e-20\t80\t80\t300\n"
                "TE\tu1\t100\t80\t0\t0\t1\t80\t101\t180\t1e-20\t80\t80\t300\n")
            summary = annotate(root / "hits.bln", root / "names.tsv", root / "reads.tsv",
                               root / "extract.json", root / "matches.tsv", root / "hsps.tsv", root / "summary.json")
            self.assertEqual(summary["te_matches"], 2)
            self.assertEqual(summary["components"], 2)


if __name__ == "__main__":
    unittest.main()
