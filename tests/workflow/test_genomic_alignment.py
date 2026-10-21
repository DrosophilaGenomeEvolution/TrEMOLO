import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from genomic_alignment import alignment_metadata, ALIGNMENT_COLUMNS
from te_candidate_provenance import PROVENANCE_COLUMNS, read_provenance


class GenomicAlignmentTests(unittest.TestCase):
    def test_status_uses_only_sam_bits_and_keeps_combined_flags(self):
        for flag, status in ((0, "primary"), (16, "primary"), (256, "secondary"),
                             (272, "secondary"), (2048, "supplementary"),
                             (2304, "secondary_and_supplementary"), (4, "unmapped")):
            with self.subTest(flag=flag):
                self.assertEqual(alignment_metadata(flag), dict(
                    genomic_alignment_flag=str(flag), genomic_alignment_status=status))
        self.assertEqual(alignment_metadata()["genomic_alignment_status"], "unknown")
        self.assertEqual(alignment_metadata(sequence_origin="vcf_reported")["genomic_alignment_status"], "not_applicable")
        with self.assertRaises(ValueError):
            alignment_metadata(0, "vcf_reported")

    def test_old_sidecars_are_unknown_or_not_applicable(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "origins.tsv"
            path.write_text("\t".join(PROVENANCE_COLUMNS) + "\nq1\te\tr\tread_sequence\nq2\te\t.\tvcf_reported\n")
            rows = read_provenance(path)
            self.assertEqual(rows["q1"]["genomic_alignment_status"], "unknown")
            self.assertEqual(rows["q2"]["genomic_alignment_status"], "not_applicable")

    def test_inconsistent_or_invalid_flags_cannot_be_reported(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "origins.tsv"
            for flag, status in (("256", "primary"), ("-1", "primary"), ("65536", "primary")):
                with path.open("w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(PROVENANCE_COLUMNS + ALIGNMENT_COLUMNS)
                    writer.writerow(("q", "e", "r", "read_sequence", flag, status))
                with self.assertRaises(ValueError):
                    read_provenance(path)


if __name__ == "__main__":
    unittest.main()
