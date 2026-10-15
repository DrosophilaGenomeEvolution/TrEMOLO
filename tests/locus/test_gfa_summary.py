import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/locus"))
from gfa_summary import refresh
from build_te_gfa import build
import test_build_te_gfa as fixtures


class GfaSummaryTests(unittest.TestCase):
    def test_call_counts_are_not_representative_counts_and_refresh_agrees(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = fixtures.GfaTests().fixture(root)
            out = root / "out"
            manifest = build(ref, db, inputs, out, min_aligned_bp=1)
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["calls"]["total"], 2)
            self.assertEqual(summary["calls"]["represented"], 2)
            self.assertEqual(summary["alternatives"]["representatives"], 4)
            self.assertEqual(manifest["counts"]["represented_calls"], 2)
            self.assertEqual(refresh(out, inputs), summary)
            self.assertTrue((out / "summary.html").is_file())

    def test_excluded_call_reasons_and_missing_source_bed_are_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = fixtures.GfaTests().fixture(root)
            out = root / "out"
            build(ref, db, inputs, out, min_aligned_bp=1, selected_tes={"C"})
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["calls"]["by_status"], {"represented": 1, "excluded_by_te_filter": 1})
            unavailable = refresh(out)
            self.assertFalse(unavailable["calls"]["available"])
            self.assertIsNone(unavailable["calls"]["total"])
            self.assertIn("unavailable", (out / "summary.html").read_text())

    def test_same_family_call_shares_representative_and_alternative_primary_is_visible(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, queries, blast = fixtures.GfaTests().fixture(root)
            infos = Path(inputs[0]["te_infos"])
            text = infos.read_text().replace("\t6\t7\tA|e2", "\t2\t3\tD|e2").replace("\t6\t50", "\t2\t50")
            infos.write_text(text)
            out = root / "out"
            build(ref, db, inputs, out, min_aligned_bp=1)
            summary = json.loads((out / "summary.json").read_text())
            self.assertEqual(summary["calls"]["by_representation"], {"own_sequence": 1, "same_locus_family": 1})
            self.assertEqual(summary["calls"]["alternative_family_only"], 1)
            with (out / "call-representation.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            second = next(r for r in rows if r["event_id"] == "e2")
            self.assertEqual(second["reported_te"], "D")
            self.assertEqual(second["represented_families"], "B")

    def test_html_escapes_original_family_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, blast = fixtures.GfaTests().fixture(root)
            family = "A<script>alert(1)</script>"
            db.write_text(db.read_text().replace(">A\n", ">" + family + "\n"))
            blast.write_text(blast.read_text().replace("\tA\t", "\t" + family + "\t"))
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            text = (root / "out/summary.html").read_text()
            self.assertNotIn("<script>", text)
            self.assertIn("&lt;script&gt;", text)


if __name__ == "__main__":
    unittest.main()
