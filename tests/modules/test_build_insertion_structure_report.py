#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/modules"))
from build_insertion_structure_report import build, build_events  # noqa: E402


TE_HEADER = (
    "#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\t"
    "NEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\n"
)
CANDIDATE_HEADER = (
    "candidate_group_id\tsource\tchrom\tstart\tend\tevent_id\ttremolo_id\t"
    "event_type\treported_te\tcandidate_te\tassignment\tcandidate_rank\t"
    "candidate_count\tevidence_count\tevidence_fraction\tevidence_channels\t"
    "ambiguity_type\n"
)


class BuildInsertionStructureReportTests(unittest.TestCase):
    def fixture(self, root):
        query = "2L:<INS>:100:101:event1:1:IMPRECISE:1"
        blast = root / "matches.bln"
        blast.write_text(
            f"{query}\troo\t95\t400\t10\t1\t1\t400\t1\t400\t1e-100\t600\n"
            f"{query}\tcopia\t96\t300\t5\t0\t501\t800\t300\t1\t1e-90\t500\n"
            f"{query}\tlow\t50\t200\t80\t4\t20\t219\t1\t200\t1e-5\t50\n"
        )
        te_sizes = root / "TE_SIZE.tsv"
        te_sizes.write_text("roo\t500\ncopia\t400\nlow\t300\n")
        query_sizes = root / "SV_SIZE.tsv"
        query_sizes.write_text(f"{query}\t1000\n")
        te_infos = root / "TE_INFOS.bed"
        te_infos.write_text(
            TE_HEADER
            + "2L\t100\t101\troo|event1\t+\tNONE\t95\t80\t500\t100\t20\t20\t1000\tTE_ID_OUTSIDER.1\tINS\n"
        )
        candidates = root / "TE_CALL_CANDIDATES.tsv"
        candidates.write_text(
            CANDIDATE_HEADER
            + "TCA1\tOUTSIDER\t2L\t100\t101\tevent1\tTE_ID_OUTSIDER.1\tINS\troo\troo\treported_primary\t1\t2\t4\t0.8\tdirect\tunresolved_family_candidates\n"
            + "TCA1\tOUTSIDER\t2L\t100\t101\tevent1\tTE_ID_OUTSIDER.1\tINS\troo\tcopia\talternative\t2\t2\t1\t0.2\tdirect\tunresolved_family_candidates\n"
        )
        return blast, te_sizes, query_sizes, te_infos, candidates

    def run_build(self, root, **overrides):
        blast, te_sizes, query_sizes, te_infos, candidates = self.fixture(root)
        if overrides.pop("empty_blast", False):
            blast.write_text("")
        values = {
            "blast": blast,
            "te_sizes": te_sizes,
            "query_sizes": query_sizes,
            "output": root / "report",
            "te_infos": te_infos,
            "te_call_candidates": candidates,
            "min_pident": 90,
            "min_aligned_bp": 50,
            "min_consensus_coverage": 50,
            "max_component_overlap": 0.2,
            "template": ROOT / "modules/2-MODULE_TE_BLAST/report.qmd",
            "style": ROOT / "modules/2-MODULE_TE_BLAST/report.css",
            "script": ROOT / "modules/2-MODULE_TE_BLAST/report.js",
            "title": "test",
        }
        values.update(overrides)
        return build(**values)

    def test_preserves_rejected_hsps_and_proposes_non_overlapping_components(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = self.run_build(root)
            self.assertEqual(data["summary"]["hsps"], 3)
            self.assertEqual(data["summary"]["retained_hsps"], 2)
            self.assertEqual(data["summary"]["retained_matches"], 2)
            self.assertEqual(data["structures"][0]["classification"], "composite_candidate")
            self.assertEqual(data["events"][0]["classification"], "multi_te_candidate")
            self.assertEqual([row["te_name"] for row in data["components"]], ["roo", "copia"])
            matches = {row["subject_te"]: row for row in data["matches"]}
            self.assertEqual(matches["roo"]["assignment"], "reported_primary")
            self.assertEqual(matches["copia"]["assignment"], "known_alternative")
            self.assertEqual(matches["copia"]["subject_strand"], "-")
            self.assertEqual(matches["copia"]["query_start"], 500)
            self.assertEqual(matches["low"]["status"], "rejected")
            with (root / "report/insertion-hsps.tsv").open(newline="") as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 3)
            report_source = (root / "report/report.qmd").read_text()
            self.assertIn('id="trm-structure-gallery"', report_source)
            self.assertIn('id="trm-structure-te"', report_source)
            self.assertIn("Evidence details", report_source)

    def test_thresholds_change_matches_but_not_raw_hsp_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = self.run_build(root, min_pident=99)
            self.assertEqual(data["summary"]["hsps"], 3)
            self.assertEqual(data["summary"]["retained_hsps"], 0)
            self.assertEqual(data["summary"]["retained_matches"], 0)
            self.assertEqual(data["summary"]["components"], 0)
            self.assertEqual(data["structures"][0]["classification"], "unresolved")

    def test_empty_blast_produces_header_only_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = self.run_build(root, empty_blast=True)
            self.assertEqual(data["summary"]["queries"], 1)
            self.assertEqual(data["structures"][0]["classification"], "no_match")
            with (root / "report/insertion-structures.tsv").open(newline="") as handle:
                self.assertEqual(len(list(csv.reader(handle, delimiter="\t"))), 2)

    def test_two_query_sequences_promote_multi_te_event_support(self):
        structures = []
        components = []
        for index in (1, 2):
            structure_id = f"S{index}"
            structures.append(
                {
                    "structure_id": structure_id,
                    "event_id": "event",
                    "chrom": "2L",
                    "locus_start": 100,
                    "locus_end": 101,
                    "final_call": "yes",
                    "reported_te": "roo",
                    "retained_matches": 2,
                    "component_count": 2,
                    "classification": "composite_candidate",
                }
            )
            components.extend(
                [
                    {"structure_id": structure_id, "te_name": "roo"},
                    {"structure_id": structure_id, "te_name": "copia"},
                ]
            )
        event = build_events(structures, components)[0]
        self.assertEqual(event["multi_family_queries"], 2)
        self.assertEqual(event["classification"], "multi_te_supported")
        self.assertEqual(event["interpretation"], "coordinate_supported_across_queries")


if __name__ == "__main__":
    unittest.main()
