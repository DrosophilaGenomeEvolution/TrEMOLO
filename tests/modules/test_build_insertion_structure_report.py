#!/usr/bin/env python3

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/modules"))
from build_insertion_structure_report import build, build_events, write_browser_data  # noqa: E402


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
    def read_chunk(self, output, relative):
        text = (output / relative).read_text()
        prefix = f"window.trmStructureChunks[{json.dumps(relative)}]="
        self.assertTrue(text.startswith(prefix))
        return json.loads(text[len(prefix):-2])

    def test_browser_chunks_preserve_evidence_and_keep_it_out_of_quarto(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data = self.run_build(root)
            output = root / "report"
            self.assertEqual(json.loads((output / "report-data.json").read_text()), data)
            qmd = (output / "report.qmd").read_text()
            bootstrap = json.loads(qmd.split('<script id="structure-report-data" type="application/json">')[1].split('</script>')[0])
            self.assertNotIn("hsps", bootstrap)
            self.assertNotIn("structures", bootstrap)
            self.assertLess(len(qmd), 100_000)
            self.assertNotIn(data["structures"][0]["query_id"], qmd)
            rows = [row for path in bootstrap["index_files"]
                    for row in self.read_chunk(output, path)["rows"]]
            index = dict(zip(bootstrap["index_columns"], rows[0]))
            self.assertEqual(index["families"], ["copia", "low", "roo"])
            payload = self.read_chunk(output, index["evidence_file"])[index["query_id"]]
            self.assertEqual(payload["hsps"], data["hsps"])
            self.assertEqual(payload["matches"], data["matches"])
            self.assertEqual(payload["components"], data["components"])
            self.assertEqual(payload["structure"], data["structures"][0])
            # A rebuild cannot overwrite chunks still referenced by older HTML.
            rebuilt = write_browser_data(output, data, chunk_size=1)
            self.assertNotEqual(rebuilt["index_files"], bootstrap["index_files"])
            self.assertTrue((output / bootstrap["index_files"][0]).is_file())

    @unittest.skipUnless(shutil.which("node"), "Node.js is needed for the browser harness")
    def test_offline_loading_pagination_global_filters_and_missing_chunks(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast, sizes, queries, infos, candidates = self.fixture(root)
            base_blast = blast.read_text()
            base_query = queries.read_text().split("\t")[0]
            blast.write_text(base_blast * 20 + "".join(base_blast.replace(base_query, f"2L:<INS>:{100+i}:{101+i}:event{i+1}:1:IMPRECISE:1") for i in range(1, 251)))
            queries.write_text("".join(f"2L:<INS>:{100+i}:{101+i}:event{i+1}:1:IMPRECISE:1\t1000\n" for i in range(251)))
            build(blast, sizes, queries, root / "report", infos, candidates,
                  90, 50, 50, 0.2,
                  ROOT / "modules/2-MODULE_TE_BLAST/report.qmd",
                  ROOT / "modules/2-MODULE_TE_BLAST/report.css",
                  ROOT / "modules/2-MODULE_TE_BLAST/report.js", "test")
            result = subprocess.run(["node", str(ROOT / "tests/modules/test_structure_report_browser.js"),
                                     str(root / "report")], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_index_chunks_keep_all_queries_without_alignments(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast, sizes, queries, infos, candidates = self.fixture(root)
            blast.write_text("")
            queries.write_text("".join(f"2L:<INS>:{i}:{i+1}:event{i}:1:IMPRECISE:1\t1000\n"
                                       for i in range(2001)))
            data = build(blast, sizes, queries, root / "report", infos, candidates,
                         90, 50, 50, 0.2,
                         ROOT / "modules/2-MODULE_TE_BLAST/report.qmd",
                         ROOT / "modules/2-MODULE_TE_BLAST/report.css",
                         ROOT / "modules/2-MODULE_TE_BLAST/report.js", "test")
            output = root / "report"
            bootstrap = write_browser_data(output, data)
            self.assertEqual(len(bootstrap["index_files"]), 2)
            self.assertEqual(bootstrap["event_files"], [])
            rows = [dict(zip(bootstrap["index_columns"], row))
                    for path in bootstrap["index_files"]
                    for row in self.read_chunk(output, path)["rows"]]
            self.assertEqual({row["query_id"] for row in rows},
                             {row["query_id"] for row in data["structures"]})
            self.assertTrue(all(row["classification"] == "no_match" for row in rows))
            self.assertTrue(all(row["match_count"] == 0 for row in rows))
            for row in (rows[0], rows[-1]):
                evidence = self.read_chunk(output, row["evidence_file"])[row["query_id"]]
                self.assertEqual(evidence["hsps"], [])
                self.assertEqual(evidence["structure"]["query_id"], row["query_id"])

    def test_restored_names_match_raw_blast_subjects(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            blast, sizes, queries, infos, candidates = self.fixture(root)
            original = 'roo|variant: "name"'
            with infos.open(newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            rows[1][3] = original + '|event1'
            with infos.open('w', newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
            with candidates.open(newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            for row in rows[1:]:
                row[8] = original
                if row[9] == 'roo':
                    row[9] = original
            with candidates.open('w', newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
            mapping = root / 'INPUT/te_headers.tsv'
            mapping.parent.mkdir()
            with mapping.open('w', newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows([
                    ['original', 'prepared'], [original, 'roo'], ['copia', 'copia'], ['low', 'low']])
            data = build(blast, sizes, queries, root / 'report', infos, candidates,
                         80, 100, 20, 0.2,
                         ROOT / 'modules/2-MODULE_TE_BLAST/report.qmd',
                         ROOT / 'modules/2-MODULE_TE_BLAST/report.css',
                         ROOT / 'modules/2-MODULE_TE_BLAST/report.js', 'test')
            self.assertTrue(any(match['subject_te'] == original and match['assignment'] == 'reported_primary'
                                for match in data['matches']))

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
