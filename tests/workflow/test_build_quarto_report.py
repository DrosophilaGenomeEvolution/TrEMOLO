#!/usr/bin/env python3

import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from build_quarto_report import (  # noqa: E402
    build_report_data,
    proximity_groups,
    read_te_infos,
    render_source,
)


HEADER = (
    "#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\t"
    "NEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\n"
)


class BuildQuartoReportTests(unittest.TestCase):
    def write(self, root, name, content):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def fixture(self, root):
        te_infos = self.write(
            root,
            "TE_INFOS.bed",
            HEADER
            + "chr1\t100\t101\troo|out.1\t+\tATAT\t99\t98\t20\t100\t40\t50\t22\tTE_ID_OUTSIDER.1.INS.1\tINS\n"
            + "chr1\t180\t181\tcopia|in.1\t-\tNONE\t97\t96\t19\t180\tINSIDER\t25\t21\tTE_ID_INSIDER.2.INSERTION\tINSERTION\n"
            + "chr1\t250\t251\tblood|out.2\t+\tNONE\t95\t94\t18\t250\tNONE\tNONE\t20\tTE_ID_OUTSIDER.3.HARD.1\tHARD\n",
        )
        index = self.write(root, "genome.fasta.fai", "chr1\t1000\t6\t80\t81\n")
        manifest = self.write(
            root,
            "input_manifest.json",
            json.dumps({"inputs": {"genome": {"path": "genome.fa", "sha256": "abc"}}}),
        )
        mapping = self.write(
            root,
            "stats.txt",
            "raw total sequences:\t10\nreads mapped:\t9\nerror rate:\t1.2e-02\n",
        )
        vcf = self.write(
            root,
            "SV.vcf",
            "##fileformat=VCFv4.1\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t1\tN\t<INS>\t.\tPASS\tSVTYPE=INS;END=100\n"
            "chr1\t200\t2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=220\n",
        )
        return te_infos, index, manifest, mapping, vcf

    def test_builds_scientific_summaries_without_changing_calls(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = self.fixture(root)
            data = build_report_data(
                te_infos=paths[0],
                genome_index=paths[1],
                manifest_path=paths[2],
                mapping_stats_path=paths[3],
                sv_vcf_path=paths[4],
                enabled_sources=["INSIDER", "OUTSIDER"],
                title="Example",
                author="Tester",
                work_directory="work",
                locus_window=100,
            )

            self.assertEqual(data["summary"]["calls"], 3)
            self.assertEqual(data["summary"]["families"], 3)
            self.assertEqual(data["summary"]["sources"], {"INSIDER": 1, "OUTSIDER": 2})
            self.assertEqual(data["summary"]["tsd_confirmed"], 1)
            self.assertEqual(data["summary"]["frequency_available"], 2)
            self.assertEqual([call["display_frequency"] for call in data["calls"]], [50, 25, None])
            self.assertEqual(data["mapping_stats"]["reads mapped"], 9)
            self.assertEqual(data["sv_counts"], {"DEL": 1, "INS": 1})
            self.assertEqual(len(data["proximity_groups"]), 1)
            self.assertEqual(data["proximity_groups"][0]["families"], ["copia", "roo"])

    def test_proximity_grouping_does_not_single_linkage_chain(self):
        calls = [
            {"row": 1, "chrom": "chr1", "anchor": 100, "family": "A", "source": "OUTSIDER", "event_id": "a"},
            {"row": 2, "chrom": "chr1", "anchor": 180, "family": "B", "source": "OUTSIDER", "event_id": "b"},
            {"row": 3, "chrom": "chr1", "anchor": 250, "family": "C", "source": "OUTSIDER", "event_id": "c"},
        ]
        groups = proximity_groups(calls, 100)
        self.assertEqual(len(groups), 1)
        self.assertEqual(groups[0]["event_ids"], ["a", "b"])
        self.assertEqual(groups[0]["span"], 80)

    def test_header_only_input_has_a_valid_zero_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "TE_INFOS.bed"
            path.write_text(HEADER)
            self.assertEqual(read_te_infos(path), [])
            index = self.write(root, "genome.fasta.fai", "chr1\t1000\t6\t80\t81\n")
            manifest = self.write(root, "input_manifest.json", "{}\n")
            data = build_report_data(
                te_infos=path,
                genome_index=index,
                manifest_path=manifest,
                mapping_stats_path=None,
                sv_vcf_path=None,
                enabled_sources=[],
                title="Empty",
                author="",
                work_directory="work",
                locus_window=100,
            )
            self.assertEqual(data["summary"]["calls"], 0)
            self.assertEqual(data["summary"]["families"], 0)
            self.assertEqual(data["calls"], [])
            self.assertEqual(data["proximity_groups"], [])

    def test_rejects_malformed_rows(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "TE_INFOS.bed"
            path.write_text(HEADER + "chr1\tbad\t101\troo|id\t+\tNONE\t1\t2\t3\t4\t5\t6\t7\tid\tINS\n")
            with self.assertRaisesRegex(ValueError, "non-integer coordinates"):
                read_te_infos(path)

    def test_embeds_data_css_and_javascript_in_quarto_source(self):
        template = (
            "title: %%TREMOLO_REPORT_TITLE%%\n"
            "author: %%TREMOLO_REPORT_AUTHOR%%\n"
            "<!-- TREMOLO_REPORT_STYLE -->\n"
            "<!-- TREMOLO_REPORT_DATA -->\n"
            "<!-- TREMOLO_REPORT_SCRIPT -->\n"
        )
        data = {"report": {"title": 'A "title"', "author": "A <B>"}, "calls": []}
        source = render_source(template, "body {}", "console.log('ok')", data)
        self.assertIn('title: "A \\"title\\""', source)
        self.assertIn('author: "A <B>"', source)
        self.assertIn("<style>\nbody {}\n</style>", source)
        self.assertIn("console.log('ok')", source)
        self.assertIn('"author":"A \\u003cB>"', source)


if __name__ == "__main__":
    unittest.main()
