#!/usr/bin/env python3

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from build_te_infos import (  # noqa: E402
    HEADER,
    build,
    legacy_chain_id,
    outsider_frequency_values,
)


class LegacyIdentifierTests(unittest.TestCase):
    def test_matches_chain_to_id_oracles(self):
        self.assertEqual(legacy_chain_id("2L_RaGOO_RaGOO:1408023"), "24486")
        self.assertEqual(legacy_chain_id("X_RaGOO_RaGOO:2989451"), "96331")
        self.assertEqual(legacy_chain_id(""), "0")

    def test_missing_outsider_frequency_uses_legacy_none_values(self):
        self.assertEqual(
            outsider_frequency_values("HARD.1.R", {}, {}), ("NONE", "NONE")
        )


class BuildTeInfosTests(unittest.TestCase):
    OUTSIDER_INPUTS = (
        "outsider_positions",
        "outsider_sniffles_combine",
        "outsider_direct_combine",
        "outsider_soft_calls",
        "outsider_hard_calls",
        "outsider_tsd",
        "outsider_frequency_precise",
        "outsider_frequency",
        "outsider_sv_sizes",
        "outsider_sniffles_calls",
        "outsider_direct_calls",
    )
    INSIDER_INPUTS = (
        "insider_deletion_bed",
        "insider_positions",
        "insider_insertion_combine",
        "insider_deletion_combine",
        "insider_tsd",
        "insider_frequency",
        "insider_variants",
    )

    def write(self, root, name, text=""):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def fixture(self, root):
        empty_header = "sseqid\tqseqid\tpident\tsize_per\tsize_el\n"
        args = SimpleNamespace(
            outsider_positions=self.write(
                root,
                "out-pos.bed",
                "chr1\t100\t101\troo|sniffles.INS.1\t50\t+\n",
            ),
            outsider_sniffles_combine=self.write(
                root,
                "out-combine.tsv",
                "sseqid\tqseqid\tgrain_pident\tsize_per\tsize_el\tqstart\tqend\tsstart\tsend\n"
                "roo\tchr1:<INS>:100:100:sniffles.INS.1:1:PRECISE:0:+\t99\t98\t20\t1\t21\t1\t20\n",
            ),
            outsider_direct_combine=self.write(root, "direct-combine.tsv", empty_header),
            outsider_soft_calls=self.write(root, "soft.tsv"),
            outsider_hard_calls=self.write(root, "hard.tsv"),
            outsider_tsd=self.write(
                root,
                "out-tsd.tsv",
                "roo|sniffles.INS.1\tATAT\t99\t22\tinfos\n",
            ),
            outsider_frequency_precise=self.write(
                root,
                "precise.tsv",
                "chr1:100-101\tsniffles.INS.1\troo\treadA\t0\t1\t1\t1\t50\t60\tINS\n",
            ),
            outsider_frequency=self.write(
                root,
                "frequency.tsv",
                "chr1:100-101\tsniffles.INS.1\troo\treadA\t0\t1\t1\t1\t40\t50\tINS\n",
            ),
            outsider_sv_sizes=self.write(
                root,
                "sizes.tsv",
                "chr1:<INS>:100:100:sniffles.INS.1:1:PRECISE:0\t25\n",
            ),
            outsider_sniffles_calls=self.write(
                root,
                "sniffles.tsv",
                "sseqid\tqseqid\tpident\tsize_per\tsize_el\n"
                "roo\tchr1:<INS>:100:100:sniffles.INS.1:1:PRECISE:0:+\t99\t98\t20\n",
            ),
            outsider_direct_calls=self.write(root, "direct.tsv", empty_header),
            insider_deletion_bed=self.write(root, "in-deletion.bed"),
            insider_positions=self.write(
                root,
                "in-pos.bed",
                "chr1\t100\t105\tcopia|Assemblytics_w_1\t5\t-\n",
            ),
            insider_insertion_combine=self.write(
                root,
                "in-combine.tsv",
                "sseqid\tqseqid\tgrain_pident\tsize_per\tsize_el\tqstart\tqend\tsstart\tsend\n"
                "copia\tAssemblytics_w_1:-:Repeat_expansion::ref:1-6:+\t97\t96\t5\t1\t6\t1\t5\n",
            ),
            insider_deletion_combine=self.write(root, "in-del-combine.tsv", empty_header),
            insider_tsd=self.write(root, "in-tsd.tsv"),
            insider_frequency=self.write(
                root,
                "in-freq.tsv",
                "chrom\tposition\ttotal_depth\tdepth_empty_site\tread_support\tread_support_percent\tTE\tinfo_TE\n"
                "chr1\t100\t2\t0\t2\t100.0000\tcopia\tcopia|Assemblytics_w_1\n",
            ),
            insider_variants=self.write(
                root,
                "variants.bed",
                "ref\t1\t1\tAssemblytics_w_1\t5\t+\tRepeat_expansion\t0\t7\tchr1:100-105:-\n",
            ),
            output=root / "TE_INFOS.bed",
        )
        return args

    def test_preserves_separate_calls_at_one_locus(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = self.fixture(root)

            self.assertEqual(build(args), 2)
            rows = [line.split("\t") for line in args.output.read_text().splitlines()]
            self.assertEqual(tuple(rows[0]), HEADER)
            self.assertEqual([row[3] for row in rows[1:]], [
                "roo|sniffles.INS.1",
                "copia|Assemblytics_w_1",
            ])
            self.assertEqual(rows[1][5:13], ["ATAT", "99", "98", "22", "99", "50", "60", "25"])
            self.assertEqual(rows[2][5:13], ["NONE", "97", "96", "5", "100", "100.0000", "INSIDER", "7"])
            self.assertEqual(rows[2][14], "Repeat_expansion")

    def test_header_only_output_replaces_stale_result(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            args = self.fixture(root)
            args.outsider_positions.write_text("")
            args.insider_positions.write_text("")
            args.output.write_text("stale\n")

            self.assertEqual(build(args), 0)
            self.assertEqual(args.output.read_text(), "\t".join(HEADER) + "\n")

    def test_disabled_pipeline_inputs_can_use_dev_null(self):
        for disabled, expected_name in (
            (self.OUTSIDER_INPUTS, "copia|Assemblytics_w_1"),
            (self.INSIDER_INPUTS, "roo|sniffles.INS.1"),
        ):
            with self.subTest(expected_name=expected_name):
                with tempfile.TemporaryDirectory() as directory:
                    args = self.fixture(Path(directory))
                    for attribute in disabled:
                        setattr(args, attribute, Path("/dev/null"))

                    self.assertEqual(build(args), 1)
                    rows = args.output.read_text().splitlines()
                    self.assertEqual(len(rows), 2)
                    self.assertEqual(rows[1].split("\t")[3], expected_name)


if __name__ == "__main__":
    unittest.main()
