#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/modules"))
from build_population_report import build, read_sample_manifest  # noqa: E402


HEADER = (
    "#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\t"
    "NEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\n"
)


def call(position, family, call_id, frequency, strand="+"):
    return (
        f"2L\t{position}\t{position + 1}\t{family}|{call_id}\t{strand}\tNONE\t"
        f"98\t99\t9000\t{position}\t{frequency}\t{frequency}\t9010\t"
        f"TE_ID_OUTSIDER.{call_id}\tINS\n"
    )


class BuildPopulationReportTests(unittest.TestCase):
    def test_manifest_supports_replicates_and_relative_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("a", "b"):
                work = root / name
                work.mkdir()
                (work / "TE_INFOS.bed").write_text(HEADER)
            manifest = root / "samples.tsv"
            manifest.write_text(
                "sample_id\ttimepoint\treplicate\tte_infos\n"
                "A\tG10\t1\ta\n"
                "B\t20\t2\tb/TE_INFOS.bed\n"
            )
            samples = read_sample_manifest(manifest)
            self.assertEqual([sample.sample_id for sample in samples], ["A", "B"])
            self.assertEqual([sample.timepoint for sample in samples], [10.0, 20.0])
            self.assertEqual(samples[1].replicate, "2")
            self.assertEqual(samples[0].te_infos, root / "a/TE_INFOS.bed")

    def test_build_never_converts_missing_calls_to_zero(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "a.bed"
            second = root / "b.bed"
            first.write_text(HEADER + call(100, "roo", "a1", 20))
            second.write_text(
                HEADER
                + call(103, "roo", "b1", 60)
                + call(102, "copia", "b2", 40, "-")
            )
            manifest = root / "samples.tsv"
            manifest.write_text(
                "sample_id\ttimepoint\treplicate\tte_infos\n"
                f"A\tG1\t1\t{first}\n"
                f"B\tG2\t1\t{second}\n"
            )
            output = root / "report"
            data = build(
                manifest=manifest,
                output=output,
                reference_id="ref-test",
                locus_window=10,
                trend_epsilon=0.001,
                template=ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.qmd",
                style=ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.css",
                script=ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.js",
                title="test",
                genome=None,
            )
            self.assertEqual(data["summary"]["loci"], 1)
            self.assertEqual(data["summary"]["alleles"], 2)
            self.assertEqual(data["summary"]["observation_status"]["missing"], 1)
            with (output / "population-observations.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            missing = [row for row in rows if row["observation_status"] == "missing"]
            self.assertEqual(len(missing), 1)
            self.assertEqual(missing[0]["frequency"], ".")
            roo = next(row for row in data["trajectories"] if row["te_family"] == "roo")
            copia = next(row for row in data["trajectories"] if row["te_family"] == "copia")
            self.assertEqual(roo["trend"], "increasing")
            self.assertEqual(copia["trend"], "insufficient_data")
            self.assertIn('id="population-report-data"', (output / "report.qmd").read_text())

    def test_legacy_manifest_orders_explicit_generations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("old", "new"):
                work = root / name
                work.mkdir()
                (work / "TE_INFOS.bed").write_text(HEADER)
            manifest = root / "legacy.txt"
            manifest.write_text(f"{root / 'new'}:G20\n{root / 'old'}:G5\n")
            samples = read_sample_manifest(manifest)
            self.assertEqual([sample.sample_id for sample in samples], ["G5", "G20"])


if __name__ == "__main__":
    unittest.main()
