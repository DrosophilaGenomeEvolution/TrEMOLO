#!/usr/bin/env python3

import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = ROOT / "lib/python/locus/build_locus_model.py"
SPEC = importlib.util.spec_from_file_location("build_locus_model", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


HEADER = "\t".join(
    [
        "#chrom", "start", "end", "TE|ID", "strand", "TSD", "pident",
        "psize_TE", "SIZE_TE", "NEW_POS", "FREQ", "FREQ_WITH_CLIPPED",
        "SV_SIZE", "ID_TrEMOLO", "TYPE",
    ]
)


class BuildLocusModelTests(unittest.TestCase):
    def test_conservative_projection_keeps_nearby_calls_separate(self):
        rows = [
            "2L\t100\t101\troo|call1\t+\tATAT\t98\t99\t9000\t100\t20\t25\t9010\tTE_ID_OUTSIDER.1\tINS",
            "2L\t105\t106\tcopia|call2\t-\tNONE\t97\t96\t5100\t105\t40\tNONE\t5110\tTE_ID_OUTSIDER.2\tINS",
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "TE_INFOS.bed"
            source.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
            output = root / "LOCUS_MODEL"
            counts = MODULE.build(
                [MODULE.InputRun("population_A", source)], output, "reference-test"
            )
            self.assertEqual(counts, {key: 2 for key in counts})
            manifest = json.loads((output / "manifest.json").read_text())
            self.assertEqual(manifest["construction_mode"], "one-call-per-locus")
            with (output / "loci.tsv").open() as handle:
                loci = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len({row["locus_id"] for row in loci}), 2)
            with (output / "observations.tsv").open() as handle:
                observations = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["frequency"] for row in observations], ["0.25", "0.4"])

    def test_ids_are_deterministic(self):
        self.assertEqual(
            MODULE.stable_id("L", "ref", "2L", 100, 101, "call"),
            MODULE.stable_id("L", "ref", "2L", 100, 101, "call"),
        )


if __name__ == "__main__":
    unittest.main()
