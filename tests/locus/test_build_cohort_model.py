#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/locus"))
from build_cohort_model import build_cohort, cluster_calls  # noqa: E402
from build_locus_model import InputRun  # noqa: E402


HEADER = "\t".join(["#chrom", "start", "end", "TE|ID", "strand", "TSD", "pident", "psize_TE", "SIZE_TE", "NEW_POS", "FREQ", "FREQ_WITH_CLIPPED", "SV_SIZE", "ID_TrEMOLO", "TYPE"])


def row(position, te, call_id, strand="+"):
    return f"2L\t{position}\t{position + 1}\t{te}|{call_id}\t{strand}\tNONE\t98\t99\t9000\t{position}\t20\t25\t9010\tTE_ID_OUTSIDER.{call_id}\tINS"


class BuildCohortModelTests(unittest.TestCase):
    def test_same_site_different_te_become_two_alleles(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            first = root / "sample_a.bed"
            second = root / "sample_b.bed"
            first.write_text(HEADER + "\n" + row(100, "roo", "call1") + "\n")
            second.write_text(HEADER + "\n" + row(103, "copia", "call2", "-") + "\n")
            output = root / "cohort"
            counts = build_cohort(
                [InputRun("A", first), InputRun("B", second)], output, "ref", 10
            )
            self.assertEqual(counts["loci"], 1)
            self.assertEqual(counts["alleles"], 2)
            with (output / "components.tsv").open() as handle:
                components = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({item["te_name"] for item in components}, {"roo", "copia"})

    def test_complete_span_prevents_chaining(self):
        calls = [
            {"chrom": "2L", "start": str(pos), "new_position": str(pos), "sample_id": "S", "tremolo_id": str(pos)}
            for pos in (100, 109, 118)
        ]
        clusters = cluster_calls(calls, 10)
        self.assertEqual([len(cluster) for cluster in clusters], [2, 1])

    def test_incompatible_nearby_calls_can_be_kept_separate_by_window(self):
        calls = [
            {"chrom": "2L", "start": str(pos), "new_position": str(pos), "sample_id": "S", "tremolo_id": str(pos)}
            for pos in (100, 105)
        ]
        self.assertEqual(len(cluster_calls(calls, 0)), 2)


if __name__ == "__main__":
    unittest.main()
