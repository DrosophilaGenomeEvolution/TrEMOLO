#!/usr/bin/env python3

import csv
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from build_te_call_candidates import build  # noqa: E402


TE_INFOS_HEADER = (
    "#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\t"
    "NEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\n"
)


class BuildTeCallCandidatesTests(unittest.TestCase):
    def write(self, root, name, content):
        path = root / name
        path.write_text(content)
        return path

    def test_writes_only_reported_calls_with_multiple_distinct_te_labels(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            infos = self.write(
                root,
                "TE_INFOS.bed",
                TE_INFOS_HEADER
                + "chr1\t100\t101\troo|event.1\t+\tNONE\t90\t80\t100\t100\t10\t10\t100\tTE_ID_OUTSIDER.1\tINS\n"
                + "chr1\t200\t201\tblood|event.2\t-\tNONE\t90\t80\t100\t200\t20\t20\t100\tTE_ID_OUTSIDER.2\tINS\n",
            )
            first = self.write(
                root,
                "first.txt",
                "event.1:roo\t2\n"
                "event.1:copia\t1\n"
                "event.2:blood\t5\n"
                "discarded:Doc\t3\n"
                "discarded:roo\t2\n",
            )
            second = self.write(root, "second.txt", "event.1:copia\t2\n")
            output = root / "TE_CALL_CANDIDATES.tsv"
            groups, rows = build(
                infos,
                [("direct", first), ("sniffles", second)],
                output,
            )
            self.assertEqual((groups, rows), (1, 2))
            with output.open(newline="") as handle:
                values = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["candidate_te"] for row in values], ["copia", "roo"])
            self.assertEqual([row["assignment"] for row in values], ["alternative", "reported_primary"])
            self.assertEqual([row["evidence_count"] for row in values], ["3", "2"])
            self.assertEqual(values[0]["evidence_channels"], "direct;sniffles")
            self.assertEqual(values[0]["candidate_count"], "2")
            self.assertEqual(values[0]["evidence_fraction"], "0.6")
            self.assertEqual(values[1]["evidence_fraction"], "0.4")
            self.assertEqual(values[0]["candidate_group_id"], values[1]["candidate_group_id"])

    def test_header_only_and_no_ambiguity_replace_stale_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            infos = self.write(root, "TE_INFOS.bed", TE_INFOS_HEADER)
            output = self.write(root, "TE_CALL_CANDIDATES.tsv", "STALE\n")
            groups, rows = build(infos, [], output)
            self.assertEqual((groups, rows), (0, 0))
            with output.open(newline="") as handle:
                values = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(len(values), 1)
            self.assertEqual(values[0][0], "candidate_group_id")

    def test_rejects_ambiguous_evidence_that_omits_reported_primary(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            infos = self.write(
                root,
                "TE_INFOS.bed",
                TE_INFOS_HEADER
                + "chr1\t1\t2\troo|event\t+\tNONE\t.\t.\t.\t1\t.\t.\t.\tTE_ID_OUTSIDER.1\tINS\n",
            )
            evidence = self.write(root, "counts.txt", "event:copia\t1\nevent:blood\t1\n")
            with self.assertRaisesRegex(ValueError, "reported TE roo is absent"):
                build(infos, [("direct", evidence)], root / "output.tsv")


if __name__ == "__main__":
    unittest.main()
