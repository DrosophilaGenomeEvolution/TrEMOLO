#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "lib/python/workflow/annotate_resident_te.py"


def blast_row(
    te,
    chrom,
    qstart,
    qend,
    sstart,
    send,
    *,
    identity=90,
    bitscore=300,
    evalue="1e-30",
    qlen=1000,
    slen=10000,
):
    length = abs(qend - qstart) + 1
    return (
        te,
        chrom,
        str(identity),
        str(length),
        "0",
        "0",
        str(qstart),
        str(qend),
        str(sstart),
        str(send),
        str(evalue),
        str(bitscore),
        str(qlen),
        str(slen),
    )


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class AnnotateResidentTeTests(unittest.TestCase):
    def run_cli(self, root, rows, extra=()):
        blast = root / "hits.tsv"
        headers = root / "te_headers.tsv"
        with blast.open("w", newline="") as handle:
            csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
        names = sorted({row[0] for row in rows} | {"A", "B", "C", "D", "E"})
        with headers.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(("original", "prepared"))
            writer.writerows((name, name) for name in names)
        outputs = {
            label: root / filename
            for label, filename in (
                ("fragments", "fragments.tsv"),
                ("matches", "matches.tsv"),
                ("copies", "copies.tsv"),
                ("relations", "relations.tsv"),
                ("bed", "copies.bed"),
                ("gff", "copies.gff3"),
            )
        }
        command = [
            sys.executable,
            str(SCRIPT),
            str(blast),
            str(headers),
            "--target",
            "GENOME",
        ]
        command.extend(extra)
        for label, path in outputs.items():
            command.extend((f"--{label}", str(path)))
        result = subprocess.run(command, text=True, capture_output=True)
        return result, outputs

    def test_distinguishes_ambiguous_matches_from_nested_components(self):
        rows = [
            blast_row("A", "chr1", 1, 400, 101, 500),
            blast_row("A", "chr1", 401, 800, 501, 900),
            blast_row("B", "chr1", 1, 760, 121, 880, bitscore=550),
            blast_row("C", "chr1", 1, 151, 300, 450, qlen=500),
            blast_row("D", "chr1", 1, 101, 2001, 2101, qlen=1000),
            blast_row("E", "chr1", 1, 200, 3001, 3200, identity=60),
        ]
        with tempfile.TemporaryDirectory() as directory:
            result, outputs = self.run_cli(Path(directory), rows)
            self.assertEqual(result.returncode, 0, result.stderr)

            fragments = read_tsv(outputs["fragments"])
            matches = read_tsv(outputs["matches"])
            copies = read_tsv(outputs["copies"])
            relations = read_tsv(outputs["relations"])

            self.assertEqual(len(fragments), 6)
            self.assertEqual(
                [row["filter_status"] for row in fragments].count("rejected"), 1
            )
            self.assertEqual(len(matches), 4)
            self.assertEqual(len(copies), 3)

            ambiguous = next(row for row in copies if row["structure"] == "ambiguous")
            self.assertEqual(ambiguous["primary_te"], "A")
            self.assertEqual(ambiguous["te_candidates"], "A;B")
            self.assertEqual(ambiguous["candidate_count"], "2")
            a_match = next(row for row in matches if row["te_name"] == "A")
            b_match = next(row for row in matches if row["te_name"] == "B")
            self.assertEqual(a_match["assignment"], "primary")
            self.assertEqual(b_match["assignment"], "alternative")
            self.assertEqual(a_match["fragment_count"], "2")
            self.assertEqual(a_match["tier"], "full_length")
            d_copy = next(row for row in copies if row["primary_te"] == "D")
            self.assertEqual(d_copy["status"], "degraded_relic")

            self.assertEqual(len(relations), 1)
            self.assertEqual(relations[0]["relation"], "nested_candidate")
            self.assertEqual(relations[0]["confidence"], "provisional")
            bed = outputs["bed"].read_text().splitlines()
            self.assertEqual(len(bed), 3)
            self.assertTrue(all(len(line.split("\t")) == 6 for line in bed))
            self.assertEqual(len(outputs["gff"].read_text().splitlines()), 8)

    def test_chains_reverse_strand_collinearly(self):
        rows = [
            blast_row("A", "chr1", 401, 800, 1200, 801),
            blast_row("A", "chr1", 1, 400, 1600, 1201),
        ]
        with tempfile.TemporaryDirectory() as directory:
            result, outputs = self.run_cli(Path(directory), rows)
            self.assertEqual(result.returncode, 0, result.stderr)
            matches = read_tsv(outputs["matches"])
            self.assertEqual(len(matches), 1)
            self.assertEqual(matches[0]["strand"], "-")
            self.assertEqual(matches[0]["fragment_count"], "2")
            self.assertEqual(matches[0]["consensus_coverage"], "80")

    def test_ambiguity_grouping_is_complete_linkage_not_single_linkage(self):
        rows = [
            blast_row("A", "chr1", 1, 100, 101, 200),
            blast_row("B", "chr1", 1, 100, 121, 220),
            blast_row("C", "chr1", 1, 100, 141, 240),
        ]
        with tempfile.TemporaryDirectory() as directory:
            result, outputs = self.run_cli(Path(directory), rows)
            self.assertEqual(result.returncode, 0, result.stderr)
            copies = read_tsv(outputs["copies"])
            self.assertEqual(len(copies), 2)
            self.assertEqual(
                sorted(row["candidate_count"] for row in copies), ["1", "2"]
            )

    def test_empty_blast_replaces_every_stale_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in (
                "fragments.tsv",
                "matches.tsv",
                "copies.tsv",
                "relations.tsv",
                "copies.bed",
                "copies.gff3",
            ):
                (root / name).write_text("STALE\n")
            result, outputs = self.run_cli(root, [])
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(read_tsv(outputs["fragments"]), [])
            self.assertEqual(read_tsv(outputs["matches"]), [])
            self.assertEqual(read_tsv(outputs["copies"]), [])
            self.assertEqual(read_tsv(outputs["relations"]), [])
            self.assertEqual(outputs["bed"].read_bytes(), b"")
            self.assertEqual(outputs["gff"].read_text(), "##gff-version 3\n")

    def test_rejects_malformed_nonempty_blast(self):
        with tempfile.TemporaryDirectory() as directory:
            result, outputs = self.run_cli(Path(directory), [("too", "short")])
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("expected 14 BLAST columns", result.stderr)
            self.assertFalse(outputs["copies"].exists())


if __name__ == "__main__":
    unittest.main()
