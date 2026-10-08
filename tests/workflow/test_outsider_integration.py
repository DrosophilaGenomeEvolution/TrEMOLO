#!/usr/bin/env python3

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from integrate_outsider_te import build as build_integration  # noqa: E402
from prepare_outsider_liftoff import prepare as prepare_liftoff  # noqa: E402
from summarize_outsider_liftoff import build as build_liftoff  # noqa: E402
from summarize_outsider_liftoff import read_lifted_features, summarize  # noqa: E402
from summarize_outsider_liftoff import read_expected_insertions  # noqa: E402


class OutsiderIntegrationTests(unittest.TestCase):
    def write(self, root, name, text):
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return path

    def fixture(self, root):
        sniffles_qseqid = "chr1:<INS>:2:2:sniffles.INS.1:1:PRECISE:0:+"
        direct_qseqid = "chr1:<INS>:6:7:TrEMOLO.INS.2:1:IMPRECISE:1:-"
        output = root / "output"
        return SimpleNamespace(
            genome=self.write(root, "genome.fa", ">chr1\nAAAA\nAAAA\n"),
            te_database=self.write(root, "te.fa", ">roo\nACGA\n"),
            merged_bed=self.write(
                root,
                "merged.bed",
                "chr1\t6\t7\troo|TrEMOLO.INS.2\n"
                "chr1\t2\t3\troo|sniffles.INS.1\n"
                "chr1\t4\t5\troo|HARD.1.R\n",
            ),
            sniffles_calls=self.write(
                root,
                "sniffles.tsv",
                "sseqid\tqseqid\nroo\t{}\n".format(sniffles_qseqid),
            ),
            direct_calls=self.write(
                root,
                "direct.tsv",
                "sseqid\tqseqid\nroo\t{}\n".format(direct_qseqid),
            ),
            sniffles_fasta=self.write(
                root,
                "sniffles.fa",
                ">{}\nTT\n".format(sniffles_qseqid[:-2]),
            ),
            direct_fasta=self.write(
                root,
                "direct.fa",
                ">{}\nCCC\n".format(direct_qseqid[:-2]),
            ),
            canonical_genome=output / "pseudo.fa",
            observed_genome=output / "neo.fa",
            canonical_bed=output / "pseudo.bed",
            observed_bed=output / "neo.bed",
            canonical_public_bed=output / "pseudo-public.bed",
            observed_public_bed=output / "neo-public.bed",
            audit=output / "audit.tsv",
        )

    def test_integrates_in_coordinate_order_without_event_ids_in_dna(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))

            self.assertEqual(build_integration(args), 2)

            self.assertEqual(
                args.canonical_genome.read_text(),
                ">chr1\nAAACGAAA AATCGTAA\n".replace(" ", ""),
            )
            self.assertEqual(
                args.observed_genome.read_text(), ">chr1\nAATTAAAACCCAA\n"
            )
            self.assertEqual(
                args.canonical_bed.read_text(),
                "chr1\t2\t6\troo:sniffles.INS.1\n"
                "chr1\t10\t14\troo:TrEMOLO.INS.2\n",
            )
            self.assertEqual(
                args.canonical_public_bed.read_text(),
                "chr1\t2\t6\troo|sniffles.INS.1\n"
                "chr1\t10\t14\troo|TrEMOLO.INS.2\n",
            )
            sequence = "".join(args.canonical_genome.read_text().splitlines()[1:])
            self.assertNotIn("sniffles", sequence)
            self.assertEqual(len(args.audit.read_text().splitlines()), 3)

    def test_clipped_only_input_produces_unchanged_genomes_and_empty_beds(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            args.merged_bed.write_text("chr1\t4\t5\troo|SOFT.1.L\n")

            self.assertEqual(build_integration(args), 0)

            expected = ">chr1\nAAAAAAAA\n"
            self.assertEqual(args.canonical_genome.read_text(), expected)
            self.assertEqual(args.observed_genome.read_text(), expected)
            self.assertEqual(args.canonical_bed.read_text(), "")
            self.assertEqual(args.audit.read_text(), "\t".join((
                "event_id", "te_family", "source", "qseqid", "chromosome",
                "breakpoint", "strand", "canonical_length", "observed_length",
                "status", "reason",
            )) + "\n")

    def test_unresolved_call_is_rejected_without_hiding_valid_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            args.direct_fasta.write_text("")
            args.allow_partial = True

            self.assertEqual(build_integration(args), 1)

            rows = args.audit.read_text().splitlines()
            self.assertEqual(len(rows), 3)
            self.assertTrue(any("\tintegrated\tintegrated" in row for row in rows[1:]))
            self.assertTrue(
                any(
                    "\trejected\tmissing_observed_sequence" in row
                    for row in rows[1:]
                )
            )
            self.assertEqual(len(args.canonical_bed.read_text().splitlines()), 1)

    def test_strict_cli_fails_and_keeps_audit_without_publishing_genomes(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            args.direct_fasta.write_text("")
            command = [sys.executable, str(ROOT / "lib/python/workflow/integrate_outsider_te.py")]
            for key, value in vars(args).items():
                command.extend(["--" + key.replace("_", "-"), str(value)])
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("2 eligible, 1 integrated, 1 rejected", result.stderr)
            self.assertIn("missing_observed_sequence", args.audit.read_text())
            self.assertFalse(args.observed_genome.exists())
            self.assertFalse(args.canonical_genome.exists())

    def test_deletions_are_excluded_without_requiring_their_sequences(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            # A native caller ID may itself contain 'INS'; the event type is
            # the second identifier field, not a substring anywhere in the ID.
            args.merged_bed.write_text("chr1\t2\t4\troo|sniffles.DEL.INS.1\n")
            args.sniffles_fasta.write_text("")
            self.assertEqual(build_integration(args), 0)
            self.assertEqual(args.observed_genome.read_text(), ">chr1\nAAAAAAAA\n")
            with args.audit.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual((rows[0]["status"], rows[0]["reason"]), ("excluded", "non_insertion_event"))
            self.assertEqual(args.observed_bed.read_text(), "")

    def test_population_alternatives_at_same_site_are_concatenated(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            args.merged_bed.write_text(args.merged_bed.read_text().replace("chr1\t6\t7", "chr1\t2\t3"))
            args.direct_calls.write_text(args.direct_calls.read_text().replace("<INS>:6:7:", "<INS>:2:3:"))
            args.direct_fasta.write_text(args.direct_fasta.read_text().replace("<INS>:6:7:", "<INS>:2:3:"))
            self.assertEqual(build_integration(args), 2)
            self.assertEqual(args.observed_genome.read_text(), ">chr1\nAACCCTTAAAAAA\n")
            self.assertEqual(args.observed_bed.read_text(),
                             "chr1\t2\t5\troo:TrEMOLO.INS.2\nchr1\t5\t7\troo:sniffles.INS.1\n")

    def test_empty_sequence_prevents_a_complete_reconstruction(self):
        with tempfile.TemporaryDirectory() as directory:
            args = self.fixture(Path(directory))
            args.direct_fasta.write_text(args.direct_fasta.read_text().replace("CCC", ""))
            with self.assertRaisesRegex(ValueError, "Incomplete OUTSIDER integration"):
                build_integration(args)
            self.assertIn("empty_insertion_sequence", args.audit.read_text())


class OutsiderLiftoffTests(unittest.TestCase):
    def test_prepares_only_insertion_flanks_and_clamps_contig_edges(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            positions = root / "positions.bed"
            positions.write_text(
                "chr1\t5\t9\troo:sniffles.INS.1\n"
                "chr1\t20\t24\troo:sniffles.DEL.2\n"
            )
            index = root / "genome.fa.fai"
            index.write_text("chr1\t30\t6\t30\t31\n")
            gff = root / "features.gff"
            feature_file = root / "feature.txt"

            self.assertEqual(prepare_liftoff(positions, index, gff, feature_file, 10), 1)

            rows = [line.split("\t") for line in gff.read_text().splitlines()]
            self.assertEqual((rows[0][3], rows[0][4]), ("1", "5"))
            self.assertEqual((rows[1][3], rows[1][4]), ("10", "19"))
            self.assertEqual(feature_file.read_text(), "repeat_element\n")

    def test_rejects_discordant_chromosomes_from_public_calls(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lifted = root / "lifted.gff3"
            lifted.write_text(
                "##gff-version 3\n"
                "chr1\tLiftoff\trepeat_element\t1\t100\t.\t+\t.\t"
                "ID=event1;NAME=roo;SIDE=L;coverage=0.9\n"
                "chr1\tLiftoff\trepeat_element\t105\t200\t.\t+\t.\t"
                "ID=event1;NAME=roo;SIDE=R;coverage=0.8\n"
                "chr1\tLiftoff\trepeat_element\t1\t300\t.\t+\t.\t"
                "ID=event2;NAME=copia;SIDE=L;coverage=0.9\n"
                "chr2\tLiftoff\trepeat_element\t400\t500\t.\t+\t.\t"
                "ID=event2;NAME=copia;SIDE=R;coverage=0.8\n"
            )
            insider = root / "insider.bed"
            insider.write_text("chr1\t1\t1\told|event\n")
            positions = root / "positions.bed"
            positions.write_text("chr1\t10\t14\troo:event1\nchr1\t20\t24\tcopia:event2\n")
            args = SimpleNamespace(
                lifted_gff=lifted,
                positions=positions,
                insider_bed=insider,
                good_bed=root / "good.bed",
                bad_bed=root / "bad.bed",
                mapped_ids=root / "ids.txt",
                public_bed=root / "public.bed",
                combined_bed=root / "combined.bed",
                audit=root / "audit.tsv",
                max_gap=20000,
            )

            self.assertEqual(build_liftoff(args), (1, 1))

            self.assertEqual(
                args.public_bed.read_text(), "chr1\t100\t104\troo|event1\n"
            )
            self.assertNotIn("event2", args.public_bed.read_text())
            self.assertIn("event2", args.bad_bed.read_text())
            self.assertEqual(len(args.combined_bed.read_text().splitlines()), 2)
            self.assertIn("discordant_chromosomes", args.audit.read_text())

    def test_edge_and_internal_flanks_contain_only_genome_bases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "positions").write_text("chr1\t0\t2\troo:TrEMOLO.INS.1\n"
                                            "chr1\t8\t10\troo:TrEMOLO.INS.2\n"
                                            "chr1\t4\t6\troo:TrEMOLO.INS.3\n")
            (root / "index").write_text("chr1\t10\n")
            self.assertEqual(prepare_liftoff(root / "positions", root / "index", root / "gff", root / "features", 3), 3)
            rows = [line.split("\t") for line in (root / "gff").read_text().splitlines()]
            self.assertEqual([(int(r[3]), int(r[4])) for r in rows], [(3, 5), (6, 8), (2, 4), (7, 9)])
            self.assertTrue(all(1 <= int(r[3]) <= int(r[4]) <= 10 for r in rows))

    def project(self, root, left, right):
        path = root / "lift.gff"
        path.write_text("".join("chr1\tLiftoff\trepeat_element\t{}\t{}\t.\t{}\t.\t"
                                "ID=event;NAME=roo;SIDE={};coverage=1\n".format(start, end, strand, side)
                                for (start, end, strand), side in [(left, "L"), (right, "R")]))
        return summarize(read_lifted_features(path), 20000)

    def test_projection_uses_inner_boundaries_on_both_strands(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for left, right in [((1, 100, "+"), (106, 200, "+")),
                                ((106, 200, "-"), (1, 100, "-"))]:
                good, _, _, audit = self.project(root, left, right)
                self.assertEqual(good[0][1:3], (100, 105))
                self.assertEqual(audit[0]["projected_gap"], "5")
                self.assertEqual(audit[0]["projected_strand"], left[2])
            for left, right in [((1, 100, "+"), (101, 200, "+")),
                                ((101, 200, "-"), (1, 100, "-"))]:
                good, _, _, _ = self.project(root, left, right)
                self.assertEqual(good[0][1:3], (100, 100))

    def test_projection_rejects_opposite_strands_and_reversed_order(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            cases = [((100, 200, "+"), (190, 290, "-"), "discordant_flank_strands"),
                     ((106, 200, "+"), (1, 100, "+"), "discordant_flank_order"),
                     ((1, 100, "-"), (106, 200, "-"), "discordant_flank_order")]
            for left, right, reason in cases:
                good, _, _, audit = self.project(root, left, right)
                self.assertEqual(good, [])
                self.assertEqual(audit[0]["reason"], reason)

    def test_small_ordered_overlap_is_retained_for_tsd(self):
        with tempfile.TemporaryDirectory() as directory:
            good, _, _, audit = self.project(Path(directory), (1, 100, "+"), (98, 200, "+"))
            self.assertEqual(good[0][1:3], (97, 100))
            self.assertEqual(audit[0]["projected_gap"], "-3")

    def test_audit_includes_fully_unmapped_and_partially_mapped_insertions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "positions").write_text("chr1\t10\t14\troo:mapped\n"
                                            "chr1\t20\t24\tcopia:unmapped\n"
                                            "chr1\t30\t34\tblood:partial\n")
            (root / "lift.gff").write_text(
                "chr1\tLiftoff\trepeat_element\t1\t100\t.\t+\t.\tID=mapped;NAME=roo;SIDE=L\n"
                "chr1\tLiftoff\trepeat_element\t101\t200\t.\t+\t.\tID=mapped;NAME=roo;SIDE=R\n"
                "chr1\tLiftoff\trepeat_element\t201\t300\t.\t+\t.\tID=partial;NAME=blood;SIDE=L\n")
            expected = read_expected_insertions(root / "positions")
            good, _, ids, audit = summarize(read_lifted_features(root / "lift.gff"), 20000, expected)
            self.assertEqual(len(good), 1)
            self.assertEqual(ids, ["mapped"])
            self.assertEqual([row["event_id"] for row in audit], ["mapped", "unmapped", "partial"])
            self.assertEqual([(r["status"], r["reason"]) for r in audit],
                             [("projected", "concordant_flanks"), ("rejected", "no_mapped_flanks"),
                              ("rejected", "missing_or_ambiguous_flank")])
            self.assertEqual(audit[1]["te_family"], "copia")
            self.assertEqual(audit[1]["projected_chromosome"], "")

    def test_empty_liftoff_gff_keeps_every_event_in_cli_audit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            # This covers total mapping failure and an insertion occupying a
            # whole contig, for which no source flank can be generated.
            (root / "positions").write_text("chr1\t0\t10\troo:TrEMOLO.INS.1\n"
                                            "chr2\t4\t8\tcopia:TrEMOLO.INS.2\n")
            (root / "lift.gff").write_text("##gff-version 3\n")
            command = [sys.executable, str(ROOT / "lib/python/workflow/summarize_outsider_liftoff.py"),
                       "--lifted-gff", str(root / "lift.gff"), "--positions", str(root / "positions")]
            outputs = {name: root / name for name in ("good-bed", "bad-bed", "mapped-ids", "public-bed", "combined-bed", "audit")}
            for name, path in outputs.items():
                command.extend(["--" + name, str(path)])
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("0 insertions; 2 mappings rejected", result.stdout)
            with outputs["audit"].open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["te_family"] for row in rows], ["roo", "copia"])
            self.assertTrue(all(row["reason"] == "no_mapped_flanks" for row in rows))
            self.assertTrue(all(path.read_text() == "" for name, path in outputs.items() if name != "audit"))

    def test_no_expected_insertions_produces_an_empty_audit(self):
        self.assertEqual(summarize({}, 20000, {}), ([], [], [], []))

    def test_unknown_lifted_event_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as directory:
            good, _, _, _ = self.project(Path(directory), (1, 100, "+"), (101, 200, "+"))
            self.assertEqual(len(good), 1)
            with self.assertRaisesRegex(ValueError, "unexpected events: event"):
                summarize(read_lifted_features(Path(directory) / "lift.gff"), 20000, {})

    def test_duplicate_expected_event_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "positions"
            path.write_text("chr1\t1\t2\troo:event\nchr1\t3\t4\troo:event\n")
            with self.assertRaisesRegex(ValueError, "duplicate integrated event event"):
                read_expected_insertions(path)


if __name__ == "__main__":
    unittest.main()
