#!/usr/bin/env python3

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
            self.assertEqual((rows[1][3], rows[1][4]), ("9", "19"))
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
            args = SimpleNamespace(
                lifted_gff=lifted,
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
                args.public_bed.read_text(), "chr1\t100\t105\troo|event1\n"
            )
            self.assertNotIn("event2", args.public_bed.read_text())
            self.assertIn("event2", args.bad_bed.read_text())
            self.assertEqual(len(args.combined_bed.read_text().splitlines()), 2)
            self.assertIn("discordant_chromosomes", args.audit.read_text())


if __name__ == "__main__":
    unittest.main()
