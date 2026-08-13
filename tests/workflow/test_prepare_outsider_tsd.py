#!/usr/bin/env python3

import sys
import tempfile
import unittest
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from prepare_outsider_tsd import (  # noqa: E402
    family_filename_component,
    prepare,
    read_fasta,
    validate_candidate,
)


def candidate(qseqid="chr1:sample:10:<INS>:event", qstart="6", qend="12"):
    return {
        "family": "roo",
        "qseqid": qseqid,
        "pident": "99",
        "size_per": "100",
        "size_el": "7",
        "qstart": qstart,
        "qend": qend,
        "sstart": "1",
        "send": "7",
        "source": "sniffles",
    }


class ReadFastaTests(unittest.TestCase):
    def test_reads_wrapped_sequences_and_uses_first_header_token(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "records.fasta"
            path.write_text(
                ">chr1 description ignored\nACGT\nTGCA\n\n>read2\naa\ncc\n"
            )

            records = read_fasta(path)

            self.assertEqual(list(records), ["chr1", "read2"])
            self.assertEqual(records["chr1"], "ACGTTGCA")
            self.assertEqual(records["read2"], "aacc")


class CandidateValidationTests(unittest.TestCase):
    def setUp(self):
        self.sequence = "A" * 24
        self.genome = {"chr1": 30}

    def validate(self, value, sequences=None, sizes=None, genome=None):
        sequences = ({value["qseqid"]: self.sequence} if sequences is None else sequences)
        sizes = ({value["qseqid"]: len(self.sequence)} if sizes is None else sizes)
        return validate_candidate(
            value,
            sequences,
            sizes,
            self.genome if genome is None else genome,
            10,
        )

    def test_reports_missing_sequence(self):
        status, reason, _ = self.validate(candidate(), sequences={})
        self.assertEqual((status, reason), ("rejected", "missing_sequence"))

    def test_reports_invalid_coordinates(self):
        status, reason, _ = self.validate(candidate(qstart="not-an-integer"))
        self.assertEqual((status, reason), ("rejected", "invalid_coordinates"))

    def test_reports_insufficient_left_flank(self):
        status, reason, _ = self.validate(candidate(qstart="4", qend="12"))
        self.assertEqual((status, reason), ("rejected", "insufficient_left_flank"))

    def test_reports_insufficient_right_flank(self):
        status, reason, _ = self.validate(candidate(qstart="6", qend="21"))
        self.assertEqual((status, reason), ("rejected", "insufficient_right_flank"))

    def test_reports_missing_genome_chromosome(self):
        value = candidate(qseqid="chrX:sample:10:<INS>:event")
        status, reason, _ = self.validate(value)
        self.assertEqual(
            (status, reason), ("rejected", "missing_genome_chromosome")
        )

    def test_reports_breakpoint_out_of_genome_bounds(self):
        value = candidate(qseqid="chr1:sample:31:<INS>:event")
        status, reason, _ = self.validate(value)
        self.assertEqual(
            (status, reason), ("rejected", "breakpoint_out_of_bounds")
        )

    def test_reports_query_interval_out_of_bounds(self):
        status, reason, _ = self.validate(candidate(qstart="6", qend="25"))
        self.assertEqual(
            (status, reason), ("rejected", "breakpoint_out_of_bounds")
        )

    def test_reports_insufficient_genome_left_flank(self):
        value = candidate(qseqid="chr1:sample:1:<INS>:event")
        status, reason, _ = self.validate(value)
        self.assertEqual(
            (status, reason), ("rejected", "insufficient_genome_left_flank")
        )

    def test_accepts_truncated_non_empty_genome_flanks_in_legacy_mode(self):
        value = candidate(qseqid="chr1:sample:2:<INS>:event")
        status, reason, _ = self.validate(value)
        self.assertEqual((status, reason), ("accepted", "accepted"))


class FakeIndexedFasta(object):
    def __init__(self, records):
        self.records = OrderedDict(records)
        self.references = tuple(self.records)
        self.lengths = tuple(len(sequence) for sequence in self.records.values())
        self.closed = False

    def fetch(self, name, start, end):
        return self.records[name][start:end]

    def close(self):
        self.closed = True


class PrepareOutsiderTsdTests(unittest.TestCase):
    def test_preparation_is_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            qseqid_dropped = "chr1:sample:15:<INS>:drop-first"
            qseqid_kept = "chr1:sample:20:<INS>:keep"
            qseqid_kept_second = "chr1:sample:25:<INS>:keep-second"
            genome = root / "genome.fasta"
            genome.write_text(">chr1\n" + "ACGT" * 10 + "\n")
            genome_index = root / "genome.fasta.fai"
            genome_index.write_text("chr1\t40\t6\t40\t41\n")
            sniffles_calls = root / "sniffles-calls.tsv"
            sniffles_calls.write_text(
                "family\tqseqid\n"
                "roo\t{}\n"
                "roo\t{}\n"
                "roo\t{}\n".format(
                    qseqid_dropped, qseqid_kept, qseqid_kept_second
                )
            )
            sniffles_combined = root / "sniffles-combined.tsv"
            sniffles_combined.write_text(
                "family\tqseqid\tpident\tsize_per\tsize_el\tqstart\tqend\tsstart\tsend\n"
                "roo\t{}\t99\t100\t5\t6\t10\t1\t5\n"
                "roo\t{}\t99\t100\t5\t6\t10\t1\t5\n"
                "roo\t{}\t99\t100\t5\t6\t10\t1\t5\n".format(
                    qseqid_dropped, qseqid_kept, qseqid_kept_second
                )
            )
            sniffles_fasta = root / "sniffles.fasta"
            sniffles_fasta.write_text(
                ">{}\nAAAAAAAAAA\nAAAAAAAAAA\n"
                ">{}\nCCCCCCCCCC\nCCCCCCCCCC\n".format(
                    qseqid_dropped, qseqid_kept
                )
            )
            with sniffles_fasta.open("a") as handle:
                handle.write(">{}\nGGGGGGGGGG\nGGGGGGGGGG\n".format(qseqid_kept_second))
            merged_bed = root / "merged.bed"
            merged_bed.write_text("")
            direct_calls = root / "direct-calls.tsv"
            direct_calls.write_text("family\tqseqid\n")
            direct_combined = root / "direct-combined.tsv"
            direct_combined.write_text(
                "family\tqseqid\tpident\tsize_per\tsize_el\tqstart\tqend\tsstart\tsend\n"
            )
            direct_fasta = root / "direct.fasta"
            direct_fasta.write_text("")
            sequence_sizes = root / "sequence-sizes.tsv"
            sequence_sizes.write_text(
                "{}\t20\n{}\t20\n{}\t20\n".format(
                    qseqid_dropped, qseqid_kept, qseqid_kept_second
                )
            )

            output = root / "output"
            args = SimpleNamespace(
                genome=genome,
                genome_index=genome_index,
                sniffles_calls=sniffles_calls,
                sniffles_combined=sniffles_combined,
                sniffles_fasta=sniffles_fasta,
                merged_bed=merged_bed,
                direct_calls=direct_calls,
                direct_combined=direct_combined,
                direct_fasta=direct_fasta,
                sequence_sizes=sequence_sizes,
                flank_size=10,
                te_fasta_dir=output / "te",
                read_artifact_dir=output / "reads",
                genome_artifact_dir=output / "genome",
                candidates=output / "candidates.tsv",
                validation=output / "validation.tsv",
                all_sequence_sizes=output / "sizes.tsv",
                read_flanks=output / "read-flanks.bed",
                genome_flanks=output / "genome-flanks.bed",
                combined_flanks=output / "combined-flanks.bed",
            )

            indexed = FakeIndexedFasta({"chr1": "ACGT" * 10})
            with patch(
                "prepare_outsider_tsd.open_indexed_genome", return_value=indexed
            ):
                self.assertEqual(prepare(args), (2, 2))
            self.assertTrue(indexed.closed)
            first = {
                path.relative_to(output): path.read_bytes()
                for path in sorted(output.rglob("*"))
                if path.is_file()
            }
            indexed = FakeIndexedFasta({"chr1": "ACGT" * 10})
            with patch(
                "prepare_outsider_tsd.open_indexed_genome", return_value=indexed
            ):
                self.assertEqual(prepare(args), (2, 2))
            self.assertTrue(indexed.closed)
            second = {
                path.relative_to(output): path.read_bytes()
                for path in sorted(output.rglob("*"))
                if path.is_file()
            }

            self.assertEqual(first, second)
            validation = args.validation.read_text().splitlines()
            self.assertEqual(len(validation), 3)
            self.assertIn("\taccepted\taccepted\t", validation[1])
            combined = args.combined_flanks.read_text().splitlines()
            self.assertEqual(len(combined), 2)
            self.assertIn("\tkeep\t", combined[0])

    def test_family_name_is_encoded_only_in_artifact_paths(self):
        self.assertEqual(family_filename_component("roo/sub family"), "roo%2Fsub%20family")


if __name__ == "__main__":
    unittest.main()
