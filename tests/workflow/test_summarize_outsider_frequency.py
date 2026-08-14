#!/usr/bin/env python3

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from summarize_outsider_frequency import (  # noqa: E402
    concatenate_support,
    main,
    precise_rows,
    read_support_counts,
    summarize_counts,
)


def evidence(
    evidence_type,
    read_name,
    event_id,
    sv_type,
    support="1",
    te="roo",
    position="chr1:100-100",
):
    fields = [
        evidence_type,
        read_name,
        event_id,
        te,
        support,
        "100",
        "100",
        ".",
        "100",
        position,
        "qseqid",
        "0",
        "90",
        "10",
        "10",
        "10",
        sv_type,
        ".",
        ".",
    ]
    return " ".join(fields) + "\n"


class FrequencyAggregationTests(unittest.TestCase):
    def test_aggregates_insertion_and_deletion_events(self):
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / "COUNT_READS.txt"
            counts.write_text(
                evidence("E", "reference", "event.INS.1", "INS")
                + evidence("H", "clipped", "event.INS.1", "INS")
                + evidence("I", "inserted", "event.INS.1", "INS")
                + evidence("E", "reference", "event.DEL.1", "DEL")
                + evidence("D", "deleted", "event.DEL.1", "DEL")
            )

            rows = summarize_counts(counts)
            by_event = {row[1]: row for row in rows}

            self.assertEqual(
                by_event["event.INS.1"],
                [
                    "chr1:100-100",
                    "event.INS.1",
                    "roo",
                    "1",
                    "1",
                    "2",
                    "2",
                    "3",
                    "50",
                    "66.6667",
                    "INS",
                ],
            )
            self.assertEqual(
                by_event["event.DEL.1"],
                [
                    "chr1:100-100",
                    "event.DEL.1",
                    "roo",
                    "1",
                    ".",
                    "2",
                    "1",
                    "2",
                    "50",
                    "50",
                    "DEL",
                ],
            )

    def test_empty_evidence_sorts_before_and_masks_insertion_for_the_same_read(self):
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / "COUNT_READS.txt"
            # I is intentionally written first: the compatibility sort still
            # moves E first for the same event/read pair.
            counts.write_text(
                evidence("I", "same-read", "event.INS.2", "INS", support="0")
                + evidence("E", "same-read", "event.INS.2", "INS", support="0")
            )

            rows = summarize_counts(counts)

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0][3], "0")
            self.assertEqual(rows[0][6:8], ["0", "1"])
            self.assertEqual(rows[0][9], "0")

    def test_pre_sorted_stream_preserves_the_same_legacy_state_machine(self):
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / "COUNT_READS.by_event.txt"
            counts.write_text(
                evidence("E", "same-read", "event.INS.2", "INS", support="0")
                + evidence("I", "same-read", "event.INS.2", "INS", support="0")
            )

            rows = summarize_counts(counts, pre_sorted=True)

            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0][6:8], ["0", "1"])

    def test_omits_hard_events_without_an_ins_or_del_type(self):
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / "COUNT_READS.txt"
            counts.write_text(evidence("H", "read", "HARD.4.R", "4"))

            self.assertEqual(summarize_counts(counts), [])

    def test_precise_correction_increments_totals_at_support_plus_one_boundary(self):
        frequency = [
            [
                "chr1:100-100",
                "sniffles.INS.1",
                "roo",
                "4",
                "2",
                "4",
                "6",
                "6",
                "100",
                "100",
                "INS",
            ]
        ]

        rows = precise_rows(frequency, {"sniffles.INS.1:roo": 5})

        self.assertEqual(
            rows,
            [
                [
                    "chr1:100-100",
                    "sniffles.INS.1",
                    "roo",
                    "5",
                    "2",
                    "5",
                    "7",
                    "7",
                    "100",
                    "100",
                    "INS",
                ]
            ],
        )


class SupportAndCliTests(unittest.TestCase):
    def test_support_files_are_concatenated_in_order_and_later_values_win(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            sniffles = root / "sniffles.txt"
            direct = root / "direct.txt"
            combined = root / "combined.txt"
            sniffles.write_bytes(b"shared 2\nsniffles.INS.1:roo 5\n")
            direct.write_bytes(b"shared 7\nTrEMOLO.INS.2:roo 3\n")

            support = read_support_counts((sniffles, direct))
            concatenate_support((sniffles, direct), combined)

            self.assertEqual(support["shared"], 7)
            self.assertEqual(
                combined.read_bytes(), sniffles.read_bytes() + direct.read_bytes()
            )

    def test_cli_writes_empty_outputs_for_empty_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            counts = root / "counts.txt"
            sniffles = root / "sniffles.txt"
            direct = root / "direct.txt"
            frequency = root / "frequency.tsv"
            precise = root / "precise.tsv"
            combined = root / "combined.txt"
            for path in (counts, sniffles, direct):
                path.write_text("")

            result = main(
                [
                    str(counts),
                    str(sniffles),
                    str(direct),
                    str(frequency),
                    str(precise),
                    "--combined-support-output",
                    str(combined),
                ]
            )

            self.assertEqual(result, 0)
            for path in (frequency, precise, combined):
                self.assertEqual(path.read_bytes(), b"")


if __name__ == "__main__":
    unittest.main()
