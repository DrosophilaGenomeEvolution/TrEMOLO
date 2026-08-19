#!/usr/bin/env python3

import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
import summarize_insider_frequency as frequency_summarizer  # noqa: E402


EXPECTED_HEADER = (
    "chrom\tposition\ttotal_depth\tdepth_empty_site\tread_support\t"
    "read_support_percent\tTE\tinfo_TE\n"
)


def write_inputs(root, insertions="", deletions="", depths=""):
    paths = (
        root / "INSERTION_TE.bed",
        root / "DEL_NB.bed",
        root / "DEPTH_FK.txt",
        root / "DEPTH_TE_INSIDER.csv",
    )
    paths[0].write_text(insertions)
    paths[1].write_text(deletions)
    paths[2].write_text(depths)
    return paths


def run_summarizer(paths, depth_margin=30):
    arguments = [str(path) for path in paths]
    arguments.extend(["--depth-margin", str(depth_margin)])
    return frequency_summarizer.main(arguments)


class SummarizeInsiderFrequencyTests(unittest.TestCase):
    def test_no_candidates_replaces_stale_output_with_header_only(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_inputs(Path(directory))
            paths[3].write_text("stale result\n")

            self.assertEqual(run_summarizer(paths), 0)

            self.assertEqual(paths[3].read_text(), EXPECTED_HEADER)

    def test_preserves_order_and_duplicates_and_uses_complete_candidate_name(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = write_inputs(
                Path(directory),
                insertions=(
                    "chr2\t100\t110\troo|same-event\n"
                    "chr1\t200\t210\tcopia|same-event\n"
                    "chr2\t100\t110\troo|same-event\n"
                ),
                deletions=(
                    "chr2\t100\t110\troo|same-event\t1\n"
                    "chr1\t200\t210\tcopia|same-event\t2\n"
                ),
                depths=(
                    "chr2:71\t3\n"
                    "chr2:140\t3\n"
                    "chr1:171\t3\n"
                    "chr1:240\t3\n"
                ),
            )

            self.assertEqual(run_summarizer(paths), 0)

            self.assertEqual(
                paths[3].read_text(),
                EXPECTED_HEADER
                + "chr2\t100\t3\t1\t2\t66.6666\troo\troo|same-event\n"
                + "chr1\t200\t3\t2\t1\t33.3333\tcopia\tcopia|same-event\n"
                + "chr2\t100\t3\t1\t2\t66.6666\troo\troo|same-event\n",
            )

    def test_bilateral_asymmetric_and_absent_depths_use_integer_mean(self):
        candidates = [
            ("chr1", 100, 110, "bilateral|event-1"),
            ("chr1", 200, 210, "asymmetric|event-2"),
            ("chr1", 300, 310, "absent|event-3"),
        ]
        depths = {
            "chr1:71": 2,
            "chr1:140": 5,
            "chr1:240": 3,
        }

        rows = frequency_summarizer.summarize(candidates, {}, depths)

        self.assertEqual([row[2] for row in rows], ["3", "1", "0"])
        self.assertEqual([row[4] for row in rows], ["3", "1", "0"])
        self.assertEqual([row[5] for row in rows], ["100.0000"] * 2 + ["0.0000"])

    def test_zero_integer_mean_resets_deletion_and_support_to_zero(self):
        candidates = [("chr1", 100, 110, "roo|one-sided")]
        deletions = {"roo|one-sided": 7}
        depths = {"chr1:140": 1}

        row = frequency_summarizer.summarize(candidates, deletions, depths)[0]

        self.assertEqual(row[2:6], ("0", "0", "0", "0.0000"))

    def test_legacy_percentage_truncates_positive_and_negative_ratios(self):
        cases = (
            (1, 3, "33.3333"),
            (2, 3, "66.6666"),
            (-1, 3, "-33.3333"),
            (-2, 3, "-66.6666"),
            (-1, 6, "-16.6666"),
            (0, 3, "0.0000"),
        )

        for support, total_depth, expected in cases:
            with self.subTest(support=support, total_depth=total_depth):
                self.assertEqual(
                    frequency_summarizer.legacy_percentage(support, total_depth),
                    expected,
                )

    def test_deletion_support_above_depth_is_not_clamped(self):
        candidates = [("chr1", 100, 110, "roo|negative-support")]
        deletions = {"roo|negative-support": 4}
        depths = {"chr1:71": 3, "chr1:140": 3}

        row = frequency_summarizer.summarize(candidates, deletions, depths)[0]

        self.assertEqual(row[2:6], ("3", "4", "-1", "-33.3333"))


class InsiderFrequencyInputValidationTests(unittest.TestCase):
    def test_malformed_records_are_rejected(self):
        readers_and_values = (
            (
                frequency_summarizer.read_insertions,
                "chr1\t100\t110\n",
                "at least four",
            ),
            (
                frequency_summarizer.read_insertions,
                "chr1\t110\t100\troo|event\n",
                "smaller than end",
            ),
            (
                frequency_summarizer.read_insertions,
                "chr1\t100\t110\tmissing-separator\n",
                "TE\\|event",
            ),
            (
                frequency_summarizer.read_deletion_counts,
                "chr1\t100\t110\troo|event\n",
                "at least five",
            ),
            (
                frequency_summarizer.read_deletion_counts,
                "chr1\t100\t110\troo|event\tnot-a-count\n",
                "must be an integer",
            ),
            (
                frequency_summarizer.read_depths,
                "malformed-key\t1\n",
                "malformed depth position",
            ),
            (
                frequency_summarizer.read_depths,
                "chr1:20\tnot-a-depth\n",
                "must be an integer",
            ),
        )

        for reader, content, message in readers_and_values:
            with self.subTest(reader=reader.__name__, content=content):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "input.tsv"
                    path.write_text(content)
                    with self.assertRaisesRegex(ValueError, message):
                        reader(path)

    def test_conflicting_duplicate_lookup_keys_are_rejected(self):
        cases = (
            (
                frequency_summarizer.read_deletion_counts,
                "chr1\t1\t2\troo|event\t1\n"
                "chr1\t1\t2\troo|event\t2\n",
                "conflicting deletion counts",
            ),
            (
                frequency_summarizer.read_depths,
                "chr1:10\t1\nchr1:10\t2\n",
                "conflicting depths",
            ),
        )

        for reader, content, message in cases:
            with self.subTest(reader=reader.__name__):
                with tempfile.TemporaryDirectory() as directory:
                    path = Path(directory) / "input.tsv"
                    path.write_text(content)
                    with self.assertRaisesRegex(ValueError, message):
                        reader(path)

    def test_identical_duplicate_lookup_keys_are_unambiguous(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            deletions = root / "deletions.bed"
            depths = root / "depths.tsv"
            deletions.write_text(
                "chr1\t1\t2\troo|event\t3\n"
                "chr1\t1\t2\troo|event\t3\n"
            )
            depths.write_text("chr1:10\t4\nchr1:10\t4\n")

            self.assertEqual(
                frequency_summarizer.read_deletion_counts(deletions),
                {"roo|event": 3},
            )
            self.assertEqual(
                frequency_summarizer.read_depths(depths), {"chr1:10": 4}
            )


if __name__ == "__main__":
    unittest.main()
