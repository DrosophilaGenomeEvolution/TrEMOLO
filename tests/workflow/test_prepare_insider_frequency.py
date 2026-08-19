#!/usr/bin/env python3

import argparse
import io
import sys
import tempfile
import unittest
from contextlib import redirect_stderr
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
import prepare_insider_frequency as frequency_preparer  # noqa: E402


HEADER = "sseqid\tqseqid\n"


def output_paths(root):
    return (
        root / "INSERTION_TE.bed",
        root / "FLANK_TE.bed",
        root / "FLANK_TE_IN.bed",
        root / "DEPTH_FK.bed",
    )


def run_preparer(candidates, outputs, extra_arguments=None):
    arguments = [str(candidates)] + [str(path) for path in outputs]
    arguments.extend(extra_arguments or [])
    return frequency_preparer.main(arguments)


class PrepareInsiderFrequencyTests(unittest.TestCase):
    def test_preserves_order_duplicates_and_exact_coordinates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "INSERTION.csv"
            candidates.write_text(
                HEADER
                + "roo\tAssemblytics_w_2:+:INSERTION::chr2:200-240:+\n"
                + "copia\tAssemblytics_b_1:-:Repeat_expansion::chr1:500-560:-\n"
                + "roo\tAssemblytics_w_2:+:INSERTION::chr2:200-240:+\n"
            )
            outputs = output_paths(root)

            self.assertEqual(run_preparer(candidates, outputs), 0)

            self.assertEqual(
                outputs[0].read_text(),
                "chr2\t200\t240\troo|Assemblytics_w_2\n"
                "chr1\t500\t560\tcopia|Assemblytics_b_1\n"
                "chr2\t200\t240\troo|Assemblytics_w_2\n",
            )
            self.assertEqual(
                outputs[1].read_text(),
                "chr2\t100\t101\n"
                "chr2\t340\t341\n"
                "chr1\t400\t401\n"
                "chr1\t660\t661\n"
                "chr2\t100\t101\n"
                "chr2\t340\t341\n",
            )
            self.assertEqual(
                outputs[2].read_text(),
                "chr2\t205\t206\n"
                "chr2\t235\t236\n"
                "chr1\t505\t506\n"
                "chr1\t555\t556\n"
                "chr2\t205\t206\n"
                "chr2\t235\t236\n",
            )
            self.assertEqual(
                outputs[3].read_text(),
                "chr2\t170\t171\n"
                "chr2\t269\t270\n"
                "chr1\t470\t471\n"
                "chr1\t589\t590\n"
                "chr2\t170\t171\n"
                "chr2\t269\t270\n",
            )

    def test_header_only_truncates_all_four_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "INSERTION.csv"
            candidates.write_text(HEADER)
            outputs = output_paths(root)
            for path in outputs:
                path.write_text("stale result\n")

            self.assertEqual(run_preparer(candidates, outputs), 0)

            self.assertEqual([path.read_bytes() for path in outputs], [b""] * 4)

    def test_invalid_qseqid_or_coordinates_fail_before_writing(self):
        invalid_values = (
            "malformed-qseqid",
            "Assemblytics_w_1:+:INSERTION::chr1:not-an-integer-200:+",
            "Assemblytics_w_1:+:INSERTION::chr1:200-200:+",
        )
        for qseqid in invalid_values:
            with self.subTest(qseqid=qseqid):
                with tempfile.TemporaryDirectory() as directory:
                    root = Path(directory)
                    candidates = root / "INSERTION.csv"
                    candidates.write_text(HEADER + "roo\t{}\n".format(qseqid))
                    outputs = output_paths(root)
                    for path in outputs:
                        path.write_text("keep me\n")

                    with self.assertRaises(ValueError):
                        run_preparer(candidates, outputs)

                    self.assertEqual(
                        [path.read_text() for path in outputs], ["keep me\n"] * 4
                    )

    def test_probes_left_of_the_contig_are_omitted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            candidates = root / "INSERTION.csv"
            candidates.write_text(
                HEADER + "roo\tAssemblytics_w_1:+:INSERTION::chr1:20-200:+\n"
            )
            outputs = output_paths(root)

            self.assertEqual(run_preparer(candidates, outputs), 0)

            self.assertEqual(
                [path.read_text() for path in outputs],
                [
                    "chr1\t20\t200\troo|Assemblytics_w_1\n",
                    "chr1\t300\t301\n",
                    "chr1\t25\t26\nchr1\t195\t196\n",
                    "chr1\t229\t230\n",
                ],
            )


class PrepareInsiderFrequencyCliTests(unittest.TestCase):
    def test_margin_options_accept_zero_and_reject_negative_values(self):
        parser = frequency_preparer.build_parser()
        positional = [
            "calls.tsv",
            "insertions.bed",
            "outer.bed",
            "inner.bed",
            "depth.bed",
        ]
        arguments = parser.parse_args(
            positional
            + [
                "--outer-margin",
                "0",
                "--inner-margin",
                "0",
                "--depth-margin",
                "0",
            ]
        )
        self.assertEqual(
            (arguments.outer_margin, arguments.inner_margin, arguments.depth_margin),
            (0, 0, 0),
        )

        for option in ("--outer-margin", "--inner-margin", "--depth-margin"):
            with self.subTest(option=option):
                with redirect_stderr(io.StringIO()):
                    with self.assertRaises(SystemExit):
                        parser.parse_args(positional + [option, "-1"])

        with self.assertRaises(argparse.ArgumentTypeError):
            frequency_preparer.non_negative_int("not-an-integer")


if __name__ == "__main__":
    unittest.main()
