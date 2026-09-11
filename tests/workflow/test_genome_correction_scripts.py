#!/usr/bin/env python3

import os
import subprocess
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SAME_CHROMOSOME = ROOT / "lib/bash/correction_genome_transloc_sim_chrom.sh"
DIFFERENT_CHROMOSOMES = ROOT / "lib/bash/correction_genome_transloc_diff_chrom.sh"


SAMTOOLS_STUB = r"""#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) != 3 or sys.argv[1] != "faidx":
    raise SystemExit("samtools stub only supports: faidx FASTA")

fasta = Path(sys.argv[2])
records = []
name = None
sequence_lines = []
for raw_line in fasta.read_text().splitlines():
    if raw_line.startswith(">"):
        if name is not None:
            records.append((name, sequence_lines))
        name = raw_line[1:].split()[0]
        sequence_lines = []
    else:
        sequence_lines.append(raw_line.strip())
if name is not None:
    records.append((name, sequence_lines))
with Path(str(fasta) + ".fai").open("w") as output:
    for name, sequence_lines in records:
        line_lengths = [len(line) for line in sequence_lines if line]
        if len(line_lengths) > 1 and (
            any(length != line_lengths[0] for length in line_lengths[:-1])
            or line_lengths[-1] > line_lengths[0]
        ):
            raise SystemExit(f"Different line length in sequence '{name}'")
        sequence = "".join(sequence_lines)
        output.write(f"{name}\t{len(sequence)}\t0\t0\t0\n")
"""


BEDTOOLS_STUB = r"""#!/usr/bin/env python3
import sys
from pathlib import Path

if len(sys.argv) < 2 or sys.argv[1] != "getfasta":
    raise SystemExit("bedtools stub only supports getfasta")

arguments = sys.argv[2:]
fasta = Path(arguments[arguments.index("-fi") + 1])
bed = Path(arguments[arguments.index("-bed") + 1])
sequences = {}
name = None
parts = []
for raw_line in fasta.read_text().splitlines():
    if raw_line.startswith(">"):
        if name is not None:
            sequences[name] = "".join(parts)
        name = raw_line[1:].split()[0]
        parts = []
    else:
        parts.append(raw_line.strip())
if name is not None:
    sequences[name] = "".join(parts)

translation = str.maketrans("ACGTNacgtn", "TGCANtgcan")
for raw_line in bed.read_text().splitlines():
    chromosome, start, end, label, _score, strand = raw_line.split("\t")[:6]
    sequence = sequences[chromosome][int(start):int(end)]
    if strand == "-":
        sequence = sequence.translate(translation)[::-1]
    print(f">{label}::{chromosome}:{start}-{end}({strand})")
    print(sequence)
"""


MINIMAP2_STUB = r"""#!/usr/bin/env python3
import os
import sys
from pathlib import Path

arguments = sys.argv[1:]
marker = os.environ.get("MINIMAP_ARGUMENTS_FILE")
if marker:
    Path(marker).write_text("\n".join(arguments) + "\n")

def first_record(path):
    name = None
    sequence = []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                break
            name = line[1:].split()[0]
        elif name is not None:
            sequence.append(line.strip())
    return name, len("".join(sequence))

reference_name, reference_size = first_record(arguments[-2])
query_name, query_size = first_record(arguments[-1])
alignment_size = min(reference_size, query_size)
print(
    f"{query_name}\t{query_size}\t0\t{alignment_size}\t+\t"
    f"{reference_name}\t{reference_size}\t0\t{alignment_size}\t"
    f"{alignment_size}\t{alignment_size}\t60\ttp:A:P\tdv:f:0"
)
"""


def read_fasta(path):
    records = {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            records[name] = ""
        elif name is not None:
            records[name] += line.strip()
    return records


class GenomeCorrectionScriptTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        binary_directory = self.root / "bin"
        binary_directory.mkdir()
        for name, content in (
            ("samtools", SAMTOOLS_STUB),
            ("bedtools", BEDTOOLS_STUB),
            ("minimap2", MINIMAP2_STUB),
        ):
            path = binary_directory / name
            path.write_text(content)
            path.chmod(0o755)
        self.environment = os.environ.copy()
        self.environment["PATH"] = f"{binary_directory}:{self.environment['PATH']}"

    def tearDown(self):
        self.temporary_directory.cleanup()

    def write(self, name, content):
        path = self.root / name
        path.write_text(textwrap.dedent(content).lstrip())
        return path

    def run_script(self, script, *arguments):
        return subprocess.run(
            ["bash", str(script), *map(str, arguments)],
            cwd=self.root,
            env=self.environment,
            check=False,
            capture_output=True,
            text=True,
        )

    def common_options(self, pairs):
        return (
            "-c", pairs,
            "-s", "4",
            "--min-mapq", "10",
            "--min-identity", "90",
            "--max-anchor-overlap", "0",
        )

    def test_same_chromosome_reorders_reverses_and_preserves_unmapped_contig(self):
        reference = self.write("reference.fa", ">ref1\nCCCCGGAAAACC\n")
        query = self.write("query.fa", ">chr1\nAAAACCCCGGGG\n>unmapped\nTTTT\n")
        pairs = self.write("pairs.txt", "chr1:ref1\n")
        paf = self.write(
            "mapping.paf",
            """
            chr1\t12\t0\t4\t+\tref1\t12\t8\t12\t4\t4\t60\ttp:A:P
            chr1\t12\t8\t12\t-\tref1\t12\t0\t4\t4\t4\t60\ttp:A:P
            """,
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            SAME_CHROMOSOME,
            *self.common_options(pairs),
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            read_fasta(output),
            {"chr1": "CCCCGGAAAACC", "unmapped": "TTTT"},
        )
        audit = Path(str(output) + ".correction.corrections.tsv").read_text()
        self.assertIn("MOVE_AND_REVERSE", audit)
        self.assertIn("\tMOVE\t", audit)
        self.assertIn("unmapped\t0\t4\tunmapped", audit)

    def test_thread_default_and_g_option_are_forwarded_to_minimap2(self):
        reference = self.write("reference.fa", ">ref1\nAAAACCCC\n")
        query = self.write("query.fa", ">chr1\nAAAACCCC\n")
        pairs = self.write("pairs.txt", "chr1:ref1\n")

        for script in (SAME_CHROMOSOME, DIFFERENT_CHROMOSOMES):
            with self.subTest(script=script.name, threads="default"):
                marker = self.root / f"{script.stem}.default.args"
                self.environment["MINIMAP_ARGUMENTS_FILE"] = str(marker)
                result = self.run_script(
                    script,
                    "--report-only",
                    "-c", pairs,
                    "-s", "4",
                    reference,
                    query,
                    self.root / f"{script.stem}.default.fa",
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                arguments = marker.read_text().splitlines()
                self.assertEqual(arguments[arguments.index("-t") + 1], "8")

        marker = self.root / "same.custom.args"
        self.environment["MINIMAP_ARGUMENTS_FILE"] = str(marker)
        result = self.run_script(
            SAME_CHROMOSOME,
            "--report-only",
            "-g", "3",
            "-c", pairs,
            "-s", "4",
            reference,
            query,
            self.root / "same.custom.fa",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        arguments = marker.read_text().splitlines()
        self.assertEqual(arguments[arguments.index("-t") + 1], "3")

    def test_minimum_size_is_applied_before_correction(self):
        reference = self.write("reference.fa", ">ref1\nAAAACCCC\n")
        query = self.write("query.fa", ">chr1\nCCCCAAAA\n")
        pairs = self.write("pairs.txt", "chr1:ref1\n")
        paf = self.write(
            "mapping.paf",
            "chr1\t8\t0\t4\t+\tref1\t8\t4\t8\t4\t4\t60\ttp:A:P\n",
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            SAME_CHROMOSOME,
            "-c", pairs,
            "-s", "5",
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(read_fasta(output), {"chr1": "CCCCAAAA"})
        anchor_audit = Path(str(output) + ".correction.anchors.tsv").read_text()
        self.assertIn("query_span_below_minimum", anchor_audit)

    def test_minimap_divergence_tag_controls_identity_across_large_gaps(self):
        reference = self.write("reference.fa", ">ref1\nAAAACCCC\n")
        query = self.write("query.fa", ">chr1\nAAAACCCC\n")
        pairs = self.write("pairs.txt", "chr1:ref1\n")
        paf = self.write(
            "mapping.paf",
            # matches/block gives 50%, but minimap2 reports only 1%
            # sequence divergence; the anchor must pass a 90% identity filter.
            "chr1\t8\t0\t8\t+\tref1\t8\t0\t8\t4\t8\t60\t"
            "tp:A:P\tdv:f:0.01\n",
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            SAME_CHROMOSOME,
            "-c", pairs,
            "-s", "4",
            "--min-identity", "90",
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        anchor_rows = Path(
            str(output) + ".correction.anchors.tsv"
        ).read_text().splitlines()
        self.assertIn("\t99.0000\taccepted\taccepted", anchor_rows[1])
        selected_rows = Path(
            str(output) + ".correction.selected_anchors.tsv"
        ).read_text().splitlines()
        self.assertIn("\t99.0000\tselected\tselected", selected_rows[1])

    def test_unequal_source_blocks_are_rewrapped_as_valid_fasta(self):
        sequence = "AAACCCGGGTTTAAA"
        reference = self.write("reference.fa", f">ref1\n{sequence}\n")
        query = self.write("query.fa", f">chr1\n{sequence}\n")
        pairs = self.write("pairs.txt", "chr1:ref1\n")
        paf = self.write(
            "mapping.paf",
            """
            chr1\t15\t0\t2\t+\tref1\t15\t0\t2\t2\t2\t60\ttp:A:P
            chr1\t15\t4\t6\t+\tref1\t15\t4\t6\t2\t2\t60\ttp:A:P
            chr1\t15\t10\t15\t+\tref1\t15\t10\t15\t5\t5\t60\ttp:A:P
            """,
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            SAME_CHROMOSOME,
            "-c", pairs,
            "-s", "2",
            "--max-anchor-overlap", "0",
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(read_fasta(output), {"chr1": sequence})
        self.assertEqual(output.read_text().splitlines(), [">chr1", sequence])

    def test_different_chromosomes_swap_blocks_without_losing_sequence(self):
        reference = self.write(
            "reference.fa",
            ">refA\nAAAATTTT\n>refB\nGGGGCCCC\n",
        )
        query = self.write(
            "query.fa",
            ">chrA\nAAAACCCC\n>chrB\nGGGGTTTT\n>unmapped\nNNNN\n",
        )
        pairs = self.write("pairs.txt", "chrA:refA\nchrB:refB\n")
        paf = self.write(
            "mapping.paf",
            """
            chrA\t8\t0\t4\t+\trefA\t8\t0\t4\t4\t4\t60\ttp:A:P
            chrA\t8\t4\t8\t+\trefB\t8\t4\t8\t4\t4\t60\ttp:A:P
            chrB\t8\t0\t4\t+\trefB\t8\t0\t4\t4\t4\t60\ttp:A:P
            chrB\t8\t4\t8\t+\trefA\t8\t4\t8\t4\t4\t60\ttp:A:P
            """,
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            DIFFERENT_CHROMOSOMES,
            *self.common_options(pairs),
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(
            read_fasta(output),
            {"chrA": "AAAATTTT", "chrB": "GGGGCCCC", "unmapped": "NNNN"},
        )
        audit = Path(str(output) + ".correction.corrections.tsv").read_text()
        self.assertEqual(audit.count("\tTRANSFER\t"), 2)
        self.assertIn("unmapped\t0\t4\tunmapped", audit)

    def test_different_chromosomes_can_follow_reference_fasta_order(self):
        reference = self.write(
            "reference.fa",
            ">refA\nAAAAAAAA\n>refB\nCCCCCCCC\n",
        )
        query = self.write(
            "query.fa",
            ">chrB\nCCCCCCCC\n>unpaired\nNNNN\n>chrA\nAAAAAAAA\n",
        )
        pairs = self.write("pairs.txt", "chrA:refA\nchrB:refB\n")
        paf = self.write(
            "mapping.paf",
            """
            chrA\t8\t0\t8\t+\trefA\t8\t0\t8\t8\t8\t60\ttp:A:P\tdv:f:0
            chrB\t8\t0\t8\t+\trefB\t8\t0\t8\t8\t8\t60\ttp:A:P\tdv:f:0
            """,
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            DIFFERENT_CHROMOSOMES,
            "--reorder-chromosomes",
            "-c", pairs,
            "-s", "4",
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(list(read_fasta(output)), ["chrA", "chrB", "unpaired"])
        order_audit = Path(
            str(output) + ".correction.output_chromosomes.tsv"
        ).read_text().splitlines()
        self.assertEqual(
            order_audit[1:],
            [
                "1\tchrA\trefA\tchromosome_pair\t3",
                "2\tchrB\trefB\tchromosome_pair\t1",
                "3\tunpaired\t.\tunpaired\t2",
            ],
        )

    def test_large_cross_chromosome_anchor_beats_short_home_repeat(self):
        reference = self.write(
            "reference.fa",
            ">refA\nAAAACCCCGGGG\n>refB\nTTTTAAAACCCC\n",
        )
        query = self.write(
            "query.fa",
            ">chrA\nAAAACCCCGGGG\n>chrB\nTTTTGGGGCCCC\n",
        )
        pairs = self.write("pairs.txt", "chrA:refA\nchrB:refB\n")
        paf = self.write(
            "mapping.paf",
            """
            chrA\t12\t0\t4\t+\trefA\t12\t0\t4\t4\t4\t60\ttp:A:P\tdv:f:0
            chrA\t12\t8\t10\t+\trefA\t12\t8\t10\t2\t2\t60\ttp:A:P\tdv:f:0
            chrB\t12\t0\t4\t+\trefB\t12\t0\t4\t4\t4\t60\ttp:A:P\tdv:f:0
            chrB\t12\t4\t12\t+\trefA\t12\t4\t12\t8\t8\t60\ttp:A:P\tdv:f:0
            """,
        )
        output = self.root / "corrected.fa"

        result = self.run_script(
            DIFFERENT_CHROMOSOMES,
            "--report-only",
            "-c", pairs,
            "-s", "2",
            "--max-anchor-overlap", "0",
            reference,
            query,
            output,
            paf,
        )

        self.assertEqual(result.returncode, 0, result.stderr)
        correction_audit = Path(
            str(output) + ".correction.corrections.tsv"
        ).read_text()
        self.assertIn("chrB\t4\t12\tchrA", correction_audit)
        self.assertIn("\tTRANSFER\t8\t100.0000\t60", correction_audit)
        selection_audit = Path(
            str(output) + ".correction.selected_anchors.tsv"
        ).read_text()
        self.assertIn(
            "chrA\t8\t10\t+\trefA\t8\t10\t2\t60\t100.0000\t"
            "rejected\toverlapping_reference_anchor",
            selection_audit,
        )

    def test_duplicate_reference_pair_is_rejected(self):
        reference = self.write("reference.fa", ">ref1\nAAAAAAAA\n")
        query = self.write("query.fa", ">chr1\nAAAA\n>chr2\nAAAA\n")
        pairs = self.write("pairs.txt", "chr1:ref1\nchr2:ref1\n")
        paf = self.write(
            "mapping.paf",
            "chr1\t4\t0\t4\t+\tref1\t8\t0\t4\t4\t4\t60\ttp:A:P\n",
        )

        result = self.run_script(
            DIFFERENT_CHROMOSOMES,
            *self.common_options(pairs),
            reference,
            query,
            self.root / "corrected.fa",
            paf,
        )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("reference chromosome has two query pairs", result.stderr)


if __name__ == "__main__":
    unittest.main()
