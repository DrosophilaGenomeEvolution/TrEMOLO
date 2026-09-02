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
sequence = []
for raw_line in fasta.read_text().splitlines():
    if raw_line.startswith(">"):
        if name is not None:
            records.append((name, "".join(sequence)))
        name = raw_line[1:].split()[0]
        sequence = []
    else:
        sequence.append(raw_line.strip())
if name is not None:
    records.append((name, "".join(sequence)))
with Path(str(fasta) + ".fai").open("w") as output:
    for name, sequence in records:
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
        for name, content in (("samtools", SAMTOOLS_STUB), ("bedtools", BEDTOOLS_STUB)):
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
