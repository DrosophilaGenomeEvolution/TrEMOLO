import io
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from chromosome_filter import compile_chromosome_filter, filter_rows


class ChromosomeFilterTests(unittest.TestCase):
    expression = r"^([23][RL]|[X4Y])$"

    def test_names_and_bed_first_column_have_identical_selection(self):
        names = ["2L", "2R", "3L", "3R", "X", "4", "Y", "2L_extra", "X_random", "scaffoldX", "2R2"]
        pattern = compile_chromosome_filter(self.expression)
        expected = names[:7]
        self.assertEqual([name for name in names if pattern.search(name)], expected)
        # The TE, read name and other columns cannot rescue an excluded chromosome.
        rows = [f"{name}\t1\t2\t2L|X\t2R\n" for name in names]
        output = io.StringIO()
        filter_rows(iter(rows), output, pattern)
        self.assertEqual(output.getvalue(), "".join(rows[:7]))

    def test_python_regex_and_legacy_comma_lists_use_the_same_engine(self):
        for expression in (self.expression, "^2L$,^2R$,^3L$,^3R$,^X$,^4$,^Y$"):
            pattern = compile_chromosome_filter(expression)
            self.assertIsNotNone(pattern.search("2L"))
            self.assertIsNone(pattern.search("2L_extra"))
        pattern = compile_chromosome_filter(r"^chr[1,2]{1,2}$,^X$")
        self.assertIsNotNone(pattern.search("chr12"))
        self.assertIsNotNone(pattern.search("X"))
        self.assertIsNone(pattern.search("chr123"))
        self.assertIsNotNone(compile_chromosome_filter(r"(?i)^x$,^2L$").search("X"))
        self.assertIsNotNone(compile_chromosome_filter(r"^a\,b$").search("a,b"))
        for expression in ("[", "", ",", "^X$,", None, [], 3):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                compile_chromosome_filter(expression)

    def test_streaming_cli_filters_stdin_and_files_preserving_rows(self):
        script = str(ROOT / "lib/python/workflow/chromosome_filter.py")
        rows = "2L\t1\t2\tTE\n2L_extra\t1\t2\t2L\n\n"
        result = subprocess.run([sys.executable, script, "--regex", self.expression],
                                input=rows, text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "2L\t1\t2\tTE\n")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.bed"
            path.write_text(rows)
            file_result = subprocess.run([sys.executable, script, "--regex", self.expression,
                                          "--input", str(path)], text=True, capture_output=True)
            self.assertEqual(file_result.returncode, 0, file_result.stderr)
            self.assertEqual(file_result.stdout, result.stdout)

    def write_vcf(self, root, caller, contigs=True):
        path = root / (caller + ".vcf")
        header = "##fileformat=VCFv4.1\n" if caller == "sniffles1" else "##fileformat=VCFv4.2\n##source=Sniffles2_2.8.1\n"
        if contigs:
            header += "##contig=<ID=2L>\n##contig=<ID=2L_extra,length=1000>\n"
        header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        for number, chrom in enumerate(("2L", "2L_extra"), 1):
            identifier = str(number) if caller == "sniffles1" else f"Sniffles2.INS.{number}"
            info = "PRECISE;SVTYPE=INS;END=10;SVLEN=4;RNAMES=r1;" + ("RE=2" if caller == "sniffles1" else "SUPPORT=2")
            header += f"{chrom}\t10\t{identifier}\tN\tACGT\t.\tPASS\t{info}\n"
        path.write_text(header)
        return path

    def test_both_sniffles_backends_use_names_for_sequences_and_read_names(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for caller in ("sniffles1", "sniffles2"):
                for contigs in (True, False):
                    path = self.write_vcf(root, caller, contigs)
                    output = root / "sequences.fa"
                    result = subprocess.run([sys.executable, str(ROOT / "lib/python/parsing/get_seq_vcf.py"),
                        "--caller", caller, "-c", self.expression, "-m", "0", str(path), str(output)],
                        capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    self.assertEqual(sum(line.startswith(">") for line in output.read_text().splitlines()), 1)
                    self.assertTrue(output.read_text().startswith(">2L:<INS>"))
                    read_dir = root / "reads"
                    result = subprocess.run([sys.executable, str(ROOT / "lib/python/parsing/extract_region_reads_vcf.py"),
                        "--caller", caller, "-c", self.expression, "-d", str(read_dir), str(path)],
                        capture_output=True, text=True)
                    self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                    paths = list(read_dir.iterdir())
                    self.assertEqual(len(paths), 1)
                    self.assertTrue(paths[0].name.startswith("reads_2L:"))
                    self.assertEqual(paths[0].read_text(), "r1\n")

    def test_no_matching_declared_contig_is_an_explicit_error_for_both_backends(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for caller in ("sniffles1", "sniffles2"):
                path = self.write_vcf(root, caller)
                for script in ("get_seq_vcf.py", "extract_region_reads_vcf.py"):
                    command = [sys.executable, str(ROOT / "lib/python/parsing" / script),
                               "--caller", caller, "-c", r"^2L\s"]
                    command += [str(path), str(root / "seq.fa")] if script == "get_seq_vcf.py" else ["-d", str(root / "reads"), str(path)]
                    result = subprocess.run(command, capture_output=True, text=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("matches no chromosome name", result.stderr)
                    self.assertIn("without whitespace", result.stderr)

    @unittest.skipUnless(shutil.which("snakemake") and shutil.which("bedtools"), "requires Snakemake and bedtools")
    def test_actual_merge_rule_filters_each_channel_by_name_only(self):
        import json
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome.fa"; genome.write_text(">2L\nACGT\n")
            reads = root / "reads.fastq"; reads.write_text("@r\nACGT\n+\nIIII\n")
            te = root / "te.fa"; te.write_text(">TE\nACGT\n")
            signatures = []
            for number, expression in enumerate((self.expression, "^2L$,^2R$,^3L$,^3R$,^X$,^4$,^Y$")):
                work = root / f"work{number}"
                files = ["OUTSIDER/TE_DETECTION/POSITION_START_TE.bed",
                         "OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.bed",
                         "OUTSIDER/TrEMOLO_SV_TE/SOFT/SOFT_TE.bed",
                         "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.bed"]
                for index, name in enumerate(files):
                    path = work / name; path.parent.mkdir(parents=True, exist_ok=True)
                    position = 100 + 500 * index
                    path.write_text(f"2L\t{position}\t{position+1}\tTE|event{index}\t+\t40\n"
                                    f"2L_extra\t{position}\t{position+1}\t2L|excluded{index}\t+\t40\n"
                                    f"unplaced\t{position}\t{position+1}\t2L|also_excluded{index}\t+\t40\n")
                config = root / f"config{number}.json"
                config.write_text(json.dumps({"DATA": {"GENOME": str(genome), "SAMPLE": str(reads),
                    "TE_DB": str(te), "WORK_DIRECTORY": str(work)},
                    "CHOICE": {"PIPELINE": {"OUTSIDER_VARIANT": True},
                               "OUTSIDER_VARIANT": {"CLIPPED_READS": True}},
                    "PARAMS": {"OUTSIDER_VARIANT": {"TE_DETECTION": {"CHROM_KEEP": expression}}}}))
                target = work / "INTERNAL/POSITION_TE_OUTSIDER.bed"
                result = subprocess.run(["snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
                    "--configfile", str(config), "--directory", str(root), "--cores", "1", str(target),
                    "--allowed-rules", "merge_outsider_te_evidence"], text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                rows = target.read_text().splitlines()
                self.assertEqual(len(rows), 4)
                self.assertTrue(all(row.startswith("2L\t") for row in rows))
                signatures.append(rows)
            self.assertEqual(signatures[0], signatures[1])


if __name__ == "__main__":
    unittest.main()
