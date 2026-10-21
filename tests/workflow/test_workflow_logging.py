import os
import json
from pathlib import Path
import subprocess
import shutil
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from workflow_log_error import diagnostic


class WorkflowLoggingTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("snakemake"), "requires Snakemake")
    def test_actual_sniffles_job_reports_filter_error_before_read_recovery(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            work = root / "analysis"
            vcf = work / "OUTSIDER/VARIANT_CALLING/SV.vcf"
            vcf.parent.mkdir(parents=True)
            vcf.write_text("##fileformat=VCFv4.1\n##contig=<ID=2L,length=1000>\n"
                           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for name in ("OUTSIDER/MAPPING/SAMPLE_mapping_GENOME_MD.sorted.bam",
                         "OUTSIDER/MAPPING/SAMPLE_mapping_GENOME_MD.sorted.bam.bai",
                         "OUTSIDER/TrEMOLO_SV_TE/INS/SV_SIZE.direct.tsv"):
                path = work / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("")
            genome = root / "genome.fa"; genome.write_text(">2L\nACGT\n")
            reads = root / "reads.fq"; reads.write_text("@r\nACGT\n+\nIIII\n")
            te = root / "te.fa"; te.write_text(">TE\nACGT\n")
            config = root / "config.json"
            config.write_text(json.dumps({"DATA": {"GENOME": str(genome), "SAMPLE": str(reads),
                "TE_DB": str(te), "WORK_DIRECTORY": str(work)},
                "CHOICE": {"PIPELINE": {"OUTSIDER_VARIANT": True}},
                "PARAMS": {"OUTSIDER_VARIANT": {"TE_DETECTION": {"CHROM_KEEP": r"^2L\s"}}}}))
            result = subprocess.run(["snakemake", "--snakefile", str(ROOT / "workflow/Snakefile"),
                "--configfile", str(config), "--directory", str(root), "--cores", "1",
                str(work / "OUTSIDER/VARIANT_CALLING/SEQUENCE_PROVENANCE.tsv"),
                "--allowed-rules", "extract_outsider_sniffles_sequences"], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            progress = (work / "log/progress.log").read_text()
            self.assertIn("[STEP] Reading sniffles1 VCF", progress)
            self.assertIn("failed (exit 2)", progress)
            self.assertIn("CHROM_KEEP matches no chromosome name", progress)
            self.assertNotIn("Recovering read insertions", progress)
            self.assertNotIn("[DONE]", progress)
            self.assertEqual(len(progress.splitlines()), 4)

    def execute(self, root, command, name="test_stage"):
        root.mkdir(parents=True, exist_ok=True)
        (root / "log").mkdir(exist_ok=True)
        environment = dict(os.environ, TRM_HELPER=str(ROOT / "lib/bash/workflow_logging.sh"),
                           TRM_WORK=str(root), TRM_LOG=str(root / "log/tool.log"), TRM_RULE=name)
        script = 'source "$TRM_HELPER"\ntrm_begin "$TRM_RULE" "Preparing sequences" "$TRM_WORK" "$TRM_LOG"\n' + command
        return subprocess.run(["bash", "-euo", "pipefail", "-c", script],
                              env=environment, capture_output=True, text=True)

    def test_success_keeps_data_stdout_clean_and_records_three_messages(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "analysis with spaces"
            result = self.execute(root, 'trm_step "Reading VCF"\nprintf ">query\\nACGT\\n"\n')
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout, ">query\nACGT\n")
            self.assertEqual(len(result.stderr.splitlines()), 3)
            self.assertIn("[analysis with spaces][test_stage][STEP] Reading VCF", result.stderr)
            self.assertIn("[DONE]", result.stderr)
            self.assertEqual((root / "log/progress.log").read_text(), result.stderr)

    def test_failure_in_redirected_subshell_surfaces_cause_and_keeps_exit_code(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = self.execute(root, 'trm_step "Reading VCF"\n'
                '(printf "Traceback (most recent call last):\\n  internal detail\\nValueError: CHROM_KEEP matches no chromosome name\\n"; exit 7) > "$TRM_LOG" 2>&1\n')
            self.assertEqual(result.returncode, 7, result.stderr)
            self.assertEqual(result.stdout, "")
            self.assertIn("Reading VCF failed (exit 7)", result.stderr)
            self.assertIn("ValueError: CHROM_KEEP matches no chromosome name", result.stderr)
            self.assertNotIn("internal detail", result.stderr)
            self.assertNotIn("[DONE]", result.stderr)
            self.assertEqual(len(result.stderr.splitlines()), 4)
            self.assertEqual((root / "log/progress.log").read_text(), result.stderr)

    def test_stale_or_missing_tool_log_does_not_invent_a_diagnostic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / "log").mkdir()
            log = root / "log/tool.log"; log.write_text("ValueError: stale error\n")
            previous = time.time() - 100
            os.utime(log, (previous, previous))
            result = self.execute(root, "exit 13\n")
            self.assertEqual(result.returncode, 13)
            self.assertNotIn("stale error", result.stderr)
            self.assertIn("Detailed log:", result.stderr)
            log.unlink()
            result = self.execute(root, "exit 14\n")
            self.assertEqual(result.returncode, 14)
            self.assertNotIn("[DONE]", result.stderr)

    def test_diagnostic_is_bounded_and_ignores_old_or_missing_files(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "huge.log"
            path.write_text("previous output\n" * 10000 + "\x1b[31mTypeError: " + "x" * 2000 + "\x1b[0m\n")
            result = diagnostic([path])
            self.assertTrue(result.startswith("TypeError: "))
            self.assertLessEqual(len(result), 400)
            self.assertNotIn("\x1b", result)
            self.assertEqual(diagnostic([root / "missing"]), "")
            self.assertEqual(diagnostic([path], time.time() + 10), "")

    def test_parallel_jobs_have_distinct_labels_in_shared_progress_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); (root / "log").mkdir()
            children = []
            for number in range(6):
                environment = dict(os.environ, TRM_HELPER=str(ROOT / "lib/bash/workflow_logging.sh"),
                                   TRM_WORK=str(root), TRM_RULE=f"stage{number}")
                children.append(subprocess.Popen(["bash", "-euo", "pipefail", "-c",
                    'source "$TRM_HELPER"; trm_begin "$TRM_RULE" "Running phase" "$TRM_WORK"; true'],
                    env=environment, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True))
            for child in children:
                stdout, stderr = child.communicate()
                self.assertEqual(child.returncode, 0, stderr)
                self.assertEqual(stdout, "")
            lines = (root / "log/progress.log").read_text().splitlines()
            self.assertEqual(len(lines), 12)
            for number in range(6):
                self.assertEqual(sum(f"[stage{number}][START]" in line for line in lines), 1)
                self.assertEqual(sum(f"[stage{number}][DONE]" in line for line in lines), 1)


if __name__ == "__main__":
    unittest.main()
