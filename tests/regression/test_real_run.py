import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location("real_run", Path(__file__).with_name("real_run.py"))
real_run = importlib.util.module_from_spec(spec)
spec.loader.exec_module(real_run)


class RealRunRegressionTests(unittest.TestCase):
    def test_cli_fails_on_changed_calls_and_refuses_to_replace_baseline(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            calls = root / "TE_INFOS.bed"
            calls.write_text("call\tstrand\na\t+\n")
            baseline = root / "baseline.json"
            command = [sys.executable, str(Path(__file__).with_name("real_run.py"))]
            created = subprocess.run(command + ["snapshot", "--run", str(root), "--output", str(baseline)], capture_output=True)
            self.assertEqual(created.returncode, 0, created.stderr)
            frozen = baseline.read_bytes()
            calls.write_text("call\tstrand\na\t-\n")
            checked = subprocess.run(command + ["check", "--baseline", str(baseline)], capture_output=True)
            self.assertEqual(checked.returncode, 1)
            replaced = subprocess.run(command + ["snapshot", "--run", str(root), "--output", str(baseline)], capture_output=True)
            self.assertNotEqual(replaced.returncode, 0)
            self.assertEqual(baseline.read_bytes(), frozen)

    def test_pipeline_revision_is_provenance_but_tool_versions_are_checked(self):
        expected = {"tool_versions": {"data": {"pipeline": {"git_commit": "old"}, "tools": {"samtools": "1.18"}}}}
        actual = {"tool_versions": {"data": {"pipeline": {"git_commit": "new"}, "tools": {"samtools": "1.18"}}}}
        self.assertEqual(real_run.differences(expected, actual), [])
        actual["tool_versions"]["data"]["tools"]["samtools"] = "1.19"
        self.assertTrue(real_run.differences(expected, actual))

    def test_table_order_is_ignored_but_call_sequence_and_orientation_changes_are_detected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "calls.tsv"
            path.write_text("call\tsequence\tstrand\na\tACGT\t+\nb\tTTTT\t-\n")
            before = real_run.signature(path, True)
            path.write_text("call\tsequence\tstrand\nb\tTTTT\t-\na\tACGT\t+\n")
            reordered = real_run.signature(path, True)
            self.assertEqual(real_run.differences(before, reordered), [])
            self.assertTrue(real_run.differences(before, reordered, strict=True))
            path.write_text("call\tsequence\tstrand\nb\tTTTT\t+\na\tACGT\t+\n")
            self.assertTrue(real_run.differences(before, real_run.signature(path, True)))

    def test_full_metadata_and_missing_metadata_are_distinguished(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "TE_INFOS.bed").write_text("a\tb\n1\t2\n")
            partial = real_run.capture(run=root)
            self.assertFalse(partial["reproducibility"]["complete"])
            config, versions, reads = (root / name for name in ("config.yaml", "versions.json", "reads.fastq"))
            config.write_text("DATA: {}\n")
            versions.write_text(json.dumps({"tools": {"samtools": {"returncode": 0, "stdout": "samtools 1.18"}}}))
            reads.write_text("@read\nACGT\n+\nIIII\n")
            full = real_run.capture(run=root, config=config, versions=versions, inputs=[("reads", reads)], command="snakemake all")
            self.assertTrue(full["reproducibility"]["complete"])
            self.assertEqual(full["configuration"]["text"], config.read_text())
            self.assertIn("TE_INFOS.bed", full["run"]["outputs"])

    def test_input_and_configuration_changes_are_detected_and_bad_versions_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "TE_INFOS.bed").write_text("a\tb\n1\t2\n")
            config = root / "config.yaml"
            config.write_text("threshold: 80\n")
            before = real_run.capture(run=root, config=config)
            config.write_text("threshold: 90\n")
            self.assertTrue(real_run.differences(before, real_run.capture(run=root, config=config)))
            versions = root / "versions.json"
            versions.write_text('{"tools":{"samtools":{"returncode":1}}}')
            with self.assertRaises(ValueError):
                real_run.capture(run=root, versions=versions)


if __name__ == "__main__":
    unittest.main()
