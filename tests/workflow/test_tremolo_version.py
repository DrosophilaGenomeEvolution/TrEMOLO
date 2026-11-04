import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('tremolo_version', ROOT / 'lib/python/tremolo_version.py')
version = importlib.util.module_from_spec(spec)
spec.loader.exec_module(version)


class SourceVersionTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which('git'), 'requires Git')
    def test_commit_and_dirty_state_belong_to_selected_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'VERSION').write_text('3.0.0-beta.1\n')
            def git(*args):
                return subprocess.check_output(['git', '-C', str(root)] + list(args), stderr=subprocess.DEVNULL, text=True).strip()
            git('init')
            git('add', 'VERSION')
            git('-c', 'user.name=Version test', '-c', 'user.email=test@example.invalid',
                '-c', 'commit.gpgsign=false', 'commit', '-m', 'fixture')
            identity = version.version_info(root)
            self.assertEqual(identity['commit'], git('rev-parse', 'HEAD'))
            self.assertFalse(identity['dirty'])
            (root / 'untracked.txt').write_text('uncommitted scientific change')
            self.assertTrue(version.version_info(root)['dirty'])
            archive = root / 'archive'
            archive.mkdir()
            (archive / 'VERSION').write_text('3.0.0-beta.1\n')
            # Do not attribute the parent repository's commit to an archive.
            self.assertIsNone(version.version_info(archive)['commit'])

    def test_archive_reports_unknown_commit_without_failing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'VERSION').write_text('3.0.0-beta.1\n')
            identity = version.version_info(root)
            self.assertIsNone(identity['commit'])
            self.assertIsNone(identity['dirty'])
            self.assertIn('unavailable', version.banner(identity))

    def test_cli_module_and_gfa_agree_without_scientific_inputs(self):
        record = json.loads(subprocess.check_output(
            [sys.executable, str(ROOT / 'lib/python/tremolo_version.py'), '--json'], text=True))
        for command in ([str(ROOT / 'module'), '--version'],
                        [sys.executable, str(ROOT / 'lib/python/locus/build_te_gfa.py'), '--version']):
            result = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(record['version'], result.stdout)
            self.assertIn(record['commit'] or 'unavailable', result.stdout)
