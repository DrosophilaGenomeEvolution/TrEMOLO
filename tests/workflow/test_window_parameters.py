import sys
import subprocess
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib/python/workflow'))
from window_parameters import DEFAULT_WINDOWS, resolve_windows
from build_te_infos import intersects_same_family_deletion


class WindowParametersTests(unittest.TestCase):
    def test_old_configs_keep_defaults_and_overrides_are_independent(self):
        self.assertEqual(resolve_windows(), DEFAULT_WINDOWS)
        values = resolve_windows({'OUTSIDER_FREQUENCY_WINDOW': 250})
        self.assertEqual(values['OUTSIDER_FREQUENCY_WINDOW'], 250)
        self.assertEqual(values['OUTSIDER_EVIDENCE_MATCH_WINDOW'], 100)
        self.assertEqual(DEFAULT_WINDOWS['OUTSIDER_FREQUENCY_WINDOW'], 50)
        self.assertIsNot(values, DEFAULT_WINDOWS)

    def test_invalid_distances_and_misspelled_options_fail(self):
        for value in [-1, True, 12.5, '50', None]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'non-negative integer'):
                resolve_windows({'OUTSIDER_FREQUENCY_WINDOW': value})
        with self.assertRaisesRegex(ValueError, 'Unknown PARAMS.WINDOWS'):
            resolve_windows({'OUTSIDER_FREQUENCY_WINDO': 50})
        with self.assertRaisesRegex(ValueError, 'must be a mapping'):
            resolve_windows([50])

    def test_zero_distance_is_supported(self):
        self.assertEqual(resolve_windows({'OUTSIDER_EVIDENCE_MATCH_WINDOW': 0})['OUTSIDER_EVIDENCE_MATCH_WINDOW'], 0)

    def test_deletion_exclusion_changes_at_configured_distance(self):
        position = {'chromosome': 'chr1', 'family': 'roo', 'start_value': 100, 'end_value': 101}
        deletions = {'chr1': [{'family': 'roo', 'start_value': 120, 'end_value': 121}]}
        self.assertFalse(intersects_same_family_deletion(position, deletions, 10))
        self.assertTrue(intersects_same_family_deletion(position, deletions, 30))
        deletions['chr1'][0]['family'] = 'copia'
        self.assertFalse(intersects_same_family_deletion(position, deletions, 100))

try:
    import pysam
except ImportError:
    pysam = None


@unittest.skipIf(pysam is None, 'requires pysam')
class BamWindowTests(unittest.TestCase):
    def make_bam(self, root):
        bam = root / 'reads.bam'
        with pysam.AlignmentFile(str(bam), 'wb', header={'HD': {'VN': '1.6', 'SO': 'coordinate'},
                                                       'SQ': [{'SN': 'chr1', 'LN': 1000}]}) as handle:
            read = pysam.AlignedSegment()
            read.query_name = 'support'
            read.query_sequence = 'A' * 20 + 'G' * 30 + 'A' * 20
            read.reference_id, read.reference_start = 0, 115
            read.mapping_quality, read.flag = 60, 0
            read.cigartuples = [(0, 20), (1, 30), (0, 20)]
            handle.write(read)
        pysam.index(str(bam))
        return bam

    def test_sniffles_fetch_matching_and_flank_windows_control_read_recovery(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = self.make_bam(root)
            bed = root / 'variants.bed'
            query = 'chr1:<INS>:100:101:sniffles.INS.1:1:PRECISE'
            bed.write_text('chr1\t100\t101\t' + query + '\t20\t' + 'G' * 20 + '\n')
            script = Path(__file__).resolve().parents[2] / 'lib/python/parsing/parse_bam_found_ins.py'
            base = [sys.executable, str(script), str(bam), str(bed)]
            narrow = subprocess.run(base + ['--fetch-window', '1', '--max_distance', '45'],
                                    check=True, capture_output=True, text=True).stdout
            self.assertEqual(sum(line.startswith('>') for line in narrow.splitlines()), 1)  # VCF consensus only
            distant = subprocess.run(base + ['--fetch-window', '40', '--max_distance', '30'],
                                     check=True, capture_output=True, text=True).stdout
            self.assertEqual(sum(line.startswith('>') for line in distant.splitlines()), 1)
            wide = subprocess.run(base + ['--fetch-window', '40', '--max_distance', '45', '--flank-size', '10'],
                                  check=True, capture_output=True, text=True).stdout
            self.assertEqual(sum(line.startswith('>') for line in wide.splitlines()), 2)
            self.assertEqual(wide.splitlines()[1], 'A' * 10 + 'G' * 30 + 'A' * 10)
            bed.write_text('chr1\t1\t2\t' + query + '\t20\t' + 'G' * 20 + '\n')
            subprocess.run(base + ['--fetch-window', '40'], check=True, capture_output=True, text=True)

    def test_frequency_window_controls_breakpoint_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bam = self.make_bam(root)
            query = 'chr1:<INS>:100:101:event.INS.1:1:PRECISE:0'
            candidates, sizes, sv_sizes = root / 'calls', root / 'sizes', root / 'sv_sizes'
            candidates.write_text('sseqid\tqseqid\nroo\t' + query + '\n')
            sizes.write_text('roo\t20\n')
            sv_sizes.write_text(query + '\t30\n')
            script = Path(__file__).resolve().parents[2] / 'lib/python/workflow/count_outsider_frequency.py'
            types = []
            for window in [20, 50]:
                output = root / ('counts' + str(window))
                subprocess.run([sys.executable, str(script), '--threads', '1', '--window', str(window),
                                str(bam), str(candidates), str(sizes), str(sv_sizes), str(output)],
                               check=True, capture_output=True, text=True)
                types.append([line.split()[0] for line in output.read_text().splitlines()])
            self.assertEqual(types, [['E'], ['I']])
