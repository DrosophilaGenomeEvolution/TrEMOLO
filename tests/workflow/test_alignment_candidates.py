"""Boundary, grouping, and worker-error tests independent of the historical oracle."""
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

try:
    import pysam
except ImportError:
    pysam = None

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'lib/python/workflow'))
from cluster_alignment_candidates import format_insertions, format_hard


class ClusterFormattingTests(unittest.TestCase):
    def test_first_singleton_last_cluster_and_unique_support(self):
        rows = [
            ['chr1', '10', '11', 'a', 'AAA', '4', '4', '3', '1'],
            ['chr1', '20', '21', 'b', 'CCC', '4', '4', '3', '2'],
            ['chr1', '20', '21', 'b', 'GGG', '4', '4', '3', '2'],
            ['chr1', '20', '21', 'c', 'TTT', '4', '4', '3', '2'],
        ]
        fasta, positions, sizes = io.StringIO(), io.StringIO(), io.StringIO()
        format_insertions(iter(rows), fasta, positions, sizes)
        self.assertEqual(sum(line.startswith('>') for line in fasta.getvalue().splitlines()), 3)
        self.assertIn('TrEMOLO.INS.1:1:IMPRECISE:1\nAAA', fasta.getvalue())
        self.assertIn('TrEMOLO.INS.2:2:IMPRECISE:2\nTTT', fasta.getvalue())
        self.assertEqual(len(positions.getvalue().splitlines()), 3)
        self.assertEqual(len(sizes.getvalue().splitlines()), 3)

    def test_hard_singletons_first_and_last_are_retained(self):
        fasta = io.StringIO()
        format_hard(iter([
            ['chr1', '10', '11', 'a', 'AAA', 'L', '3', '1'],
            ['chr1', '20', '21', 'b', 'CCC', 'R', '3', '2'],
        ]), fasta)
        self.assertEqual(fasta.getvalue(), '>chr1:<HARD>:10:11:HARD.1.L:1:IMPRECISE:0\nAAA\n>chr1:<HARD>:20:21:HARD.2.R:1:IMPRECISE:0\nCCC\n')

    def test_empty_clusters(self):
        streams = [io.StringIO() for _ in range(3)]
        format_insertions(iter([]), *streams)
        self.assertTrue(all(stream.getvalue() == '' for stream in streams))


@unittest.skipIf(pysam is None, 'requires scientific container with pysam')
class AlignmentExtractionTests(unittest.TestCase):
    def run_extraction(self, root, records, threads=2, no_clipped=False, chunk_size=0):
        bam = root / 'input.bam'
        header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 100000}]}
        with pysam.AlignmentFile(str(bam), 'wb', header=header) as output:
            for name, start, cigar, sequence in records:
                read = pysam.AlignedSegment()
                read.query_name = name
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigartuples = cigar
                read.query_sequence = sequence
                output.write(read)
        pysam.index(str(bam))
        command = [sys.executable, str(ROOT / 'lib/python/parsing/find_all_type_ins.py'), str(bam),
                   '--threads', str(threads), '--chunk-size', str(chunk_size), '--output-ins', str(root / 'ins'),
                   '--output-soft', str(root / 'soft'), '--output-hard', str(root / 'hard'),
                   '--output-seq-tsd', str(root / 'tsd')]
        if no_clipped:
            command.append('--no-clipped')
        return subprocess.run(command, capture_output=True, text=True, timeout=20)

    def test_soft_grouping_hard_and_insertion_survive(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = self.run_extraction(root, [
                ('soft1', 100, [(4, 40), (0, 100)], 'A'*140),
                ('soft2', 110, [(4, 50), (0, 100)], 'C'*150),
                ('ins', 500, [(0, 50), (1, 40), (0, 50)], 'G'*140),
                ('hard', 700, [(5, 40), (0, 100)], 'T'*100),
            ])
            self.assertEqual(result.returncode, 0, result.stderr)
            soft = (root / 'soft').read_text().splitlines()
            self.assertEqual(len(soft), 1)
            self.assertIn('BEST_L_SIZE=50', soft[0])
            self.assertIn('NB_RS=2', soft[0])
            self.assertEqual(len((root / 'soft.bis').read_text().splitlines()), 3)
            self.assertEqual(len((root / 'hard').read_text().splitlines()), 2)
            self.assertEqual(len((root / 'ins').read_text().splitlines()), 1)
            self.assertFalse(list(root.glob('*.tmp.*')))

    def test_parallel_order_is_deterministic(self):
        records = [('a', 100, [(4, 40), (0, 100)], 'A'*140),
                   ('b', 110, [(4, 50), (0, 100)], 'C'*150),
                   ('c', 150, [(0, 90), (1, 40), (0, 20)], 'G'*150),
                   ('d', 500, [(5, 40), (0, 100)], 'T'*100)]
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            roots = [Path(first), Path(second)]
            for root, threads in zip(roots, [1, 4]):
                result = self.run_extraction(root, records, threads, chunk_size=105 if threads == 4 else 0)
                self.assertEqual(result.returncode, 0, result.stderr)
            for name in ['ins', 'soft', 'soft.bis', 'hard', 'tsd']:
                self.assertEqual((roots[0]/name).read_bytes(), (roots[1]/name).read_bytes())

    def test_no_clipped_truncates_previous_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ['soft', 'soft.bis', 'hard']:
                (root/name).write_text('stale evidence')
            result = self.run_extraction(root, [], no_clipped=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root/'ins').read_text(), '')
            for name in ['soft', 'soft.bis', 'hard']:
                self.assertEqual((root/name).read_text(), '')

    def test_timeout_fails_without_publishing_partial_output(self):
        from unittest.mock import patch
        import importlib.util
        spec = importlib.util.spec_from_file_location('alignment_extractor', ROOT/'lib/python/parsing/find_all_type_ins.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            self.run_extraction(root, [('a', 100, [(0, 100)], 'A'*100)])
            (root/'ins').write_text('previous successful output')
            command = ['extractor', str(root/'input.bam'), '--output-ins', str(root/'ins'), '--time-limit', '1']
            with patch.object(sys, 'argv', command), patch.object(module.multiprocessing, 'Pool') as pool:
                pool.return_value.__enter__.return_value.map_async.return_value.get.side_effect = module.multiprocessing.TimeoutError()
                with self.assertRaises(TimeoutError):
                    module.main()
            self.assertEqual((root/'ins').read_text(), 'previous successful output')
            self.assertFalse(list(root.glob('.alignment-candidates-*')))

    def test_worker_error_fails_without_publishing_partial_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = self.run_extraction(root, [('broken', 100, [(5, 40), (0, 100)], None)])
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse((root/'ins').exists())
            self.assertFalse(list(root.glob('*.tmp.*')))


@unittest.skipIf(pysam is None, 'requires scientific container with pysam')
class HardSequencePreparationTests(unittest.TestCase):
    def test_direct_fastq_extraction_preserves_order_slices_and_audit_indexes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'reads.fastq').write_text('@b\nTTTTCCCC\n+\nIIIIIIII\n@a\nAAAAGGGG\n+\nIIIIIIII\n')
            (root/'candidates').write_text('#header\nchr1\t10\t11\tHARD.0_0\ta\t6\t0\tL\t8\tSEQ\nchr1\t20\t21\tHARD.0_1\tb\t6\t0\tR\t8\tSEQ\n')
            result = subprocess.run([sys.executable, str(ROOT/'lib/python/workflow/prepare_hard_clip_sequences.py'),
                '--candidates', str(root/'candidates'), '--fastq', str(root/'reads.fastq'),
                '--reads-fasta', str(root/'reads.fasta'), '--reads-index', str(root/'index'),
                '--variants', str(root/'variants'), '--max-size', '4'], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root/'variants').read_text(), 'chr1\t10\t11\ta\tAAGG\tL\t4\nchr1\t20\t21\tb\tTTCC\tR\t4\n')
            self.assertEqual((root/'index').read_text(), '1:b\n3:a\n')
            self.assertEqual((root/'reads.fasta').read_text(), '>b\nTTTTCCCC\n>a\nAAAAGGGG\n')


if __name__ == '__main__':
    unittest.main()
