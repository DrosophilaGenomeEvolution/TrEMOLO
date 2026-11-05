import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'lib/python/workflow'))
from report_frequency_depth import attach_frequency_depth
from summarize_outsider_frequency import FREQUENCY_COLUMNS


class FrequencyDepthTests(unittest.TestCase):
    def call(self, event='e', source='OUTSIDER', frequency=80, clipped=80, family='TE'):
        return dict(event_id=event, source=source, chrom='chr1', family=family,
                    display_frequency=frequency, frequency_with_clipped=clipped)

    def table(self, root, name, rows):
        path = root / name
        path.write_text('\t'.join(FREQUENCY_COLUMNS) + '\n' + ''.join('\t'.join(map(str, row)) + '\n' for row in rows))
        return path

    def test_builder_attaches_depth_without_changing_calls_or_frequency(self):
        from tests.workflow.test_build_quarto_report import BuildQuartoReportTests
        from build_quarto_report import build_report_data
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            infos, index, manifest, mapping, vcf = BuildQuartoReportTests().fixture(root)
            original = infos.read_bytes()
            raw = self.table(root, 'raw', [['chr1:100-101','out.1','roo',4,1,10,5,10,40,50,'INS']])
            args = dict(te_infos=infos, genome_index=index, manifest_path=manifest,
                        mapping_stats_path=mapping, sv_vcf_path=vcf,
                        enabled_sources=['OUTSIDER','INSIDER'], title='depth', author='',
                        work_directory=str(root), locus_window=100)
            before = build_report_data(**args)
            after = build_report_data(**args, outsider_frequency_path=raw)
            self.assertEqual(after['calls'][0]['frequency_depth'], 10)
            self.assertEqual(after['calls'][0]['frequency_support_reads'], 5)
            self.assertEqual(before['summary'], after['summary'])
            for left, right in zip(before['calls'], after['calls']):
                self.assertEqual({k:v for k,v in left.items() if not k.startswith('frequency_depth') and k != 'frequency_support_reads'},
                                 {k:v for k,v in right.items() if not k.startswith('frequency_depth') and k != 'frequency_support_reads'})
            self.assertEqual(infos.read_bytes(), original)

    def test_precise_selection_raw_fallback_and_te_presence_for_deletion(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.table(root, 'raw', [
                ['chr1:10-11','e','TE',7,1,10,8,10,70,80,'INS'],
                ['chr1:20-21','d','TE',2,'.',10,8,10,80,80,'DEL']])
            precise = self.table(root, 'precise', [
                ['chr1:10-11','e','TE',70,10,100,80,100,70,80,'INS'],
                ['chr1:20-21','d','TE',12,'.',10,12,10,120,120,'DEL']])
            calls = [self.call(), self.call('d')]
            attach_frequency_depth(calls, raw, precise)
            self.assertEqual((calls[0]['frequency_depth'], calls[0]['frequency_support_reads'], calls[0]['frequency_depth_table']), (100,80,'precise'))
            self.assertEqual((calls[1]['frequency_depth'], calls[1]['frequency_support_reads'], calls[1]['frequency_depth_table']), (10,8,'raw'))
            self.assertEqual([c['frequency_depth_status'] for c in calls], ['matched','matched'])

    def test_non_clipped_depth_and_no_matching_estimate_are_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            raw = self.table(root, 'raw', [['chr1:10-11','e','TE',4,4,5,8,10,80,80,'INS']])
            calls = [self.call(clipped=None), self.call(frequency=50), self.call('missing'), self.call(frequency=None)]
            attach_frequency_depth(calls, raw)
            self.assertEqual(calls[0]['frequency_depth'], 5)
            self.assertEqual(calls[0]['frequency_support_reads'], 4)
            self.assertEqual([c['frequency_depth_status'] for c in calls], ['matched','frequency_mismatch','unavailable','frequency_unavailable'])
            self.assertIsNone(calls[1]['frequency_depth'])

    def test_insider_names_zero_depth_and_negative_support(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            headers = root/'headers';headers.write_text('original\tprepared\nTE| original label\tte_1\n')
            insider = root/'insider';insider.write_text('chrom\tposition\ttotal_depth\tdepth_empty_site\tread_support\tread_support_percent\tTE\tinfo_TE\n'
                'chr1\t10\t20\t4\t16\t80.0000\tte_1\tte_1|e\n'
                'chr1\t11\t0\t0\t0\t0.0000\tte_1\tte_1|zero\n'
                'chr1\t12\t2\t3\t-1\t-50.0000\tte_1\tte_1|bad\n')
            calls=[self.call(source='INSIDER',clipped=None,family='TE| original label'),
                   self.call('zero','INSIDER',0,None,'TE| original label'),
                   self.call('bad','INSIDER',-50,None,'TE| original label')]
            attach_frequency_depth(calls, insider=insider, te_headers=headers)
            self.assertEqual(calls[0]['frequency_depth'],20)
            self.assertEqual(calls[0]['frequency_depth_method'],'mean_flanking_depth')
            self.assertEqual([c['frequency_depth_status'] for c in calls],['matched','zero_depth','invalid_counts'])

    def test_mismatched_chromosome_and_invalid_counts_are_not_plotted(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            raw=self.table(root,'raw',[
                ['chr2:10-11','e','TE',8,0,10,8,10,80,80,'INS'],
                ['chr1:10-11','bad','TE',8,0,-10,8,-10,80,80,'INS']])
            calls=[self.call(),self.call('bad')]
            attach_frequency_depth(calls,raw)
            self.assertEqual([c['frequency_depth_status'] for c in calls],['chromosome_mismatch','invalid_counts'])
