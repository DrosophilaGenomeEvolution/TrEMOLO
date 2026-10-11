import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'lib/python/workflow'))
from restore_te_names import publish, read_names
from prepare_inputs import prepare
from build_quarto_report import read_te_infos


class RestoreNamesTests(unittest.TestCase):
    def test_exact_names_and_compound_labels_preserve_event_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'internal.tsv', root / 'public.tsv'
            original = 'family|variant: "quoted" & name'
            source.write_text('chrTrEMOLOTE1\t1\t2\tTrEMOLOTE1|event.TrEMOLOTE1\nchr2\t2\t3\tTrEMOLOTE10|event2\n')
            before = source.read_bytes()
            names = {'TrEMOLOTE1': original, 'TrEMOLOTE10': 'TrEMOLOTE1'}
            self.assertEqual(publish(source, output, 'compound:3', names), 2)
            with output.open(newline='') as handle:
                rows = list(csv.reader(handle, delimiter='\t'))
            self.assertEqual(rows[0], ['chrTrEMOLOTE1', '1', '2', original + '|event.TrEMOLOTE1'])
            self.assertEqual(rows[1][3], 'TrEMOLOTE1|event2')
            self.assertEqual(source.read_bytes(), before)
            first = output.read_bytes()
            publish(source, output, 'compound:3', names)
            self.assertEqual(first, output.read_bytes())

    def test_only_family_columns_change_in_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'input', root / 'output'
            source.write_text('reported_te\tcandidate_te\tqseqid\tevent_id\nTrEMOLOTE1\tTrEMOLOTE10\tTrEMOLOTE1:read\tTrEMOLOTE1\n')
            publish(source, output, 'table', {'TrEMOLOTE1': 'A|B', 'TrEMOLOTE10': 'C:D'})
            self.assertEqual(output.read_text(), 'reported_te\tcandidate_te\tqseqid\tevent_id\nA|B\tC:D\tTrEMOLOTE1:read\tTrEMOLOTE1\n')

    def test_documentation_copy_does_not_transform_text(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'input', root / 'output'
            source.write_bytes(b'Documentation\nTrEMOLOTE1\tvalue\n')
            publish(source, output, 'copy', {'TrEMOLOTE1': 'family|name'})
            self.assertEqual(output.read_bytes(), source.read_bytes())

    def test_identity_mapping_preserves_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'input', root / 'output'
            source.write_bytes(b'family\tvalue\r\nroo\t1.000\r\n')
            self.assertEqual(publish(source, output, 'table', {'roo': 'roo'}), 0)
            self.assertEqual(output.read_bytes(), source.read_bytes())

    def test_failed_export_preserves_previous_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'input', root / 'output'
            source.write_text('one\ttwo\n')
            output.write_text('previous\n')
            with self.assertRaisesRegex(ValueError, 'missing TE column'):
                publish(source, output, 'compound:3', {'TrEMOLOTE1': 'roo'})
            self.assertEqual(output.read_text(), 'previous\n')
            self.assertEqual(sorted(p.name for p in root.iterdir()), ['input', 'output'])
            with self.assertRaisesRegex(ValueError, 'must differ'):
                publish(source, source, 'table', {})

    def test_empty_table_and_duplicate_mapping(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, output = root / 'input', root / 'output'
            source.write_text('')
            publish(source, output, 'table', {'TrEMOLOTE1': 'roo'})
            self.assertEqual(output.read_bytes(), b'')
            source.write_text('original\tprepared\na\tTrEMOLOTE1\nb\tTrEMOLOTE1\n')
            with self.assertRaisesRegex(ValueError, 'duplicate'):
                read_names(source)

    def test_prepare_mapping_roundtrip_and_report_with_pipes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome, database = root / 'genome.fa', root / 'te.fa'
            genome.write_text('>chr1\nACGT\n')
            original = 'family|variant: "quoted" & name'
            database.write_text('>' + original + '\nACGT\n>TrEMOLOTE1\nTGCA\n')
            prepare(genome, database, root / 'work')
            names = read_names(root / 'work/INPUT/te_headers.tsv')
            self.assertEqual(names, {'TrEMOLOTE1': original, 'TrEMOLOTE2': 'TrEMOLOTE1'})
            source, output = root / 'input', root / 'output'
            source.write_text('#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\nchr1\t1\t2\tTrEMOLOTE1|out.1\t+\tNONE\t99\t98\t20\t1\t40\t50\t22\tTE_ID_OUTSIDER.1.INS.1\tINS\n')
            publish(source, output, 'compound:3', names)
            calls = read_te_infos(output)
            self.assertEqual(calls[0]['family'], original)
            self.assertEqual(calls[0]['event_id'], 'out.1')

    def test_prepare_normalizes_whitespace_and_rejects_tabs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome, database = root / 'genome.fa', root / 'te.fa'
            genome.write_text('>chr1\nACGT\n')
            database.write_text('>family with spaces\nACGT\n')
            self.assertTrue(prepare(genome, database, root / 'work')['pseudonyms_used'])
            database.write_text('>family\twith tab\nACGT\n')
            with self.assertRaisesRegex(ValueError, 'tab or newline'):
                prepare(genome, database, root / 'work')
