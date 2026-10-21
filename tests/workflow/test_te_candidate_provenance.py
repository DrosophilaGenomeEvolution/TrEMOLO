import csv
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from te_candidate_provenance import load_evidence, summarize, read_provenance
from build_te_call_candidates import build
from build_quarto_report import read_ambiguous_call_candidates


class CandidateProvenanceTests(unittest.TestCase):
    def fixture(self, root, queries, hits):
        fasta, provenance, blast = (root / name for name in ("seq.fa", "origin.tsv", "hits.bln"))
        provenance.write_text("query_id\tevent_id\tread_id\tsequence_origin\n" + "".join(
            f"chr1:<INS>:10:10:event.INS.1:2:IMPRECISE:{i}\tevent.INS.1\t{read}\t{origin}\n"
            for i, (sequence, read, origin) in enumerate(queries)))
        fasta.write_text("".join(f">chr1:<INS>:10:10:event.INS.1:2:IMPRECISE:{i}\n{sequence}\n"
                                for i, (sequence, _, _) in enumerate(queries)))
        blast.write_text("".join(f"chr1:<INS>:10:10:event.INS.1:2:IMPRECISE:{i}\t{family}\t99\t{end-start+1}\t0\t0\t{start}\t{end}\t1\t{end-start+1}\t1e-20\t100\n"
                                for i, family, start, end in hits))
        rows = load_evidence([("direct", blast)], [("direct", provenance)], [("direct", fasta)])["event.INS.1"]
        return summarize(rows, {"A", "B"}), (fasta, provenance, blast)

    def test_same_query_competing_annotations_and_distinct_regions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for second_start, expected in ((1, "same_sequence_annotation_ambiguity"), (11, "same_sequence_multiple_regions")):
                (summary, rows), _ = self.fixture(root, [("A"*30, "read:one", "read_sequence")],
                                                [(0,"A",1,10),(0,"B",second_start,second_start+9)])
                self.assertEqual(summary["ambiguity_type"], expected)
                self.assertEqual((summary["sequence_count"], summary["read_count"]), (1,1))
                self.assertEqual(rows[0]["read_id"], "read:one")
                self.assertEqual(rows[0]["query_start"], 0)

    def test_different_reads_can_have_identical_or_different_insertions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for second_sequence in ("A"*30, "C"*30):
                (summary, _), _ = self.fixture(root, [("A"*30,"r1","read_sequence"),(second_sequence,"r2","read_sequence")],
                                               [(0,"A",1,10),(1,"B",1,10)])
                self.assertEqual(summary["read_count"], 2)
                self.assertEqual(summary["query_count"], 2)
                self.assertEqual(summary["sequence_count"], 1 if second_sequence[0]=="A" else 2)
                self.assertIn("different_reads", summary["relation_types"])
                self.assertEqual(summary["ambiguity_type"], "unresolved_family_candidates" if second_sequence[0]=="A" else "different_sequence_evidence")

    def test_mixed_evidence_and_vcf_origin_are_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            (summary, rows), _ = self.fixture(Path(directory), [("A"*30,"r1","read_sequence"),("C"*30,".","vcf_reported")],
                                            [(0,"A",1,10),(0,"B",1,10),(1,"B",1,10)])
            self.assertEqual(summary["ambiguity_type"], "mixed_sequence_evidence")
            self.assertEqual(summary["vcf_query_count"], 1)
            self.assertEqual(summary["read_count"], 1)
            self.assertEqual(rows[-1]["read_id"], ".")

    def test_provenance_cannot_assign_a_vcf_allele_to_a_read(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"origin.tsv"
            path.write_text("query_id\tevent_id\tread_id\tsequence_origin\nq\te\tlast-read\tvcf_reported\n")
            with self.assertRaisesRegex(ValueError,"VCF allele"):
                read_provenance(path)

    def test_builder_and_report_preserve_counts_and_expose_sequence_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, (fasta, origin, blast) = self.fixture(root, [("A"*30,"r1","read_sequence")],[(0,"A",1,10),(0,"B",1,10)])
            infos = root/"TE_INFOS.bed"
            infos.write_text("#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE\nchr1\t10\t11\tA|event.INS.1\t+\tNONE\t99\t100\t30\t10\t50\t50\t30\tTE_ID_OUTSIDER.1\tINS\n")
            counts = root/"counts.tsv"
            counts.write_text("event.INS.1:A\t1\nevent.INS.1:B\t1\n")
            output, evidence = root/"candidates.tsv", root/"evidence.tsv"
            self.assertEqual(build(infos,[("direct",counts)],output,[("direct",blast)],[("direct",origin)],[("direct",fasta)],evidence),(1,2))
            report = read_ambiguous_call_candidates(output,evidence)
            self.assertEqual(len(report["evidence"]),2)
            self.assertEqual(report["candidates"][0]["ambiguity_type"],"same_sequence_annotation_ambiguity")
            self.assertEqual(report["candidates"][0]["evidence_fraction"],0.5)
            self.assertEqual(report["candidates"][0]["provenance_status"],"complete")
            self.assertEqual(report["evidence"][0]["genomic_alignment_status"], "unknown")
            first = output.read_bytes(), evidence.read_bytes()
            build(infos,[("direct",counts)],output,[("direct",blast)],[("direct",origin)],[("direct",fasta)],evidence)
            self.assertEqual(first,(output.read_bytes(),evidence.read_bytes()))
            # Upgrade only the metadata; annotation decisions/counts stay equal.
            origin.write_text(origin.read_text().replace('sequence_origin\n',
                'sequence_origin\tgenomic_alignment_flag\tgenomic_alignment_status\n').replace(
                '\tread_sequence\n', '\tread_sequence\t256\tsecondary\n'))
            build(infos,[("direct",counts)],output,[("direct",blast)],[("direct",origin)],[("direct",fasta)],evidence)
            self.assertEqual(first[0], output.read_bytes())
            updated = read_ambiguous_call_candidates(output,evidence)
            self.assertTrue(all(row['genomic_alignment_status'] == 'secondary' for row in updated['evidence']))
            # Old evidence exports without the two added columns remain readable.
            with evidence.open() as handle:
                rows = list(csv.reader(handle, delimiter='\t'))
            with evidence.open('w', newline='') as handle:
                csv.writer(handle, delimiter='\t').writerows(row[:-2] for row in rows)
            self.assertEqual(len(read_ambiguous_call_candidates(output,evidence)['evidence']), 2)


if __name__ == "__main__":
    unittest.main()
