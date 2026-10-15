import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/locus"))
from build_te_gfa import best_chain, build, reverse_complement, run_input


class GfaTests(unittest.TestCase):
    def test_chromosome_filter_removes_reference_only_contigs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            inputs[0]["reference"] = str(ref)
            manifest = build(ref, db, inputs, root / "out", min_aligned_bp=1,
                             chrom_regex="^chr1$", selected_tes={"C"})
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual({chrom for sample, chrom in walks}, {"chr1"})
            self.assertEqual(manifest["chromosomes"], ["chr1"])
            self.assertEqual(manifest["chrom_regex"], "^chr1$")
            self.assertEqual(manifest["counts"]["representatives"], 3)

    def test_chromosome_filter_excludes_calls_and_candidates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            manifest = build(ref, db, inputs, root / "out", min_aligned_bp=1, chrom_regex="^chr2$")
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual(set(walks), {("reference", "chr2")})
            self.assertEqual(manifest["counts"]["loci"], 0)
            with (root / "out/sequence-alternatives.tsv").open() as handle:
                self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), [])
            with (root / "out/excluded-events.tsv").open() as handle:
                exclusions = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({r["event_id"] for r in exclusions}, {"e1", "e2"})
            self.assertEqual({r["reason"] for r in exclusions}, {"excluded_by_chrom_filter"})

    def test_chromosome_filter_rejects_invalid_or_empty_matches_before_writing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            for pattern, reason in (("[", "Invalid chromosome regex"), ("^missing$", "matches no reference")):
                with self.subTest(pattern=pattern), self.assertRaisesRegex(ValueError, reason):
                    build(ref, db, inputs, root / "out", chrom_regex=pattern)
                self.assertFalse((root / "out").exists())

    def test_cli_chromosome_filter(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            input_file = root / "inputs.json"
            input_file.write_text(json.dumps(inputs))
            result = subprocess.run([sys.executable, str(ROOT / "lib/python/locus/build_te_gfa.py"),
                "--reference", str(ref), "--te-database", str(db), "--inputs", str(input_file),
                "--output", str(root / "out"), "--min-aligned-bp", "1", "--chrom-regex", "^chr1$"],
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual({chrom for sample, chrom in walks}, {"chr1"})

    def fixture(self, root):
        ref, db, infos, seq, origin, blast = [root / n for n in
            ("ref.fa", "db.fa", "TE_INFOS.bed", "ins.fa", "origin.tsv", "hits.bln")]
        ref.write_text(">chr1\nAACCGGTT\n>chr2\nTTAA\n")
        db.write_text(">A\nAAAA\n>B\nCCCC\n>C\nGGGG\n>D\nTTTT\n")
        infos.write_text("\t".join("h" for _ in range(15)) + "\n" + "\n".join(
            f"chr1\t{pos}\t{pos+1}\tA|{event}\t+\tNONE\t99\t100\t4\t{pos}\t50\t50\t4\tTE_ID_OUTSIDER.{i}\tINS"
            for i, (pos, event) in enumerate([(2, "e1"), (6, "e2")])) + "\n")
        queries = [("e1", "AAAA", "r1", "A", "+", 40),
                   ("e1", "AAAT", "r2", "A", "+", 20),
                   ("e1", "AGTC", "r3", "B", "-", 30),
                   ("e1", "GGGG", "r4", "C", "+", 10),
                   ("e2", "CCCC", ".", "B", "+", 25)]
        names = [f"chr1:<INS>:2:2:{event}:1:IMPRECISE:{i}" for i, (event, *_) in enumerate(queries)]
        seq.write_text("".join(f">{q}\n{row[1]}\n" for q, row in zip(names, queries)))
        origin.write_text("query_id\tevent_id\tread_id\tsequence_origin\n" + "".join(
            f"{q}\t{event}\t{read}\t{'vcf_reported' if read == '.' else 'read_sequence'}\n"
            for q, (event, dna, read, *_) in zip(names, queries)))
        blast.write_text("".join(f"{q}\t{family}\t99\t4\t0\t0\t1\t4\t{1 if strand == '+' else 4}\t{4 if strand == '+' else 1}\t1e-5\t{score}\n"
            for q, (event, dna, read, family, strand, score) in zip(names, queries)))
        inputs = [dict(sample_id="pool", te_infos=str(infos), channels=[dict(name="direct",
            alignment=str(blast), provenance=str(origin), insertion_fasta=str(seq))])]
        return ref, db, inputs, names, blast

    def inspect(self, path):
        segments, walks, links = {}, {}, set()
        for line in path.read_text().splitlines():
            row = line.split("\t")
            if row[0] == "S":
                self.assertGreater(int(row[1]), 0)
                self.assertNotIn(row[1], segments)
                segments[row[1]] = row[2]
            elif row[0] == "L":
                self.assertEqual(row[5], "0M")
                links.add(tuple(row[1:5]))
            elif row[0] == "W":
                import re
                tokens = [(n, "+" if o == ">" else "-") for o, n in re.findall(r"([<>])(\d+)", row[6])]
                dna = "".join(segments[n] if o == "+" else reverse_complement(segments[n]) for n, o in tokens)
                self.assertEqual(int(row[5]), len(dna))
                walks[(row[1], row[3])] = (dna, tokens)
        for dna, tokens in walks.values():
            for (a, ao), (b, bo) in zip(tokens, tokens[1:]):
                self.assertIn((a, ao, b, bo), links)
        return segments, walks

    def test_rank_levels_orientation_reference_and_repeatability(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            output = root / "out"
            build(ref, db, inputs, output, min_aligned_bp=1)
            segments, walks = self.inspect(output / "te-graph.gfa")
            self.assertEqual(walks[("reference", "chr1")][0], "AACCGGTT")
            self.assertEqual(walks[("reference", "chr2")][0], "TTAA")
            self.assertEqual(walks[("level_1", "chr1")][0], "AAAAAACCGGCCCC TT".replace(" ", ""))
            self.assertEqual(walks[("level_2", "chr1")][0], "AAAGTCCCGGCCCCTT")
            self.assertEqual(walks[("level_3", "chr1")][0], "AAGGGGCCGGCCCCTT")
            self.assertIn("GACT", segments.values())
            self.assertTrue(any(o == "-" for n, o in walks[("level_2", "chr1")][1]))
            with (output / "sequence-alternatives.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(sum(r["status"] == "representative" for r in rows), 4)
            self.assertEqual(sum(r["status"] == "same_family_alternative" for r in rows), 1)
            self.assertEqual(next(r for r in rows if r["sequence_origin"] == "vcf_reported")["support_reads"], "0")
            before = {p.name: p.read_bytes() for p in output.iterdir()}
            build(ref, db, inputs, output, min_aligned_bp=1)
            self.assertEqual(before, {p.name: p.read_bytes() for p in output.iterdir()})

    def test_same_sequence_multiple_families_are_excluded(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, queries, blast = self.fixture(root)
            with blast.open("a") as handle:
                handle.write(f"{queries[0]}\tD\t99\t4\t0\t0\t1\t4\t1\t4\t1e-5\t100\n")
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            with (root / "out/sequence-alternatives.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({r["status"] for r in rows if r["query_id"] == queries[0]}, {"ambiguous_or_composite_annotation"})
            self.assertFalse(any(r["family"] == "D" and r["status"] == "representative" for r in rows))

    def test_bitscore_chain_is_optimal_collinear_and_nonoverlapping(self):
        def h(qs, qe, ts, te, score, strand="+"):
            return dict(query_start=qs, query_end=qe, te_start=ts, te_end=te, bitscore=score, strand=strand)
        rows = [h(0, 4, 0, 4, 20), h(4, 8, 4, 8, 20), h(0, 8, 0, 8, 35), h(8, 12, 0, 4, 100)]
        self.assertEqual(sum(r["bitscore"] for r in best_chain(rows[:3])), 40)
        self.assertEqual(sum(r["bitscore"] for r in best_chain(rows)), 100)
        self.assertEqual(sum(r["bitscore"] for r in best_chain([h(0, 4, 4, 8, 20, "-"), h(4, 8, 0, 4, 20, "-")])), 40)

    def test_invalid_reference_anchor_and_thresholds(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            infos = Path(inputs[0]["te_infos"])
            infos.write_text(infos.read_text().replace("\t2\t50", "\t99\t50"))
            with self.assertRaisesRegex(ValueError, "Invalid reference anchor"):
                build(ref, db, inputs, root / "out", min_aligned_bp=1)
            with self.assertRaisesRegex(ValueError, "thresholds"):
                build(ref, db, inputs, root / "out", min_identity=float("nan"))

    def test_terminal_insertions_have_no_empty_nodes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            infos = Path(inputs[0]["te_infos"])
            infos.write_text(infos.read_text().replace("\t2\t50", "\t0\t50").replace("\t6\t50", "\t8\t50"))
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            segments, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertTrue(all(segments.values()))
            self.assertEqual(walks[("level_1", "chr1")][0], "AAAAAACCGGTTCCCC")

    def test_run_adapter_reports_old_provenance_and_reference_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            directory = root / "OUTSIDER/TrEMOLO_SV_TE/INS"
            directory.mkdir(parents=True)
            (directory / "SV_INS_CLUST.insertions.fasta").write_text(">q\nAAAA\n")
            with self.assertRaisesRegex(ValueError, "Incomplete channel"):
                run_input("pool=" + str(root))
            different = root / "different.fa"
            different.write_text(">chr1\nTTTT\n")
            inputs[0]["reference"] = str(different)
            with self.assertRaisesRegex(ValueError, "different reference"):
                build(ref, db, inputs, root / "out", min_aligned_bp=1)

    def test_higher_bitscore_wins_even_with_lower_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, queries, blast = self.fixture(root)
            blast.write_text(blast.read_text().replace("\t99\t4\t", "\t90\t4\t", 1))
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            with (root / "out/sequence-alternatives.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(next(r for r in rows if r["query_id"] == queries[0])["status"], "representative")

    def test_original_te_headers_with_spaces_are_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, blast = self.fixture(root)
            db.write_text(db.read_text().replace(">A\n", ">A original label\n"))
            blast.write_text(blast.read_text().replace("\tA\t", "\tTrEMOLOTE1\t"))
            mapping = root / "names.tsv"
            mapping.write_text("original\tprepared\nA original label\tTrEMOLOTE1\n")
            inputs[0]["te_headers"] = str(mapping)
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            with (root / "out/sequence-alternatives.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertIn("A original label", {r["family"] for r in rows})

    def test_empty_channel_produces_reference_and_explicit_exclusions(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, blast = self.fixture(root)
            blast.write_text("")
            Path(inputs[0]["channels"][0]["insertion_fasta"]).write_text("")
            build(ref, db, inputs, root / "out", min_aligned_bp=1)
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual(set(walks), {("reference", "chr1"), ("reference", "chr2")})
            with (root / "out/excluded-events.tsv").open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 2)

    def test_te_filter_keeps_entire_locus_when_only_an_alternative_matches(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            manifest = build(ref, db, inputs, root / "out", min_aligned_bp=1, selected_tes={"C"})
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual(walks[("level_1", "chr1")][0], "AAAAAACCGGTT")
            self.assertEqual(walks[("level_2", "chr1")][0], "AAAGTCCCGGTT")
            self.assertEqual(walks[("level_3", "chr1")][0], "AAGGGGCCGGTT")
            self.assertEqual(manifest["counts"]["selected_loci"], 1)
            self.assertEqual(manifest["counts"]["representatives"], 3)
            with (root / "out/sequence-alternatives.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            excluded = [r for r in rows if r["event_id"] == "e2"]
            self.assertEqual({r["status"] for r in excluded}, {"excluded_by_te_filter"})
            self.assertEqual({r["node_id"] for r in excluded}, {"."})

    def test_te_filter_matches_primary_label_without_representative_of_that_family(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            infos = Path(inputs[0]["te_infos"])
            infos.write_text(infos.read_text().replace("A|e1", "D|e1"))
            manifest = build(ref, db, inputs, root / "out", min_aligned_bp=1, selected_tes={"D"})
            self.assertEqual(manifest["counts"]["representatives"], 3)
            self.assertEqual(manifest["selected_tes"], ["D"])

    def test_te_filter_no_matches_leaves_reference_and_empty_selection_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            manifest = build(ref, db, inputs, root / "out", min_aligned_bp=1, selected_tes={"D"})
            self.assertEqual(manifest["counts"]["representatives"], 0)
            _, walks = self.inspect(root / "out/te-graph.gfa")
            self.assertEqual(set(walks), {("reference", "chr1"), ("reference", "chr2")})
            with self.assertRaisesRegex(ValueError, "at least one name"):
                build(ref, db, inputs, root / "out", selected_tes=set())

    def test_cli_accepts_file_or_comma_list_with_identical_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            ref, db, inputs, _, _ = self.fixture(root)
            input_file = root / "inputs.json"
            input_file.write_text(json.dumps(inputs))
            names_file = root / "TE_names.txt"
            names_file.write_text("C\nD\n\n")
            for suffix, option in (("csv", ["--te-names", "C,D"]),
                                   ("file", ["--te-names-file", str(names_file)])):
                result = subprocess.run([sys.executable, str(ROOT / "lib/python/locus/build_te_gfa.py"),
                    "--reference", str(ref), "--te-database", str(db), "--inputs", str(input_file),
                    "--output", str(root / suffix), "--min-aligned-bp", "1", *option], capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((root / "csv/te-graph.gfa").read_bytes(), (root / "file/te-graph.gfa").read_bytes())


if __name__ == "__main__":
    unittest.main()
