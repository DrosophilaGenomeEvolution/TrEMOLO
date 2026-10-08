#!/usr/bin/env python3
"""Validate corrected population genomes and CIGAR sequence provenance."""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "lib/python/workflow"))
from integrate_outsider_te import read_fasta, reverse_complement  # noqa: E402
from summarize_outsider_liftoff import read_lifted_features  # noqa: E402


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def bed(path):
    return [(f[0], int(f[1]), int(f[2]), f[3]) for f in
            (line.split("\t") for line in path.read_text().splitlines())]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    args = parser.parse_args()
    work = args.workdir
    outsider = work / "OUTSIDER"
    integration = outsider / "TE_TOWARD_GENOME"
    audit = table(integration / "INTEGRATION_TE.tsv")
    merged = bed(outsider / "TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed")
    expected = {name.rsplit("|", 1)[1] for _, _, _, name in merged
                if name.rsplit("|", 1)[1].split(".", 2)[1:2] == ["INS"]}
    integrated = {r["event_id"]: r for r in audit if r["status"] == "integrated"}
    assert set(integrated) == expected, "Not all retained INS were integrated"
    assert not any(r["status"] == "rejected" for r in audit), "Fixture integration is partial"
    assert all(identifier.split(".", 2)[1:2] == ["INS"] for identifier in integrated), "Non-INS sequence inserted"
    assert all(r["reason"] == "non_insertion_event" for r in audit if r["status"] == "excluded")

    genome = read_fasta(work / "INPUT/genome.fasta")
    database = read_fasta(work / "INPUT/te_database.fasta")
    direct = outsider / "TrEMOLO_SV_TE/INS"
    sniffles = outsider / "VARIANT_CALLING"
    sequences = {"direct": read_fasta(direct / "SV_INS_CLUST.insertions.fasta"),
                 "sniffles": read_fasta(sniffles / "SEQUENCE_INDEL.insertions.fasta")}
    fragments = {"direct": read_fasta(direct / "SV_INS_CLUST.fasta"),
                 "sniffles": read_fasta(sniffles / "SEQUENCE_INDEL.fasta")}
    read_names = {}
    for line in (direct / "RD_NUMBER.txt").read_text().splitlines():
        fields = line.split(":")
        read_names[(fields[4], fields[5])] = fields[6]
    for line in (outsider / "TE_DETECTION/RD_NUMBER.txt").read_text().splitlines():
        fields = line.split(":")
        read_names[(fields[0], fields[2])] = fields[1]

    needed = {}
    observed = {}
    trimmed = 0
    for identifier, row in integrated.items():
        key = row["qseqid"][:-2]
        sequence = sequences[row["source"]][key]
        assert len(sequence) == int(row["observed_length"]) and sequence
        observed[identifier] = sequence
        # Index zero denotes the caller consensus, with no BAM flank to trim.
        number = key.split(":")[-1]
        if row["source"] == "direct" or number != "0":
            read_name = read_names[(identifier, number)]
            needed[identifier] = (read_name, row["chromosome"], sequence)
            assert sequence in fragments[row["source"]][key], "Insertion absent from classified fragment"
            trimmed += len(fragments[row["source"]][key]) - len(sequence)

    provenance = defaultdict(set)
    bam = outsider / "MAPPING/SAMPLE_mapping_GENOME_MD.sorted.bam"
    selected_names = {value[0] for value in needed.values()}
    with pysam.AlignmentFile(str(bam), "rb") as handle:
        for read in handle.fetch(until_eof=True):
            if read.query_name not in selected_names or not read.query_sequence or read.is_unmapped:
                continue
            query = 0
            for operation, length in read.cigartuples:
                if operation == 1:
                    provenance[(read.query_name, read.reference_name)].add(read.query_sequence[query:query + length].upper())
                if operation in (0, 1, 4, 7, 8):
                    query += length
    for identifier, (name, chromosome, sequence) in needed.items():
        assert sequence in provenance[(name, chromosome)], "Not an exact CIGAR insertion: " + identifier

    for fasta_name, bed_name, public_bed, sequence_key in [
        ("NEO_GENOME.fasta", "TRUE_POSITION_TE_NEO.bed", "POSITION_TE_OUTSIDER_IN_NEO_GENOME.bed", "observed"),
        ("PSEUDO_GENOME_TE_DB_ID.fasta", "TRUE_POSITION_TE_PSEUDO.bed", "POSITION_TE_OUTSIDER_IN_PSEUDO_GENOME.bed", "canonical"),
    ]:
        rebuilt = read_fasta(integration / fasta_name)
        positions = bed(integration / bed_name)
        assert set(rebuilt) == set(genome) and len(positions) == len(expected)
        assert bed(work / public_bed) == [(c, s, e, n.replace(":", "|")) for c, s, e, n in positions]
        by_chrom = defaultdict(list)
        for chromosome, start, end, name in positions:
            identifier = name.rsplit(":", 1)[1]
            row = integrated[identifier]
            sequence = observed[identifier] if sequence_key == "observed" else database[row["te_family"]].upper()
            if sequence_key == "canonical" and row["strand"] == "-":
                sequence = reverse_complement(sequence)
            assert 0 <= start < end <= len(rebuilt[chromosome])
            assert rebuilt[chromosome][start:end] == sequence, "BED does not delimit the inserted sequence"
            assert end - start == int(row[sequence_key + "_length"])
            by_chrom[chromosome].append((start, end, row))
        for chromosome, original in genome.items():
            recovered = []
            cursor = shift = 0
            for start, end, row in sorted(by_chrom[chromosome], key=lambda item: item[0]):
                assert start >= cursor, "Integrated BED intervals overlap"
                assert start - shift == int(row["breakpoint"]), "Incorrect cumulative coordinate shift"
                recovered.append(rebuilt[chromosome][cursor:start])
                cursor = end
                shift += end - start
            recovered.append(rebuilt[chromosome][cursor:])
            assert "".join(recovered) == original, "Integration altered or duplicated original genome bases"

    lift_dir = outsider / "INSIDER_VR"
    if (lift_dir / "INOUTSIDER.gff").is_file():
        canonical_positions = {name.rsplit(":", 1)[1]: (c, s, e) for c, s, e, name in bed(integration / "TRUE_POSITION_TE_PSEUDO.bed")}
        for identifier, features in read_lifted_features(lift_dir / "INOUTSIDER.gff").items():
            chromosome, start, end = canonical_positions[identifier]
            for feature in features:
                assert feature["chromosome"] == chromosome and 1 <= feature["start"] <= feature["end"]
                if feature["side"] == "L":
                    assert feature["end"] == start
                else:
                    assert feature["start"] == end + 1
        mapped = read_lifted_features(lift_dir / "output_INOUT.gff3")
        public = {name.rsplit("|", 1)[1]: (c, s, e) for c, s, e, name in bed(work / "POS_TE_OUTSIDER_ON_REF.bed")}
        for identifier, (chromosome, start, end) in public.items():
            pair = {feature["side"]: feature for feature in mapped[identifier]}
            assert len(mapped[identifier]) == 2 and set(pair) == {"L", "R"}
            left, right = pair["L"], pair["R"]
            assert left["chromosome"] == right["chromosome"] == chromosome
            assert left["strand"] == right["strand"] in ("+", "-")
            boundaries = (left["end"], right["start"] - 1) if left["strand"] == "+" else (right["end"], left["start"] - 1)
            assert (start, end) == (min(boundaries), max(boundaries)), "Incorrect projected BED conversion"
    print("Population integration: PASSED ({} INS, {} excluded, {} exact CIGAR sequences, {} flank bases removed)".format(
        len(integrated), sum(r["status"] == "excluded" for r in audit), len(needed), trimmed))


if __name__ == "__main__":
    main()
