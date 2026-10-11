#!/usr/bin/env python3
"""Extract genuinely unmapped reads and annotate TE matches in read coordinates."""
import argparse
import csv
import json
import math
import re
import sqlite3
import tempfile
from pathlib import Path

from annotate_resident_te import (
    read_fragments, read_te_names, build_matches, build_copies,
    interval_union_length, format_number, format_significant,
)

DEFAULTS = dict(MIN_PIDENT=80.0, MIN_ALIGNED_BP=80, MAX_EVALUE=1e-10,
                FULL_LENGTH_COVERAGE=80.0, MIN_CHAIN_GAP=200,
                MAX_GAP_FRACTION=0.2, MAX_OVERLAP=20, AMBIGUITY_OVERLAP=0.8)
READ_COLUMNS = ("query_id", "read_id", "read_length")
MATCH_COLUMNS = ("match_id", "component_id", "read_id", "read_length", "te_family",
                 "read_start", "read_end", "strand", "te_length", "te_covered_bp",
                 "te_coverage", "read_aligned_bp", "identity", "bitscore",
                 "best_evalue", "fragment_count", "coverage_status", "assignment",
                 "candidate_count", "read_segments")
HSP_COLUMNS = ("fragment_id", "match_id", "read_id", "read_length", "te_family",
               "read_start", "read_end", "strand", "te_start", "te_end", "te_length",
               "aligned_bp", "identity", "evalue", "bitscore", "filter_status",
               "filter_reason", "redundant_with")


def resolve_parameters(parameters=None):
    parameters = {} if parameters is None else parameters
    if not isinstance(parameters, dict):
        raise ValueError("PARAMS.UNMAPPED_TE must be a mapping")
    if set(parameters) - set(DEFAULTS):
        raise ValueError("Unknown PARAMS.UNMAPPED_TE: " + ", ".join(sorted(set(parameters) - set(DEFAULTS))))
    result = dict(DEFAULTS, **parameters)
    for name, value in result.items():
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError("PARAMS.UNMAPPED_TE.{} must be finite and numeric".format(name))
        if name in ("MIN_ALIGNED_BP", "MIN_CHAIN_GAP", "MAX_OVERLAP") and type(value) is not int:
            raise ValueError("PARAMS.UNMAPPED_TE.{} must be an integer".format(name))
        if value < 0 or (name == "MIN_ALIGNED_BP" and value == 0):
            raise ValueError("Invalid PARAMS.UNMAPPED_TE." + name)
        if name in ("MIN_PIDENT", "FULL_LENGTH_COVERAGE") and value > 100:
            raise ValueError("PARAMS.UNMAPPED_TE.{} must be <= 100".format(name))
        if name == "AMBIGUITY_OVERLAP" and not 0 < value <= 1:
            raise ValueError("PARAMS.UNMAPPED_TE.AMBIGUITY_OVERLAP must be in (0, 1]")
    return result


def write_table(path, header, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def extract(sam, fastq, fasta, reads, stats):
    # Imports stay local: Snakemake can validate configuration without pysam.
    import pysam
    fasta, reads, stats = map(Path, (fasta, reads, stats))
    for path in (fasta, reads, stats):
        path.parent.mkdir(parents=True, exist_ok=True)
    # Disk-backed names keep memory bounded on large populations. Two SAM passes
    # avoid storing all mapped names, which usually dominate the dataset.
    with tempfile.TemporaryDirectory(prefix=".unmapped-", dir=str(fasta.parent)) as directory:
        with sqlite3.connect(str(Path(directory) / "names.sqlite")) as db:
            db.execute("CREATE TABLE reads (name TEXT PRIMARY KEY, found INTEGER DEFAULT 0)")
            with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
                db.executemany("INSERT OR IGNORE INTO reads(name) VALUES (?)",
                               ((row.query_name,) for row in alignments.fetch(until_eof=True) if row.is_unmapped))
            candidates = db.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
            with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
                db.executemany("DELETE FROM reads WHERE name = ?",
                               ((row.query_name,) for row in alignments.fetch(until_eof=True) if not row.is_unmapped))
            unmapped = db.execute("SELECT COUNT(*) FROM reads").fetchone()[0]
            total, selected = 0, 0
            with pysam.FastxFile(str(fastq)) as source, fasta.open("w") as output, reads.open("w", newline="") as table:
                writer = csv.writer(table, delimiter="\t", lineterminator="\n")
                writer.writerow(READ_COLUMNS)
                for record in source:
                    total += 1
                    candidate = db.execute("SELECT found FROM reads WHERE name = ?", (record.name,)).fetchone()
                    if candidate is None:
                        continue
                    if candidate[0]:
                        raise ValueError("Duplicate unmapped read identifier in FASTQ: " + record.name)
                    if not record.sequence:
                        raise ValueError("Empty unmapped read: " + record.name)
                    selected += 1
                    query_id = "unmapped_{:012d}".format(selected)
                    output.write(">{}\n{}\n".format(query_id, record.sequence))
                    writer.writerow((query_id, record.name, len(record.sequence)))
                    db.execute("UPDATE reads SET found = 1 WHERE name = ?", (record.name,))
            missing = db.execute("SELECT COUNT(*) FROM reads WHERE found = 0").fetchone()[0]
            if missing:
                raise ValueError("{} unmapped SAM identifiers absent from original reads".format(missing))
    summary = dict(total_reads=total, unmapped_candidates=candidates,
                   excluded_with_mapped_alignment=candidates - unmapped, unmapped_reads=selected)
    stats.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def segments(fragments):
    intervals = sorted((row.start, row.end) for row in fragments)
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return ";".join("{}-{}".format(*row) for row in merged)


def annotate(blast, headers, reads, extraction_stats, output, hsps, stats, parameters=None):
    params = resolve_parameters(parameters)
    with Path(reads).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != READ_COLUMNS:
            raise ValueError("Invalid unmapped read table header")
        names = {row["query_id"]: (row["read_id"], int(row["read_length"])) for row in reader}
    fragments = read_fragments(Path(blast), read_te_names(Path(headers)), "UNMAPPED_TE", re.compile("."),
                              params["MIN_PIDENT"], params["MIN_ALIGNED_BP"], params["MAX_EVALUE"])
    for row in fragments:
        if row.chrom not in names or names[row.chrom][1] != row.target_length:
            raise ValueError("Unknown read or inconsistent BLAST read length: " + row.chrom)
    matches = build_matches(fragments, "UNMAPPED_TE", params["MIN_CHAIN_GAP"],
                            params["MAX_GAP_FRACTION"], params["MAX_OVERLAP"],
                            params["FULL_LENGTH_COVERAGE"], params["MIN_PIDENT"], 0)
    copies = build_copies(matches, "UNMAPPED_TE", params["AMBIGUITY_OVERLAP"])
    counts = {row.copy_id: len(row.matches) for row in copies}
    write_table(output, MATCH_COLUMNS, (
        (row.match_id, row.copy_id, *names[row.chrom], row.te_name, row.start, row.end,
         row.strand, row.consensus_length, row.consensus_covered_bp,
         format_number(row.consensus_coverage), interval_union_length((f.start, f.end) for f in row.fragments),
         format_number(row.identity), format_number(row.bitscore), format_significant(row.best_evalue),
         len(row.fragments), "near_full_length" if row.tier == "full_length" else "partial",
         row.assignment, counts[row.copy_id], segments(row.fragments)) for row in matches))
    write_table(hsps, HSP_COLUMNS, (
        (row.fragment_id, row.match_id, *names[row.chrom], row.te_name, row.start, row.end,
         row.strand, row.consensus_start, row.consensus_end, row.consensus_length,
         row.aligned_bp, format_number(row.pident), format_significant(row.evalue),
         format_number(row.bitscore), row.filter_status, row.filter_reason, row.redundant_with)
        for row in fragments))
    summary = json.loads(Path(extraction_stats).read_text())
    summary.update(reads_with_te=len({row.chrom for row in matches}), te_matches=len(matches),
                   components=len(copies), ambiguous_components=sum(len(row.matches) > 1 for row in copies),
                   blast_fragments=len(fragments), parameters=params,
                   coordinate_system="0-based, half-open; original read sequence; no genomic position")
    Path(stats).write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    extraction = sub.add_parser("extract")
    for name in ("sam", "fastq", "fasta", "reads", "stats"):
        extraction.add_argument("--" + name, type=Path, required=True)
    annotation = sub.add_parser("annotate")
    for name in ("blast", "headers", "reads", "extraction-stats", "output", "hsps", "stats"):
        annotation.add_argument("--" + name, type=Path, required=True)
    annotation.add_argument("--parameters", default="{}", help="JSON parameter mapping")
    args = vars(parser.parse_args())
    command = args.pop("command")
    if command == "extract":
        extract(**args)
    else:
        args["parameters"] = json.loads(args["parameters"])
        annotate(**args)


if __name__ == "__main__":
    main()
