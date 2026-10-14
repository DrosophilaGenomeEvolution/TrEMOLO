#!/usr/bin/env python3
"""Extract genuinely unmapped reads and annotate TE matches in read coordinates."""
import argparse
import csv
import json
import math
import re
import sqlite3
import subprocess
import sys
import tempfile
from pathlib import Path
from contextlib import contextmanager

from annotate_resident_te import (
    read_fragments, read_te_names, build_matches, build_copies,
    interval_union_length, format_number, format_significant,
)

DEFAULTS = dict(MIN_PIDENT=80.0, MIN_ALIGNED_BP=80, MAX_EVALUE=1e-10,
                FULL_LENGTH_COVERAGE=80.0, MIN_CHAIN_GAP=200,
                MAX_GAP_FRACTION=0.2, MAX_OVERLAP=20, AMBIGUITY_OVERLAP=0.8)
READ_COLUMNS = ("query_id", "read_id", "read_length")


def resolve_extraction_parameters(parameters=None):
    parameters = {} if parameters is None else parameters
    defaults = dict(METHOD="two-pass", NAME_MEMORY_MB=256, SAMTOOLS="samtools")
    if not isinstance(parameters, dict) or set(parameters) - set(defaults):
        raise ValueError("Invalid PARAMS.UNMAPPED_TE_EXTRACTION options")
    result = dict(defaults, **parameters)
    if result["METHOD"] not in ("two-pass", "collate"):
        raise ValueError("UNMAPPED_TE_EXTRACTION.METHOD must be two-pass or collate")
    if type(result["NAME_MEMORY_MB"]) is not int or result["NAME_MEMORY_MB"] < 0:
        raise ValueError("UNMAPPED_TE_EXTRACTION.NAME_MEMORY_MB must be a non-negative integer")
    if not isinstance(result["SAMTOOLS"], str) or not result["SAMTOOLS"].strip():
        raise ValueError("UNMAPPED_TE_EXTRACTION.SAMTOOLS must be an executable name or path")
    return result
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


class _CandidateNames:
    """Bounded in-memory index, spilling to a disposable SQLite database."""

    def __init__(self, database, memory_bytes):
        self.database = database
        self.memory_bytes = memory_bytes
        self.names = {}
        self.estimated_bytes = 0

    def add(self, names):
        if self.names is None:
            self.database.executemany("INSERT OR IGNORE INTO reads(name) VALUES (?)", ((n,) for n in names))
            return
        for index, name in enumerate(names):
            if name in self.names:
                continue
            self.names[name] = False
            # Conservative allowance for dictionary slots in addition to strings.
            self.estimated_bytes += sys.getsizeof(name) + 96
            if self.estimated_bytes > self.memory_bytes:
                self.database.execute("CREATE TABLE reads (name TEXT PRIMARY KEY, found INTEGER DEFAULT 0) WITHOUT ROWID")
                self.database.executemany("INSERT INTO reads(name) VALUES (?)", ((n,) for n in self.names))
                self.names = None
                self.database.executemany("INSERT OR IGNORE INTO reads(name) VALUES (?)",
                                          ((n,) for n in names[index + 1:]))
                break

    def count(self, missing=False):
        if self.names is not None:
            return sum(not found for found in self.names.values()) if missing else len(self.names)
        clause = " WHERE found = 0" if missing else ""
        return self.database.execute("SELECT COUNT(*) FROM reads" + clause).fetchone()[0]

    def remove(self, names):
        if self.names is not None:
            for name in names:
                self.names.pop(name, None)
            return
        for batch in _batches(names, 900):
            placeholders = ",".join("?" for _ in batch)
            self.database.execute("DELETE FROM reads WHERE name IN (" + placeholders + ")", batch)

    def lookup(self, names):
        if self.names is not None:
            return {name: self.names[name] for name in names if name in self.names}
        result = {}
        for batch in _batches(names, 900):
            placeholders = ",".join("?" for _ in batch)
            result.update(self.database.execute("SELECT name, found FROM reads WHERE name IN (" + placeholders + ")", batch))
        return result

    def mark_found(self, names):
        if self.names is not None:
            for name in names:
                self.names[name] = True
        else:
            self.database.executemany("UPDATE reads SET found = 1 WHERE name = ?", ((n,) for n in names))


def _batches(items, size=10000):
    batch = []
    for item in items:
        batch.append(item)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def _fastq_batches(source):
    """Bound record count and sequence storage, including very long reads."""
    batch, size = [], 0
    for record in source:
        name, sequence = record.name, record.sequence
        batch.append((name, sequence))
        size += sys.getsizeof(name) + sys.getsizeof(sequence) + 96
        if len(batch) >= 1000 or size >= 4 * 1024 * 1024:
            yield batch
            batch, size = [], 0
    if batch:
        yield batch


@contextmanager
def _collated_alignments(sam, directory, samtools, threads):
    """Stream unfiltered collate output, surfacing errors and reaping the child."""
    import pysam
    command = [samtools, "collate", "-O", "-u", "-@", str(max(0, threads - 2)),
               str(sam), str(Path(directory) / "collate")]
    error_path = Path(directory) / "collate.stderr"
    with error_path.open("wb") as error_file:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=error_file, stdin=subprocess.DEVNULL)
        try:
            try:
                with pysam.AlignmentFile(process.stdout, "rb", check_sq=False) as alignments:
                    yield alignments
            except BaseException as error:
                # A malformed/truncated stream may be the consequence of a
                # samtools error. Capture stderr before scratch files disappear.
                process.stdout.close()
                try:
                    code = process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    process.terminate()
                    try:
                        code = process.wait(timeout=2)
                    except subprocess.TimeoutExpired:
                        process.kill()
                        code = process.wait()
                if code != 0 and isinstance(error, Exception):
                    detail = error_path.read_text(errors="replace").strip()
                    raise RuntimeError(f"samtools collate failed (exit code {code}): {detail}") from error
                raise
            code = process.wait()
            if code != 0:
                detail = error_path.read_text(errors="replace").strip()
                raise RuntimeError(f"samtools collate failed (exit code {code}): {detail}")
        finally:
            process.stdout.close()
            if process.poll() is None:
                process.terminate()
                try:
                    process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.wait()


def _collect_collated(index, alignments):
    """Decide each QNAME once, including every secondary/supplementary record."""
    current, unmapped, mapped = None, False, False
    candidates = 0
    accepted = []
    for row in alignments.fetch(until_eof=True):
        if current is not None and row.query_name != current:
            if unmapped:
                candidates += 1
                if not mapped:
                    accepted.append(current)
                    if len(accepted) == 10000:
                        index.add(accepted)
                        accepted = []
            unmapped, mapped = False, False
        current = row.query_name
        unmapped = unmapped or row.is_unmapped
        mapped = mapped or not row.is_unmapped
    if current is not None and unmapped:
        candidates += 1
        if not mapped:
            accepted.append(current)
    index.add(accepted)
    return candidates


def extract(sam, fastq, fasta, reads, stats, name_memory_mb=256,
            selection_method="two-pass", samtools="samtools", threads=1):
    # Imports stay local: Snakemake can validate configuration without pysam.
    import pysam
    if type(name_memory_mb) is not int or name_memory_mb < 0:
        raise ValueError("name_memory_mb must be a non-negative integer")
    if selection_method not in ("two-pass", "collate"):
        raise ValueError("selection_method must be two-pass or collate")
    if type(threads) is not int or threads < 1:
        raise ValueError("threads must be a positive integer")
    fasta, reads, stats = map(Path, (fasta, reads, stats))
    for path in (fasta, reads, stats):
        path.parent.mkdir(parents=True, exist_ok=True)
    # Two SAM passes preserve arbitrary alignment ordering without storing all
    # mapped names. Only candidate identifiers occupy the bounded memory index.
    with tempfile.TemporaryDirectory(prefix=".unmapped-", dir=str(fasta.parent)) as directory:
        with sqlite3.connect(str(Path(directory) / "names.sqlite")) as db:
            # This database is scratch space, deleted even after an extraction
            # error. Durability/journaling would only add unnecessary disk I/O.
            db.execute("PRAGMA journal_mode=OFF")
            db.execute("PRAGMA synchronous=OFF")
            db.execute("PRAGMA cache_size=-32768")
            index = _CandidateNames(db, name_memory_mb * 1024 * 1024)
            if selection_method == "collate":
                with _collated_alignments(sam, directory, samtools, threads) as alignments:
                    candidates = _collect_collated(index, alignments)
            else:
                with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
                    for batch in _batches(row.query_name for row in alignments.fetch(until_eof=True) if row.is_unmapped):
                        index.add(batch)
                candidates = index.count()
                if candidates:
                    with pysam.AlignmentFile(str(sam), "r", check_sq=False) as alignments:
                        for batch in _batches(row.query_name for row in alignments.fetch(until_eof=True) if not row.is_unmapped):
                            index.remove(batch)
            unmapped = index.count()
            total, selected = 0, 0
            with pysam.FastxFile(str(fastq)) as source, fasta.open("w") as output, reads.open("w", newline="") as table:
                writer = csv.writer(table, delimiter="\t", lineterminator="\n")
                writer.writerow(READ_COLUMNS)

                def emit(name, sequence):
                    nonlocal selected
                    if not sequence:
                        raise ValueError("Empty unmapped read: " + name)
                    selected += 1
                    query_id = "unmapped_{:012d}".format(selected)
                    output.write(">{}\n{}\n".format(query_id, sequence))
                    writer.writerow((query_id, name, len(sequence)))

                if not unmapped:
                    total = sum(1 for _ in source)
                elif index.names is not None:
                    # One dictionary lookup per record. Only selected sequences
                    # are copied from the parser, with no per-record SQL or batches.
                    for record in source:
                        total += 1
                        name = record.name
                        found = index.names.get(name)
                        if found is None:
                            continue
                        if found:
                            raise ValueError("Duplicate unmapped read identifier in FASTQ: " + name)
                        emit(name, record.sequence)
                        index.names[name] = True
                else:
                    for batch in _fastq_batches(source):
                        total += len(batch)
                        found = index.lookup(name for name, _ in batch)
                        emitted = []
                        for name, sequence in batch:
                            if name not in found:
                                continue
                            if found[name]:
                                raise ValueError("Duplicate unmapped read identifier in FASTQ: " + name)
                            emit(name, sequence)
                            found[name] = True
                            emitted.append(name)
                        index.mark_found(emitted)
            missing = index.count(missing=True)
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
    extraction.add_argument("--name-memory-mb", type=int, default=256,
                            help="Candidate-name memory budget in MiB; 0 forces disk-backed batched lookup")
    extraction.add_argument("--selection-method", choices=("two-pass", "collate"), default="two-pass")
    extraction.add_argument("--samtools", default="samtools", help="Samtools executable for collate selection")
    extraction.add_argument("--threads", type=int, default=1,
                            help="CPU budget for collate plus Python; extra compression threads = max(0, N-2)")
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
