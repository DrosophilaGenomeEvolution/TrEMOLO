#!/usr/bin/env python3
"""Normalize insertion-to-TE alignments and build a standalone structure report."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Optional


PIPELINE_ROOT = Path(__file__).resolve().parents[3]
SCHEMA_VERSION = "1.0.0"
HSP_COLUMNS = (
    "hsp_id", "query_id", "event_id", "chrom", "locus_start", "locus_end",
    "subject_te", "pident", "alignment_length", "mismatch", "gapopen",
    "query_start", "query_end", "subject_start", "subject_end", "subject_strand",
    "evalue", "bitscore", "query_length", "subject_length", "retained",
    "filter_reasons",
)
MATCH_COLUMNS = (
    "match_id", "query_id", "event_id", "chrom", "locus_start", "locus_end",
    "subject_te", "subject_strand", "query_start", "query_end", "query_segments",
    "aligned_query_bp", "aligned_subject_bp", "query_coverage", "consensus_coverage",
    "weighted_identity", "bitscore_sum", "best_evalue", "hsp_count", "query_length",
    "subject_length", "status", "filter_reasons", "final_call", "reported_te",
    "assignment",
)
COMPONENT_COLUMNS = (
    "component_id", "structure_id", "query_id", "event_id", "rank", "match_id",
    "te_name", "query_start", "query_end", "query_segments", "strand",
    "consensus_coverage", "weighted_identity", "assignment", "interpretation",
)
STRUCTURE_COLUMNS = (
    "structure_id", "query_id", "event_id", "chrom", "locus_start", "locus_end",
    "query_length", "final_call", "reported_te", "retained_matches",
    "component_count", "alternative_matches", "classification", "interpretation",
)
EVENT_COLUMNS = (
    "event_structure_id", "event_id", "chrom", "locus_start", "locus_end",
    "final_call", "reported_te", "query_count", "queries_with_matches",
    "queries_with_components", "multi_component_queries", "multi_family_queries",
    "maximum_component_count", "classification", "interpretation",
)


def stable_id(prefix: str, *parts: object) -> str:
    value = "\x1f".join(map(str, parts)).encode("utf-8")
    return prefix + hashlib.sha256(value).hexdigest()[:16]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent), text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(value)
        os.replace(str(temporary), str(path))
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def write_tsv(path: Path, columns: tuple[str, ...], rows: Iterable[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent), text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        os.replace(str(temporary), str(path))
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def read_sizes(path: Path, label: str) -> dict[str, int]:
    sizes = {}
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected at least two columns")
            try:
                size = int(fields[1])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid {label} size") from error
            if size < 1:
                raise ValueError(f"{path}:{line_number}: {label} size must be positive")
            if fields[0] in sizes and sizes[fields[0]] != size:
                raise ValueError(f"{path}:{line_number}: conflicting size for {fields[0]}")
            sizes[fields[0]] = size
    return sizes


def parse_query_id(query_id: str) -> tuple[str, int, int, str]:
    fields = query_id.split(":")
    if len(fields) < 5:
        raise ValueError(f"malformed insertion query id: {query_id}")
    try:
        start = int(fields[2])
        end = int(fields[3])
    except ValueError as error:
        raise ValueError(f"malformed insertion coordinates: {query_id}") from error
    return fields[0], start, end, fields[4]


def read_hsps(
    blast: Path,
    query_sizes: dict[str, int],
    subject_sizes: dict[str, int],
    min_pident: float,
    min_aligned_bp: int,
) -> list[dict]:
    rows = []
    with blast.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 12:
                raise ValueError(f"{blast}:{line_number}: expected BLAST outfmt 6")
            query_id, subject_te = fields[:2]
            if query_id not in query_sizes:
                raise ValueError(f"{blast}:{line_number}: missing query size for {query_id}")
            if subject_te not in subject_sizes:
                raise ValueError(f"{blast}:{line_number}: missing TE size for {subject_te}")
            try:
                pident = float(fields[2])
                alignment_length = int(fields[3])
                mismatch = int(fields[4])
                gapopen = int(fields[5])
                raw_qstart, raw_qend = int(fields[6]), int(fields[7])
                raw_sstart, raw_send = int(fields[8]), int(fields[9])
                evalue = float(fields[10])
                bitscore = float(fields[11])
            except ValueError as error:
                raise ValueError(f"{blast}:{line_number}: invalid numeric BLAST value") from error
            if (
                not math.isfinite(pident)
                or not math.isfinite(evalue)
                or not math.isfinite(bitscore)
                or not 0 <= pident <= 100
                or evalue < 0
                or alignment_length < 1
            ):
                raise ValueError(f"{blast}:{line_number}: invalid alignment values")
            chrom, locus_start, locus_end, event_id = parse_query_id(query_id)
            qstart, qend = min(raw_qstart, raw_qend) - 1, max(raw_qstart, raw_qend)
            sstart, send = min(raw_sstart, raw_send) - 1, max(raw_sstart, raw_send)
            if (
                qstart < 0
                or qend > query_sizes[query_id]
                or sstart < 0
                or send > subject_sizes[subject_te]
            ):
                raise ValueError(f"{blast}:{line_number}: alignment exceeds sequence size")
            strand = "+" if raw_sstart <= raw_send else "-"
            reasons = []
            if pident < min_pident:
                reasons.append("LOW_IDENTITY")
            if alignment_length < min_aligned_bp:
                reasons.append("SHORT_HSP")
            rows.append(
                {
                    "hsp_id": stable_id(
                        "H", query_id, subject_te, line_number, qstart, qend, sstart, send
                    ),
                    "query_id": query_id,
                    "event_id": event_id,
                    "chrom": chrom,
                    "locus_start": locus_start,
                    "locus_end": locus_end,
                    "subject_te": subject_te,
                    "pident": pident,
                    "alignment_length": alignment_length,
                    "mismatch": mismatch,
                    "gapopen": gapopen,
                    "query_start": qstart,
                    "query_end": qend,
                    "subject_start": sstart,
                    "subject_end": send,
                    "subject_strand": strand,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "query_length": query_sizes[query_id],
                    "subject_length": subject_sizes[subject_te],
                    "retained": "yes" if not reasons else "no",
                    "filter_reasons": ";".join(reasons) if reasons else ".",
                }
            )
    return rows


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def interval_length(intervals: Iterable[tuple[int, int]]) -> int:
    return sum(end - start for start, end in merge_intervals(intervals))


def serialize_intervals(intervals: Iterable[tuple[int, int]]) -> str:
    return ";".join(f"{start}-{end}" for start, end in merge_intervals(intervals))


def read_final_calls(path: Optional[Path]) -> dict[str, str]:
    if path is None:
        return {}
    calls = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if len(header) != 15:
            raise ValueError(f"{path}: expected 15-column TE_INFOS.bed")
        for line_number, row in enumerate(reader, 2):
            if not row:
                continue
            if len(row) != 15:
                raise ValueError(f"{path}:{line_number}: expected 15 columns")
            family, separator, event_id = row[3].partition("|")
            if not separator:
                raise ValueError(f"{path}:{line_number}: malformed TE|ID")
            previous = calls.get(event_id)
            if previous is not None and previous != family:
                raise ValueError(f"{path}:{line_number}: conflicting final call for {event_id}")
            calls[event_id] = family
    return calls


def read_candidate_calls(path: Optional[Path]) -> dict[str, set[str]]:
    if path is None:
        return {}
    candidates: dict[str, set[str]] = defaultdict(set)
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"event_id", "candidate_te"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path}: invalid TE_CALL_CANDIDATES header")
        for row in reader:
            candidates[row["event_id"]].add(row["candidate_te"])
    return candidates


def numeric_evalue(value: object) -> float:
    result = float(value)
    return result if math.isfinite(result) else float("inf")


def build_matches(
    hsps: list[dict],
    final_calls: dict[str, str],
    candidate_calls: dict[str, set[str]],
    min_pident: float,
    min_aligned_bp: int,
    min_consensus_coverage: float,
) -> list[dict]:
    grouped: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for hsp in hsps:
        grouped[(hsp["query_id"], hsp["subject_te"], hsp["subject_strand"])].append(hsp)
    matches = []
    for (query_id, subject_te, strand), values in sorted(grouped.items()):
        first = values[0]
        usable = [value for value in values if value["retained"] == "yes"]
        aggregate = usable or values
        query_intervals = merge_intervals(
            (int(value["query_start"]), int(value["query_end"])) for value in aggregate
        )
        subject_intervals = merge_intervals(
            (int(value["subject_start"]), int(value["subject_end"])) for value in aggregate
        )
        aligned_query_bp = interval_length(query_intervals)
        aligned_subject_bp = interval_length(subject_intervals)
        weight = sum(int(value["alignment_length"]) for value in aggregate)
        identity = sum(
            float(value["pident"]) * int(value["alignment_length"])
            for value in aggregate
        ) / weight
        query_coverage = 100 * aligned_query_bp / int(first["query_length"])
        consensus_coverage = 100 * aligned_subject_bp / int(first["subject_length"])
        reasons = []
        if not usable:
            reasons.append("NO_USABLE_HSP")
        if identity < min_pident:
            reasons.append("LOW_IDENTITY")
        if aligned_subject_bp < min_aligned_bp:
            reasons.append("SHORT_MATCH")
        if consensus_coverage < min_consensus_coverage:
            reasons.append("LOW_CONSENSUS_COVERAGE")
        event_id = str(first["event_id"])
        reported_te = final_calls.get(event_id, ".")
        if subject_te == reported_te:
            assignment = "reported_primary"
        elif subject_te in candidate_calls.get(event_id, set()):
            assignment = "known_alternative"
        else:
            assignment = "unreported_match"
        matches.append(
            {
                "match_id": stable_id("M", query_id, subject_te, strand),
                "query_id": query_id,
                "event_id": event_id,
                "chrom": first["chrom"],
                "locus_start": first["locus_start"],
                "locus_end": first["locus_end"],
                "subject_te": subject_te,
                "subject_strand": strand,
                "query_start": min(start for start, _ in query_intervals),
                "query_end": max(end for _, end in query_intervals),
                "query_segments": serialize_intervals(query_intervals),
                "aligned_query_bp": aligned_query_bp,
                "aligned_subject_bp": aligned_subject_bp,
                "query_coverage": round(query_coverage, 6),
                "consensus_coverage": round(consensus_coverage, 6),
                "weighted_identity": round(identity, 6),
                "bitscore_sum": round(sum(float(value["bitscore"]) for value in aggregate), 6),
                "best_evalue": min(numeric_evalue(value["evalue"]) for value in aggregate),
                "hsp_count": len(aggregate),
                "query_length": first["query_length"],
                "subject_length": first["subject_length"],
                "status": "retained" if not reasons else "rejected",
                "filter_reasons": ";".join(reasons) if reasons else ".",
                "final_call": "yes" if event_id in final_calls else "no",
                "reported_te": reported_te,
                "assignment": assignment,
                "_query_intervals": query_intervals,
            }
        )
    return matches


def interval_overlap(
    first: Iterable[tuple[int, int]], second: Iterable[tuple[int, int]]
) -> int:
    left = merge_intervals(first)
    right = merge_intervals(second)
    i = j = total = 0
    while i < len(left) and j < len(right):
        total += max(0, min(left[i][1], right[j][1]) - max(left[i][0], right[j][0]))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return total


def compatible_components(first: dict, second: dict, overlap_fraction: float) -> bool:
    overlap = interval_overlap(first["_query_intervals"], second["_query_intervals"])
    denominator = min(int(first["aligned_query_bp"]), int(second["aligned_query_bp"]))
    return overlap / denominator <= overlap_fraction


def choose_components(matches: list[dict], overlap_fraction: float) -> list[dict]:
    ordered = sorted(
        (match for match in matches if match["status"] == "retained"),
        key=lambda match: (
            -float(match["bitscore_sum"]),
            -float(match["consensus_coverage"]),
            match["subject_te"],
        ),
    )
    selected = []
    for match in ordered:
        if all(compatible_components(match, other, overlap_fraction) for other in selected):
            selected.append(match)
    return sorted(selected, key=lambda match: (int(match["query_start"]), int(match["query_end"]), match["subject_te"]))


def classify_structure(selected: list[dict]) -> str:
    if not selected:
        return "unresolved"
    if len(selected) == 1:
        return "single_candidate"
    for index, first in enumerate(selected):
        for second in selected[index + 1 :]:
            if (
                int(first["query_start"]) <= int(second["query_start"])
                and int(first["query_end"]) >= int(second["query_end"])
            ) or (
                int(second["query_start"]) <= int(first["query_start"])
                and int(second["query_end"]) >= int(first["query_end"])
            ):
                return "nested_candidate"
    return "composite_candidate"


def build_structures(
    matches: list[dict],
    query_sizes: dict[str, int],
    final_calls: dict[str, str],
    overlap_fraction: float,
) -> tuple[list[dict], list[dict]]:
    by_query: dict[str, list[dict]] = defaultdict(list)
    for match in matches:
        by_query[match["query_id"]].append(match)
    structures = []
    components = []
    for query_id in sorted(query_sizes):
        query_matches = by_query.get(query_id, [])
        chrom, locus_start, locus_end, event_id = parse_query_id(query_id)
        selected = choose_components(query_matches, overlap_fraction)
        structure_id = stable_id("S", query_id)
        retained_count = sum(match["status"] == "retained" for match in query_matches)
        classification = classify_structure(selected) if query_matches else "no_match"
        if not query_matches:
            interpretation = "no_blast_alignment"
        elif len(selected) > 1:
            interpretation = "provisional_coordinate_supported"
        else:
            interpretation = "candidate_only"
        structures.append(
            {
                "structure_id": structure_id,
                "query_id": query_id,
                "event_id": event_id,
                "chrom": chrom,
                "locus_start": locus_start,
                "locus_end": locus_end,
                "query_length": query_sizes[query_id],
                "final_call": "yes" if event_id in final_calls else "no",
                "reported_te": final_calls.get(event_id, "."),
                "retained_matches": retained_count,
                "component_count": len(selected),
                "alternative_matches": max(0, retained_count - len(selected)),
                "classification": classification,
                "interpretation": interpretation,
            }
        )
        for rank, match in enumerate(selected, 1):
            components.append(
                {
                    "component_id": stable_id("C", structure_id, rank, match["match_id"]),
                    "structure_id": structure_id,
                    "query_id": query_id,
                    "event_id": event_id,
                    "rank": rank,
                    "match_id": match["match_id"],
                    "te_name": match["subject_te"],
                    "query_start": match["query_start"],
                    "query_end": match["query_end"],
                    "query_segments": match["query_segments"],
                    "strand": match["subject_strand"],
                    "consensus_coverage": match["consensus_coverage"],
                    "weighted_identity": match["weighted_identity"],
                    "assignment": match["assignment"],
                    "interpretation": "provisional_component",
                }
            )
    return structures, components


def build_events(structures: list[dict], components: list[dict]) -> list[dict]:
    structures_by_event: dict[str, list[dict]] = defaultdict(list)
    for structure in structures:
        structures_by_event[structure["event_id"]].append(structure)
    component_families: dict[str, set[str]] = defaultdict(set)
    for component in components:
        component_families[component["structure_id"]].add(component["te_name"])
    events = []
    for event_id, rows in sorted(structures_by_event.items()):
        first = rows[0]
        invariants = ("chrom", "locus_start", "locus_end", "final_call", "reported_te")
        if any(
            row[field] != first[field]
            for row in rows[1:]
            for field in invariants
        ):
            raise ValueError(f"conflicting query metadata for event {event_id}")
        multi_component = sum(int(row["component_count"]) > 1 for row in rows)
        multi_family = sum(
            len(component_families.get(row["structure_id"], set())) > 1
            for row in rows
        )
        maximum = max(int(row["component_count"]) for row in rows)
        if multi_family >= 2:
            classification = "multi_te_supported"
        elif multi_family == 1:
            classification = "multi_te_candidate"
        elif multi_component >= 2:
            classification = "repeated_or_rearranged_te_supported"
        elif multi_component == 1:
            classification = "repeated_or_rearranged_te_candidate"
        elif any(int(row["component_count"]) == 1 for row in rows):
            classification = "single_candidate"
        elif any(row["classification"] == "unresolved" for row in rows):
            classification = "unresolved"
        else:
            classification = "no_match"
        events.append(
            {
                "event_structure_id": stable_id("E", event_id),
                "event_id": event_id,
                "chrom": first["chrom"],
                "locus_start": min(int(row["locus_start"]) for row in rows),
                "locus_end": max(int(row["locus_end"]) for row in rows),
                "final_call": first["final_call"],
                "reported_te": first["reported_te"],
                "query_count": len(rows),
                "queries_with_matches": sum(
                    int(row["retained_matches"]) > 0 for row in rows
                ),
                "queries_with_components": sum(
                    int(row["component_count"]) > 0 for row in rows
                ),
                "multi_component_queries": multi_component,
                "multi_family_queries": multi_family,
                "maximum_component_count": maximum,
                "classification": classification,
                "interpretation": (
                    "coordinate_supported_across_queries"
                    if classification.endswith("_supported")
                    else "candidate_only"
                ),
            }
        )
    return events


def public_row(row: dict, columns: tuple[str, ...]) -> dict:
    return {column: row[column] for column in columns}


def render_qmd(template: str, style: str, script: str, data: dict) -> str:
    replacements = {
        "<!-- STRUCTURE_REPORT_STYLE -->": f"<style>\n{style}\n</style>",
        "<!-- STRUCTURE_REPORT_DATA -->": (
            '<script id="structure-report-data" type="application/json">'
            + json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace("<", "\\u003c")
            + "</script>"
        ),
        "<!-- STRUCTURE_REPORT_SCRIPT -->": f"<script>\n{script}\n</script>",
        '"STRUCTURE_REPORT_TITLE"': json.dumps(data["report"]["title"]),
    }
    for placeholder in replacements:
        if template.count(placeholder) != 1:
            raise ValueError(f"template requires exactly one {placeholder}")
    for placeholder, replacement in replacements.items():
        template = template.replace(placeholder, replacement)
    return template


def build(
    blast: Path,
    te_sizes: Path,
    query_sizes: Path,
    output: Path,
    te_infos: Optional[Path],
    te_call_candidates: Optional[Path],
    min_pident: float,
    min_aligned_bp: int,
    min_consensus_coverage: float,
    max_component_overlap: float,
    template: Path,
    style: Path,
    script: Path,
    title: str,
) -> dict:
    if not 0 <= min_pident <= 100:
        raise ValueError("minimum identity must be in [0,100]")
    if min_aligned_bp < 1:
        raise ValueError("minimum aligned bp must be positive")
    if not 0 <= min_consensus_coverage <= 100:
        raise ValueError("minimum consensus coverage must be in [0,100]")
    if not 0 <= max_component_overlap < 1:
        raise ValueError("maximum component overlap must be in [0,1)")
    output.mkdir(parents=True, exist_ok=True)
    subject_sizes = read_sizes(te_sizes, "TE")
    insertion_sizes = read_sizes(query_sizes, "query")
    final_calls = read_final_calls(te_infos)
    candidate_calls = read_candidate_calls(te_call_candidates)
    hsps = read_hsps(blast, insertion_sizes, subject_sizes, min_pident, min_aligned_bp)
    matches = build_matches(
        hsps,
        final_calls,
        candidate_calls,
        min_pident,
        min_aligned_bp,
        min_consensus_coverage,
    )
    structures, components = build_structures(
        matches, insertion_sizes, final_calls, max_component_overlap
    )
    events = build_events(structures, components)
    public_hsps = [public_row(row, HSP_COLUMNS) for row in hsps]
    public_matches = [public_row(row, MATCH_COLUMNS) for row in matches]
    write_tsv(output / "insertion-hsps.tsv", HSP_COLUMNS, public_hsps)
    write_tsv(output / "insertion-matches.tsv", MATCH_COLUMNS, public_matches)
    write_tsv(output / "insertion-components.tsv", COMPONENT_COLUMNS, components)
    write_tsv(output / "insertion-structures.tsv", STRUCTURE_COLUMNS, structures)
    write_tsv(output / "insertion-events.tsv", EVENT_COLUMNS, events)
    classification_counts: dict[str, int] = defaultdict(int)
    for row in structures:
        classification_counts[row["classification"]] += 1
    data = {
        "schema_version": SCHEMA_VERSION,
        "report": {
            "title": title,
            "thresholds": {
                "min_pident": min_pident,
                "min_aligned_bp": min_aligned_bp,
                "min_consensus_coverage": min_consensus_coverage,
                "max_component_overlap": max_component_overlap,
            },
            "inputs": {
                "blast": str(blast.resolve()),
                "blast_sha256": file_sha256(blast),
                "te_sizes": str(te_sizes.resolve()),
                "query_sizes": str(query_sizes.resolve()),
                "te_infos": str(te_infos.resolve()) if te_infos else None,
                "te_call_candidates": str(te_call_candidates.resolve()) if te_call_candidates else None,
            },
        },
        "summary": {
            "queries": len(structures),
            "events": len(events),
            "hsps": len(hsps),
            "retained_hsps": sum(row["retained"] == "yes" for row in hsps),
            "matches": len(matches),
            "retained_matches": sum(row["status"] == "retained" for row in matches),
            "components": len(components),
            "final_queries": sum(row["final_call"] == "yes" for row in structures),
            "final_events": sum(row["final_call"] == "yes" for row in events),
            "classifications": dict(sorted(classification_counts.items())),
        },
        "structures": structures,
        "events": events,
        "components": components,
        "matches": public_matches,
        "hsps": public_hsps,
    }
    manifest_data = dict(data["report"])
    manifest_data["schema_version"] = SCHEMA_VERSION
    manifest_data["summary"] = data["summary"]
    atomic_write_text(
        output / "manifest.json",
        json.dumps(manifest_data, indent=2, sort_keys=True) + "\n",
    )
    atomic_write_text(
        output / "report-data.json",
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    atomic_write_text(
        output / "report.qmd",
        render_qmd(template.read_text(), style.read_text(), script.read_text(), data),
    )
    return data


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", required=True, type=Path)
    parser.add_argument("--te-sizes", required=True, type=Path)
    parser.add_argument("--query-sizes", required=True, type=Path)
    parser.add_argument("--te-infos", type=Path)
    parser.add_argument("--te-call-candidates", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--min-pident", type=float, default=90.0)
    parser.add_argument("--min-aligned-bp", type=int, default=80)
    parser.add_argument("--min-consensus-coverage", type=float, default=20.0)
    parser.add_argument("--max-component-overlap", type=float, default=0.2)
    parser.add_argument(
        "--template",
        type=Path,
        default=PIPELINE_ROOT / "modules/2-MODULE_TE_BLAST/report.qmd",
    )
    parser.add_argument(
        "--style",
        type=Path,
        default=PIPELINE_ROOT / "modules/2-MODULE_TE_BLAST/report.css",
    )
    parser.add_argument(
        "--script",
        type=Path,
        default=PIPELINE_ROOT / "modules/2-MODULE_TE_BLAST/report.js",
    )
    parser.add_argument("--title", default="TrEMOLO insertion structure explorer")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    data = build(
        blast=args.blast,
        te_sizes=args.te_sizes,
        query_sizes=args.query_sizes,
        output=args.output,
        te_infos=args.te_infos,
        te_call_candidates=args.te_call_candidates,
        min_pident=args.min_pident,
        min_aligned_bp=args.min_aligned_bp,
        min_consensus_coverage=args.min_consensus_coverage,
        max_component_overlap=args.max_component_overlap,
        template=args.template,
        style=args.style,
        script=args.script,
        title=args.title,
    )
    print(
        "Insertion structure report prepared: "
        f"{data['summary']['queries']} queries, "
        f"{data['summary']['retained_matches']} retained matches, "
        f"{data['summary']['components']} provisional components."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
