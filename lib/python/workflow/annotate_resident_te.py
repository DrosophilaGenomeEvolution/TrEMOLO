#!/usr/bin/env python3
"""Build a lossless resident-TE annotation from one TE-versus-genome BLAST."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import re
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import quote


BLAST_COLUMNS = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
)

FRAGMENT_COLUMNS = (
    "fragment_id",
    "target",
    "chrom",
    "start",
    "end",
    "strand",
    "te_name",
    "consensus_start",
    "consensus_end",
    "consensus_length",
    "aligned_bp",
    "identity",
    "mismatches",
    "gap_opens",
    "evalue",
    "bitscore",
    "filter_status",
    "filter_reason",
    "redundant_with",
    "match_id",
)

MATCH_COLUMNS = (
    "match_id",
    "copy_id",
    "target",
    "chrom",
    "start",
    "end",
    "strand",
    "te_name",
    "consensus_length",
    "consensus_covered_bp",
    "consensus_coverage",
    "identity",
    "bitscore",
    "best_evalue",
    "fragment_count",
    "tier",
    "assignment",
)

COPY_COLUMNS = (
    "copy_id",
    "target",
    "chrom",
    "start",
    "end",
    "strand",
    "primary_te",
    "te_candidates",
    "consensus_coverage",
    "identity",
    "aligned_bp",
    "fragment_count",
    "candidate_count",
    "status",
    "structure",
)

RELATION_COLUMNS = (
    "relation_id",
    "target",
    "source_copy_id",
    "target_copy_id",
    "relation",
    "overlap_bp",
    "confidence",
)


@dataclass
class Fragment:
    query_id: str
    te_name: str
    chrom: str
    pident: float
    aligned_bp: int
    mismatches: int
    gap_opens: int
    consensus_start: int
    consensus_end: int
    start: int
    end: int
    strand: str
    evalue: float
    bitscore: float
    consensus_length: int
    target_length: int
    filter_status: str
    filter_reason: str
    fragment_id: str = "."
    redundant_with: str = "."
    match_id: str = "."

    @property
    def oriented_consensus_start(self) -> int:
        if self.strand == "+":
            return self.consensus_start
        return self.consensus_length - self.consensus_end

    @property
    def oriented_consensus_end(self) -> int:
        if self.strand == "+":
            return self.consensus_end
        return self.consensus_length - self.consensus_start


@dataclass
class Match:
    match_id: str
    target: str
    chrom: str
    start: int
    end: int
    strand: str
    te_name: str
    consensus_length: int
    consensus_covered_bp: int
    consensus_coverage: float
    identity: float
    bitscore: float
    best_evalue: float
    fragments: List[Fragment]
    tier: str
    copy_id: str = "."
    assignment: str = "."


@dataclass
class Copy:
    copy_id: str
    target: str
    chrom: str
    start: int
    end: int
    strand: str
    primary_te: str
    te_candidates: List[str]
    consensus_coverage: float
    identity: float
    aligned_bp: int
    fragment_count: int
    candidate_count: int
    status: str
    structure: str
    matches: List[Match] = field(default_factory=list)


def format_number(value: float, digits: int = 4) -> str:
    if math.isinf(value):
        return "inf"
    return f"{value:.{digits}f}".rstrip("0").rstrip(".")


def format_significant(value: float, digits: int = 6) -> str:
    if math.isinf(value):
        return "inf"
    return f"{value:.{digits}g}"


def stable_id(prefix: str, values: Iterable[object]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(str(value).encode("utf-8"))
        digest.update(b"\x1f")
    return prefix + digest.hexdigest()[:16]


def read_te_names(path: Path) -> Dict[str, str]:
    result: Dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["original", "prepared"]:
            raise ValueError(f"{path}: unexpected TE-header schema: {reader.fieldnames!r}")
        for line_number, row in enumerate(reader, 2):
            original = row["original"]
            prepared = row["prepared"]
            if not original or not prepared:
                raise ValueError(f"{path}:{line_number}: empty TE identifier")
            aliases = {prepared, prepared.split()[0]}
            for alias in aliases:
                previous = result.get(alias)
                if previous is not None and previous != original:
                    raise ValueError(
                        f"{path}:{line_number}: ambiguous prepared TE identifier {alias!r}"
                    )
                result[alias] = original
    return result


def parse_integer(value: str, label: str, path: Path, line_number: int) -> int:
    try:
        return int(value)
    except ValueError as error:
        raise ValueError(f"{path}:{line_number}: invalid {label}: {value!r}") from error


def parse_float(value: str, label: str, path: Path, line_number: int) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ValueError(f"{path}:{line_number}: invalid {label}: {value!r}") from error
    if math.isnan(result):
        raise ValueError(f"{path}:{line_number}: invalid {label}: {value!r}")
    return result


def read_fragments(
    path: Path,
    te_names: Dict[str, str],
    target: str,
    chrom_pattern: re.Pattern,
    min_pident: float,
    min_aligned_bp: int,
    max_evalue: float,
) -> List[Fragment]:
    fragments = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, 1):
            if not fields or not any(fields):
                continue
            if len(fields) != len(BLAST_COLUMNS):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(BLAST_COLUMNS)} BLAST "
                    f"columns, got {len(fields)}"
                )
            values = dict(zip(BLAST_COLUMNS, fields))
            query_id = values["qseqid"]
            if query_id not in te_names:
                raise ValueError(
                    f"{path}:{line_number}: BLAST query {query_id!r} is absent "
                    "from the prepared TE-header table"
                )
            pident = parse_float(values["pident"], "identity", path, line_number)
            aligned_bp = parse_integer(values["length"], "alignment length", path, line_number)
            mismatches = parse_integer(values["mismatch"], "mismatch count", path, line_number)
            gap_opens = parse_integer(values["gapopen"], "gap-open count", path, line_number)
            qstart = parse_integer(values["qstart"], "query start", path, line_number)
            qend = parse_integer(values["qend"], "query end", path, line_number)
            sstart = parse_integer(values["sstart"], "subject start", path, line_number)
            send = parse_integer(values["send"], "subject end", path, line_number)
            evalue = parse_float(values["evalue"], "e-value", path, line_number)
            bitscore = parse_float(values["bitscore"], "bit score", path, line_number)
            qlen = parse_integer(values["qlen"], "query length", path, line_number)
            slen = parse_integer(values["slen"], "subject length", path, line_number)
            if (
                aligned_bp <= 0
                or qlen <= 0
                or slen <= 0
                or min(qstart, qend, sstart, send) < 1
                or max(qstart, qend) > qlen
                or max(sstart, send) > slen
            ):
                raise ValueError(f"{path}:{line_number}: BLAST coordinates are out of bounds")

            reasons = []
            if chrom_pattern.search(values["sseqid"]) is None:
                reasons.append("chrom_filtered")
            if pident < min_pident:
                reasons.append("identity")
            if aligned_bp < min_aligned_bp:
                reasons.append("aligned_bp")
            if evalue > max_evalue:
                reasons.append("evalue")

            fragments.append(
                Fragment(
                    query_id=query_id,
                    te_name=te_names[query_id],
                    chrom=values["sseqid"],
                    pident=pident,
                    aligned_bp=aligned_bp,
                    mismatches=mismatches,
                    gap_opens=gap_opens,
                    consensus_start=min(qstart, qend) - 1,
                    consensus_end=max(qstart, qend),
                    start=min(sstart, send) - 1,
                    end=max(sstart, send),
                    strand="+" if sstart <= send else "-",
                    evalue=evalue,
                    bitscore=bitscore,
                    consensus_length=qlen,
                    target_length=slen,
                    filter_status="accepted" if not reasons else "rejected",
                    filter_reason=";".join(reasons) if reasons else ".",
                )
            )

    fragments.sort(
        key=lambda row: (
            row.chrom,
            row.start,
            row.end,
            row.te_name,
            row.strand,
            row.consensus_start,
            row.consensus_end,
            -row.bitscore,
            row.evalue,
            row.mismatches,
            row.gap_opens,
        )
    )
    duplicate_ids = defaultdict(int)
    for fragment in fragments:
        identity = (
            target,
            fragment.query_id,
            fragment.chrom,
            fragment.start,
            fragment.end,
            fragment.strand,
            fragment.consensus_start,
            fragment.consensus_end,
            format_number(fragment.pident),
            format_number(fragment.bitscore),
            format_significant(fragment.evalue),
            fragment.mismatches,
            fragment.gap_opens,
        )
        duplicate_ids[identity] += 1
        fragment.fragment_id = stable_id(
            "TF", identity + (duplicate_ids[identity],)
        )
    return fragments


def overlap_length(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    return max(0, min(end_a, end_b) - max(start_a, start_b))


def reciprocal_overlap(
    start_a: int, end_a: int, start_b: int, end_b: int
) -> Tuple[float, float]:
    overlap = overlap_length(start_a, end_a, start_b, end_b)
    return overlap / (end_a - start_a), overlap / (end_b - start_b)


def choose_nonredundant(fragments: Sequence[Fragment]) -> List[Fragment]:
    selected = []
    genome_bins = defaultdict(list)
    bin_size = 1000
    for fragment in sorted(
        fragments,
        key=lambda row: (-row.bitscore, -row.aligned_bp, row.start, row.consensus_start),
    ):
        duplicate = None
        candidate_indexes = set()
        for bin_index in range(fragment.start // bin_size, (fragment.end - 1) // bin_size + 1):
            candidate_indexes.update(genome_bins[bin_index])
        for selected_index in sorted(candidate_indexes):
            previous = selected[selected_index]
            genome_overlap = reciprocal_overlap(
                fragment.start, fragment.end, previous.start, previous.end
            )
            consensus_overlap = reciprocal_overlap(
                fragment.consensus_start,
                fragment.consensus_end,
                previous.consensus_start,
                previous.consensus_end,
            )
            if min(genome_overlap) >= 0.8 and min(consensus_overlap) >= 0.8:
                duplicate = previous
                break
        if duplicate is None:
            selected_index = len(selected)
            selected.append(fragment)
            for bin_index in range(
                fragment.start // bin_size, (fragment.end - 1) // bin_size + 1
            ):
                genome_bins[bin_index].append(selected_index)
        else:
            fragment.redundant_with = duplicate.fragment_id
    return selected


def compatible(
    previous: Fragment,
    current: Fragment,
    min_chain_gap: int,
    max_gap_fraction: float,
    max_overlap: int,
) -> Optional[Tuple[int, int]]:
    genome_gap = current.start - previous.end
    consensus_gap = (
        current.oriented_consensus_start - previous.oriented_consensus_end
    )
    gap_limit = max(min_chain_gap, int(round(current.consensus_length * max_gap_fraction)))
    if genome_gap < -max_overlap or consensus_gap < -max_overlap:
        return None
    if genome_gap > gap_limit or consensus_gap > gap_limit:
        return None
    return genome_gap, consensus_gap


def interval_union_length(intervals: Iterable[Tuple[int, int]]) -> int:
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    start, end = ordered[0]
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def classify_match(
    coverage: float,
    identity: float,
    full_length_coverage: float,
    high_confidence_pident: float,
    partial_coverage: float,
) -> str:
    if coverage >= full_length_coverage and identity >= high_confidence_pident:
        return "full_length"
    if coverage >= partial_coverage:
        return "partial"
    return "degraded_relic"


def build_matches(
    fragments: Sequence[Fragment],
    target: str,
    min_chain_gap: int,
    max_gap_fraction: float,
    max_overlap: int,
    full_length_coverage: float,
    high_confidence_pident: float,
    partial_coverage: float,
) -> List[Match]:
    groups = defaultdict(list)
    for fragment in fragments:
        if fragment.filter_status == "accepted":
            groups[
                (
                    fragment.chrom,
                    fragment.te_name,
                    fragment.strand,
                    fragment.consensus_length,
                )
            ].append(fragment)

    chains: List[List[Fragment]] = []
    for group_key in sorted(groups):
        nonredundant = choose_nonredundant(groups[group_key])
        family_chains: List[List[Fragment]] = []
        for fragment in sorted(
            nonredundant,
            key=lambda row: (
                row.start,
                row.oriented_consensus_start,
                row.end,
                row.oriented_consensus_end,
            ),
        ):
            candidates = []
            for chain_index, chain in enumerate(family_chains):
                gaps = compatible(
                    chain[-1],
                    fragment,
                    min_chain_gap,
                    max_gap_fraction,
                    max_overlap,
                )
                if gaps is not None:
                    candidates.append(
                        (
                            abs(max(0, gaps[0]) - max(0, gaps[1])),
                            max(gaps),
                            chain_index,
                        )
                    )
            if candidates:
                family_chains[min(candidates)[2]].append(fragment)
            else:
                family_chains.append([fragment])
        chains.extend(family_chains)

    chains.sort(
        key=lambda chain: (
            chain[0].chrom,
            min(row.start for row in chain),
            max(row.end for row in chain),
            chain[0].te_name,
            chain[0].strand,
        )
    )
    matches = []
    for chain in chains:
        covered = interval_union_length(
            (row.consensus_start, row.consensus_end) for row in chain
        )
        consensus_length = chain[0].consensus_length
        coverage = covered * 100.0 / consensus_length
        aligned_total = sum(row.aligned_bp for row in chain)
        identity = sum(row.pident * row.aligned_bp for row in chain) / aligned_total
        match_id = stable_id(
            "TM",
            (
                target,
                chain[0].chrom,
                min(row.start for row in chain),
                max(row.end for row in chain),
                chain[0].strand,
                chain[0].te_name,
                *(row.fragment_id for row in chain),
            ),
        )
        match = Match(
            match_id=match_id,
            target=target,
            chrom=chain[0].chrom,
            start=min(row.start for row in chain),
            end=max(row.end for row in chain),
            strand=chain[0].strand,
            te_name=chain[0].te_name,
            consensus_length=consensus_length,
            consensus_covered_bp=covered,
            consensus_coverage=coverage,
            identity=identity,
            bitscore=sum(row.bitscore for row in chain),
            best_evalue=min(row.evalue for row in chain),
            fragments=list(chain),
            tier=classify_match(
                coverage,
                identity,
                full_length_coverage,
                high_confidence_pident,
                partial_coverage,
            ),
        )
        for fragment in chain:
            fragment.match_id = match.match_id
        matches.append(match)

    by_id = {row.fragment_id: row for row in fragments}
    for fragment in fragments:
        if fragment.redundant_with != ".":
            fragment.match_id = by_id[fragment.redundant_with].match_id
    return matches


def primary_key(match: Match) -> Tuple[int, int, float, float, float, str]:
    tier_rank = {"degraded_relic": 0, "partial": 1, "full_length": 2}
    return (
        tier_rank[match.tier],
        match.consensus_covered_bp,
        match.consensus_coverage,
        match.identity,
        match.bitscore,
        match.te_name,
    )


def build_copies(
    matches: Sequence[Match], target: str, ambiguity_overlap: float
) -> List[Copy]:
    # Complete-linkage grouping prevents a chain of partially overlapping
    # matches from collapsing several neighbouring copies into one component.
    ordered_matches = sorted(
        matches, key=lambda row: (row.chrom, row.start, row.end, row.te_name)
    )
    ordered_components: List[List[Match]] = []
    components_by_chrom = defaultdict(list)
    for match in ordered_matches:
        compatible_groups = []
        for group_index in components_by_chrom[match.chrom]:
            group = ordered_components[group_index]
            if max(row.end for row in group) <= match.start:
                continue
            if all(
                min(reciprocal_overlap(match.start, match.end, row.start, row.end))
                >= ambiguity_overlap
                for row in group
            ):
                compatible_groups.append(group_index)
        if compatible_groups:
            ordered_components[compatible_groups[0]].append(match)
        else:
            ordered_components.append([match])
            components_by_chrom[match.chrom].append(len(ordered_components) - 1)

    copies = []
    for group in ordered_components:
        copy_id = stable_id(
            "TC",
            (
                target,
                group[0].chrom,
                min(row.start for row in group),
                max(row.end for row in group),
                *(row.match_id for row in sorted(group, key=lambda row: row.match_id)),
            ),
        )
        primary = max(group, key=primary_key)
        candidates = sorted({row.te_name for row in group})
        for match in group:
            match.copy_id = copy_id
            match.assignment = "primary" if match is primary else "alternative"
        if len(group) > 1:
            structure = "ambiguous"
        elif len(primary.fragments) > 1:
            structure = "fragmented"
        else:
            structure = "single"
        copies.append(
            Copy(
                copy_id=copy_id,
                target=target,
                chrom=primary.chrom,
                start=min(row.start for row in group),
                end=max(row.end for row in group),
                strand=primary.strand,
                primary_te=primary.te_name,
                te_candidates=candidates,
                consensus_coverage=primary.consensus_coverage,
                identity=primary.identity,
                aligned_bp=primary.consensus_covered_bp,
                fragment_count=len(primary.fragments),
                candidate_count=len(group),
                status=primary.tier,
                structure=structure,
                matches=sorted(group, key=lambda row: (row.assignment != "primary", row.te_name)),
            )
        )
    return copies


def build_relations(copies: Sequence[Copy], target: str) -> List[Tuple[str, ...]]:
    relations = []
    by_chrom = defaultdict(list)
    for copy in copies:
        by_chrom[copy.chrom].append(copy)
    for chrom in sorted(by_chrom):
        ordered = sorted(by_chrom[chrom], key=lambda row: (row.start, row.end))
        for position, left in enumerate(ordered):
            for right in ordered[position + 1 :]:
                if right.start >= left.end:
                    break
                overlap = overlap_length(left.start, left.end, right.start, right.end)
                if left.primary_te == right.primary_te:
                    source, destination, relation = left, right, "same_family_overlap"
                elif left.start <= right.start and right.end <= left.end:
                    source, destination, relation = right, left, "nested_candidate"
                elif right.start <= left.start and left.end <= right.end:
                    source, destination, relation = left, right, "nested_candidate"
                else:
                    source, destination, relation = left, right, "overlap_candidate"
                relation_id = stable_id(
                    "TR",
                    (target, source.copy_id, destination.copy_id, relation),
                )
                relations.append(
                    (
                        relation_id,
                        target,
                        source.copy_id,
                        destination.copy_id,
                        relation,
                        str(overlap),
                        "provisional",
                    )
                )
    relations.sort(key=lambda row: (row[2], row[3], row[4]))
    return relations


def fragment_rows(fragments: Sequence[Fragment], target: str) -> Iterable[Sequence[str]]:
    for row in fragments:
        yield (
            row.fragment_id,
            target,
            row.chrom,
            str(row.start),
            str(row.end),
            row.strand,
            row.te_name,
            str(row.consensus_start),
            str(row.consensus_end),
            str(row.consensus_length),
            str(row.aligned_bp),
            format_number(row.pident),
            str(row.mismatches),
            str(row.gap_opens),
            format_significant(row.evalue),
            format_number(row.bitscore),
            row.filter_status,
            row.filter_reason,
            row.redundant_with,
            row.match_id,
        )


def match_rows(matches: Sequence[Match]) -> Iterable[Sequence[str]]:
    for row in sorted(matches, key=lambda item: item.match_id):
        yield (
            row.match_id,
            row.copy_id,
            row.target,
            row.chrom,
            str(row.start),
            str(row.end),
            row.strand,
            row.te_name,
            str(row.consensus_length),
            str(row.consensus_covered_bp),
            format_number(row.consensus_coverage),
            format_number(row.identity),
            format_number(row.bitscore),
            format_significant(row.best_evalue),
            str(len(row.fragments)),
            row.tier,
            row.assignment,
        )


def copy_rows(copies: Sequence[Copy]) -> Iterable[Sequence[str]]:
    for row in copies:
        yield (
            row.copy_id,
            row.target,
            row.chrom,
            str(row.start),
            str(row.end),
            row.strand,
            row.primary_te,
            ";".join(row.te_candidates),
            format_number(row.consensus_coverage),
            format_number(row.identity),
            str(row.aligned_bp),
            str(row.fragment_count),
            str(row.candidate_count),
            row.status,
            row.structure,
        )


def bed_rows(copies: Sequence[Copy]) -> Iterable[Sequence[str]]:
    for row in copies:
        score = min(1000, max(0, int(round(row.identity * row.consensus_coverage / 10))))
        yield (
            row.chrom,
            str(row.start),
            str(row.end),
            f"{quote(row.primary_te, safe='._:-')}|{row.copy_id}",
            str(score),
            row.strand,
        )


def gff_attributes(values: Sequence[Tuple[str, object]]) -> str:
    return ";".join(f"{key}={quote(str(value), safe='._:-')}" for key, value in values)


def gff_rows(copies: Sequence[Copy]) -> Iterable[Sequence[str]]:
    for copy in copies:
        yield (
            copy.chrom,
            "TrEMOLO",
            "transposable_element",
            str(copy.start + 1),
            str(copy.end),
            ".",
            copy.strand,
            ".",
            gff_attributes(
                (
                    ("ID", copy.copy_id),
                    ("Name", copy.primary_te),
                    ("status", copy.status),
                    ("structure", copy.structure),
                    ("candidate_te", ",".join(copy.te_candidates)),
                )
            ),
        )
        for match in copy.matches:
            yield (
                match.chrom,
                "TrEMOLO",
                "nucleotide_match",
                str(match.start + 1),
                str(match.end),
                format_number(match.bitscore),
                match.strand,
                ".",
                gff_attributes(
                    (
                        ("ID", match.match_id),
                        ("Parent", copy.copy_id),
                        ("Name", match.te_name),
                        ("assignment", match.assignment),
                        ("coverage", format_number(match.consensus_coverage)),
                        ("identity", format_number(match.identity)),
                    )
                ),
            )


def write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def write_outputs(
    destinations: Dict[str, Path],
    fragments: Sequence[Fragment],
    matches: Sequence[Match],
    copies: Sequence[Copy],
    relations: Sequence[Sequence[str]],
    target: str,
) -> None:
    temporary = {}
    try:
        for label, destination in destinations.items():
            destination.parent.mkdir(parents=True, exist_ok=True)
            descriptor, name = tempfile.mkstemp(
                prefix=f".{destination.name}.", dir=str(destination.parent), text=True
            )
            os.close(descriptor)
            temporary[label] = Path(name)
        write_tsv(temporary["fragments"], FRAGMENT_COLUMNS, fragment_rows(fragments, target))
        write_tsv(temporary["matches"], MATCH_COLUMNS, match_rows(matches))
        write_tsv(temporary["copies"], COPY_COLUMNS, copy_rows(copies))
        write_tsv(temporary["relations"], RELATION_COLUMNS, relations)
        with temporary["bed"].open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(bed_rows(copies))
        with temporary["gff"].open("w", newline="") as handle:
            handle.write("##gff-version 3\n")
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerows(gff_rows(copies))
        for label, destination in destinations.items():
            os.replace(str(temporary[label]), str(destination))
    finally:
        for path in temporary.values():
            if path.exists():
                path.unlink()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("blast", type=Path)
    parser.add_argument("te_headers", type=Path)
    parser.add_argument("--target", required=True)
    parser.add_argument("--chrom-regex", default=".")
    parser.add_argument("--min-pident", type=float, default=65.0)
    parser.add_argument("--min-aligned-bp", type=int, default=80)
    parser.add_argument("--max-evalue", type=float, default=1e-10)
    parser.add_argument("--min-chain-gap", type=int, default=200)
    parser.add_argument("--max-gap-fraction", type=float, default=0.2)
    parser.add_argument("--max-overlap", type=int, default=100)
    parser.add_argument("--full-length-coverage", type=float, default=80.0)
    parser.add_argument("--high-confidence-pident", type=float, default=80.0)
    parser.add_argument("--partial-coverage", type=float, default=20.0)
    parser.add_argument("--ambiguity-overlap", type=float, default=0.8)
    parser.add_argument("--fragments", required=True, type=Path)
    parser.add_argument("--matches", required=True, type=Path)
    parser.add_argument("--copies", required=True, type=Path)
    parser.add_argument("--relations", required=True, type=Path)
    parser.add_argument("--bed", required=True, type=Path)
    parser.add_argument("--gff", required=True, type=Path)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    for label in ("min_pident", "full_length_coverage", "high_confidence_pident", "partial_coverage"):
        if not 0 <= getattr(args, label) <= 100:
            parser.error(f"{label.replace('_', '-')} must be between 0 and 100")
    if args.min_aligned_bp <= 0 or args.min_chain_gap < 0 or args.max_overlap < 0:
        parser.error("length and gap parameters must be non-negative, with min-aligned-bp > 0")
    if args.max_evalue < 0 or args.max_gap_fraction < 0:
        parser.error("max-evalue and max-gap-fraction must be non-negative")
    if not 0 < args.ambiguity_overlap <= 1:
        parser.error("ambiguity-overlap must be in (0, 1]")
    try:
        chrom_pattern = re.compile(args.chrom_regex)
    except re.error as error:
        parser.error(f"invalid chrom-regex: {error}")

    fragments = read_fragments(
        args.blast,
        read_te_names(args.te_headers),
        args.target,
        chrom_pattern,
        args.min_pident,
        args.min_aligned_bp,
        args.max_evalue,
    )
    matches = build_matches(
        fragments,
        args.target,
        args.min_chain_gap,
        args.max_gap_fraction,
        args.max_overlap,
        args.full_length_coverage,
        args.high_confidence_pident,
        args.partial_coverage,
    )
    copies = build_copies(matches, args.target, args.ambiguity_overlap)
    relations = build_relations(copies, args.target)
    write_outputs(
        {
            "fragments": args.fragments,
            "matches": args.matches,
            "copies": args.copies,
            "relations": args.relations,
            "bed": args.bed,
            "gff": args.gff,
        },
        fragments,
        matches,
        copies,
        relations,
        args.target,
    )
    accepted = sum(row.filter_status == "accepted" for row in fragments)
    ambiguous = sum(row.candidate_count > 1 for row in copies)
    print(
        f"resident TE annotation: {len(fragments)} HSPs, {accepted} accepted, "
        f"{len(matches)} matches, {len(copies)} copies, {ambiguous} ambiguous copies, "
        f"{len(relations)} provisional relations"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
