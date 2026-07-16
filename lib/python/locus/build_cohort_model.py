#!/usr/bin/env python3
"""Build a conservative multi-sample locus model from TE_INFOS.bed files."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from build_locus_model import (
    SCHEMA_VERSION,
    InputRun,
    file_sha256,
    parse_frequency,
    parse_number,
    parse_sample,
    read_calls,
    stable_id,
    write_tsv,
)


def anchor(call: dict[str, str]) -> int:
    value = call["new_position"]
    if value not in {"", ".", "NONE"}:
        try:
            return int(value)
        except ValueError:
            pass
    return int(call["start"])


def event_class(event_type: str) -> str:
    value = event_type.upper()
    if "DEL" in value or "CONTRACTION" in value:
        return "deletion"
    if "INS" in value or "EXPANSION" in value or value in {"HARD", "SOFT"}:
        return "insertion"
    return "complex"


def cluster_calls(calls: list[dict[str, str]], window: int) -> list[list[dict[str, str]]]:
    """Complete-span clustering; unlike single linkage, calls cannot chain indefinitely."""
    clusters = []
    by_chromosome: dict[str, list[dict[str, str]]] = defaultdict(list)
    for call in calls:
        by_chromosome[call["chrom"]].append(call)

    for chromosome in sorted(by_chromosome):
        ordered = sorted(
            by_chromosome[chromosome],
            key=lambda call: (anchor(call), int(call["start"]), call["sample_id"], call["tremolo_id"]),
        )
        current: list[dict[str, str]] = []
        minimum_anchor = 0
        for call in ordered:
            position = anchor(call)
            if not current or position - minimum_anchor <= window:
                if not current:
                    minimum_anchor = position
                current.append(call)
            else:
                clusters.append(current)
                current = [call]
                minimum_anchor = position
        if current:
            clusters.append(current)
    return clusters


def allele_signature(call: dict[str, str]) -> tuple[str, str, str]:
    te_name = call["te_call"].split("|", 1)[0]
    strand = call["strand"] if call["strand"] in {"+", "-"} else "."
    return event_class(call["event_type"]), te_name, strand


def build_cohort(
    runs: list[InputRun], output: Path, reference_id: str, locus_window: int
) -> dict[str, int]:
    if locus_window < 0:
        raise ValueError("locus window must be non-negative")
    calls = []
    seen_calls = set()
    for run in runs:
        if not run.te_infos.is_file():
            raise FileNotFoundError(run.te_infos)
        for call in read_calls(run):
            key = (run.sample_id, call["tremolo_id"])
            if key in seen_calls:
                raise ValueError(f"duplicate call id for sample {run.sample_id}: {call['tremolo_id']}")
            seen_calls.add(key)
            calls.append(call)

    loci_rows = []
    allele_rows = []
    component_rows = []
    observation_rows = []
    evidence_rows = []

    for members in cluster_calls(calls, locus_window):
        chromosome = members[0]["chrom"]
        anchors = [anchor(call) for call in members]
        locus_start = min(int(call["start"]) for call in members)
        locus_end = max(int(call["end"]) for call in members)
        locus_id = stable_id(
            "L", reference_id, chromosome, min(anchors), max(anchors), locus_window
        )
        grouped_alleles: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
        for call in members:
            grouped_alleles[allele_signature(call)].append(call)

        status = "resolved" if len(set(anchors)) == 1 else "provisional"
        if len(grouped_alleles) > 1:
            status = "ambiguous"
        loci_rows.append(
            {
                "locus_id": locus_id,
                "chrom": chromosome,
                "start": locus_start,
                "end": locus_end,
                "anchor_left": min(anchors),
                "anchor_right": max(anchors),
                "status": status,
                "confidence": ".",
                "method": f"breakpoint-span-v0.1;window={locus_window}",
            }
        )

        for signature in sorted(grouped_alleles):
            allele_type, te_name, strand = signature
            allele_calls = grouped_alleles[signature]
            allele_id = stable_id("A", locus_id, *signature)
            component_id = stable_id("C", allele_id, 1, te_name)
            sizes = [parse_number(call["sv_size"]) for call in allele_calls]
            numeric_sizes = [float(value) for value in sizes if value != "."]
            allele_rows.append(
                {
                    "allele_id": allele_id,
                    "locus_id": locus_id,
                    "allele_type": allele_type,
                    "structure": "single",
                    "length": f"{max(numeric_sizes):g}" if numeric_sizes else ".",
                    "sequence_id": ".",
                    "confidence": ".",
                }
            )
            identities = [parse_number(call["identity"]) for call in allele_calls]
            coverages = [parse_number(call["coverage"]) for call in allele_calls]
            component_rows.append(
                {
                    "component_id": component_id,
                    "allele_id": allele_id,
                    "rank": 1,
                    "parent_component_id": ".",
                    "te_name": te_name,
                    "te_family": te_name,
                    "te_class": ".",
                    "start_on_allele": ".",
                    "end_on_allele": ".",
                    "strand": strand,
                    "identity": max((float(v) for v in identities if v != "."), default="."),
                    "coverage": max((float(v) for v in coverages if v != "."), default="."),
                }
            )

            by_sample: dict[str, list[dict[str, str]]] = defaultdict(list)
            for call in allele_calls:
                by_sample[call["sample_id"]].append(call)
                original_call_id = call["te_call"].partition("|")[2] or call["tremolo_id"]
                source = "INSIDER" if "_INSIDER" in call["tremolo_id"] else "OUTSIDER"
                evidence_rows.append(
                    {
                        "call_id": original_call_id,
                        "sample_id": call["sample_id"],
                        "locus_id": locus_id,
                        "allele_id": allele_id,
                        "component_id": component_id,
                        "source": source,
                        "source_file": "TE_INFOS.bed",
                        "assignment": "exact" if len(members) == 1 else "compatible",
                        "reason": "breakpoint_span_and_allele_signature",
                    }
                )
            for sample_id, sample_calls in sorted(by_sample.items()):
                frequencies = [
                    parse_frequency(call["frequency_with_clipped"], call["frequency"])
                    for call in sample_calls
                ]
                numeric_frequencies = [float(v) for v in frequencies if v != "."]
                observation_rows.append(
                    {
                        "sample_id": sample_id,
                        "locus_id": locus_id,
                        "allele_id": allele_id,
                        "frequency": f"{max(numeric_frequencies):.8g}" if numeric_frequencies else ".",
                        "support_reads": ".",
                        "total_reads": ".",
                        "genotype": ".",
                        "quality": ".",
                        "filter": "PASS" if status == "resolved" else "PROVISIONAL_LOCUS",
                    }
                )

    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "loci.tsv", ("locus_id", "chrom", "start", "end", "anchor_left", "anchor_right", "status", "confidence", "method"), loci_rows)
    write_tsv(output / "alleles.tsv", ("allele_id", "locus_id", "allele_type", "structure", "length", "sequence_id", "confidence"), allele_rows)
    write_tsv(output / "components.tsv", ("component_id", "allele_id", "rank", "parent_component_id", "te_name", "te_family", "te_class", "start_on_allele", "end_on_allele", "strand", "identity", "coverage"), component_rows)
    write_tsv(output / "observations.tsv", ("sample_id", "locus_id", "allele_id", "frequency", "support_reads", "total_reads", "genotype", "quality", "filter"), observation_rows)
    write_tsv(output / "evidence.tsv", ("call_id", "sample_id", "locus_id", "allele_id", "component_id", "source", "source_file", "assignment", "reason"), evidence_rows)

    counts = {
        "calls": len(calls), "loci": len(loci_rows), "alleles": len(allele_rows),
        "components": len(component_rows), "observations": len(observation_rows),
        "evidence": len(evidence_rows),
    }
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "coordinate_system": "0-based-half-open",
        "construction_mode": "cohort-breakpoint-span",
        "locus_window": locus_window,
        "reference_id": reference_id,
        "inputs": [
            {"sample_id": run.sample_id, "te_infos": str(run.te_infos.resolve()), "sha256": file_sha256(run.te_infos)}
            for run in runs
        ],
        "counts": counts,
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return counts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", action="append", required=True, type=parse_sample)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--locus-window", type=int, default=20)
    args = parser.parse_args()
    counts = build_cohort(args.sample, args.output, args.reference_id, args.locus_window)
    print("Cohort locus model written:", ", ".join(f"{key}={value}" for key, value in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
