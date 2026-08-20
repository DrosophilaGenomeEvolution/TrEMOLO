#!/usr/bin/env python3
"""Build deterministic data and a self-contained Quarto source for TrEMOLO."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Optional


REPORT_SCHEMA_VERSION = "1.0.0"
TE_INFO_COLUMNS = (
    "chrom",
    "start",
    "end",
    "te_call",
    "strand",
    "tsd",
    "identity",
    "coverage",
    "te_size",
    "new_position",
    "frequency",
    "frequency_with_clipped",
    "sv_size",
    "tremolo_id",
    "event_type",
)
MISSING_VALUES = {"", ".", "NONE", "NONE-FREQA", "NONE-FREQB", "INSIDER"}
PLACEHOLDER_STYLE = "<!-- TREMOLO_REPORT_STYLE -->"
PLACEHOLDER_DATA = "<!-- TREMOLO_REPORT_DATA -->"
PLACEHOLDER_SCRIPT = "<!-- TREMOLO_REPORT_SCRIPT -->"
PLACEHOLDER_TITLE = "%%TREMOLO_REPORT_TITLE%%"
PLACEHOLDER_AUTHOR = "%%TREMOLO_REPORT_AUTHOR%%"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def optional_number(value: str) -> Optional[float]:
    if value in MISSING_VALUES:
        return None
    try:
        number = float(value.replace(",", "."))
    except ValueError:
        return None
    if number != number or number in {float("inf"), float("-inf")}:
        return None
    return number


def display_number(value: Optional[float]) -> Optional[float | int]:
    if value is None:
        return None
    return int(value) if value.is_integer() else value


def classify_source(tremolo_id: str) -> str:
    if "_INSIDER" in tremolo_id:
        return "INSIDER"
    if "_OUTSIDER" in tremolo_id:
        return "OUTSIDER"
    return "UNKNOWN"


def normalize_event_type(event_type: str) -> str:
    upper = event_type.upper()
    if "DEL" in upper or "CONTRACTION" in upper:
        return "deletion"
    if "INS" in upper or "EXPANSION" in upper or "HARD" in upper:
        return "insertion"
    return "other"


def read_te_infos(path: Path) -> list[dict]:
    calls = []
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"empty TE_INFOS file: {path}") from error
        normalized_header = tuple(value.lstrip("#") for value in header)
        if len(header) != len(TE_INFO_COLUMNS):
            raise ValueError(
                f"{path}: expected {len(TE_INFO_COLUMNS)} columns, got {len(header)}"
            )
        if normalized_header[:4] != ("chrom", "start", "end", "TE|ID"):
            raise ValueError(f"{path}: unexpected TE_INFOS header")

        for line_number, fields in enumerate(reader, 2):
            if not fields or not any(fields):
                continue
            if len(fields) != len(TE_INFO_COLUMNS):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(TE_INFO_COLUMNS)} "
                    f"columns, got {len(fields)}"
                )
            row = dict(zip(TE_INFO_COLUMNS, fields))
            try:
                start = int(row["start"])
                end = int(row["end"])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-integer coordinates") from error
            if start < 0 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval {start}-{end}")

            family, separator, event_id = row["te_call"].partition("|")
            if not separator or not family or not event_id:
                raise ValueError(f"{path}:{line_number}: malformed TE|ID value")

            frequency_raw = optional_number(row["frequency"])
            frequency_clipped = optional_number(row["frequency_with_clipped"])
            selected_frequency = (
                frequency_clipped if frequency_clipped is not None else frequency_raw
            )
            new_position = optional_number(row["new_position"])
            call = {
                "row": len(calls) + 1,
                "chrom": row["chrom"],
                "start": start,
                "end": end,
                "anchor": int(new_position) if new_position is not None else start,
                "family": family,
                "event_id": event_id,
                "strand": row["strand"] if row["strand"] in {"+", "-"} else ".",
                "tsd": None if row["tsd"] in MISSING_VALUES else row["tsd"],
                "identity": display_number(optional_number(row["identity"])),
                "coverage": display_number(optional_number(row["coverage"])),
                "te_size": display_number(optional_number(row["te_size"])),
                "new_position": display_number(new_position),
                "frequency": display_number(frequency_raw),
                "frequency_with_clipped": display_number(frequency_clipped),
                "display_frequency": display_number(selected_frequency),
                "sv_size": display_number(optional_number(row["sv_size"])),
                "tremolo_id": row["tremolo_id"],
                "event_type": row["event_type"],
                "event_class": normalize_event_type(row["event_type"]),
                "source": classify_source(row["tremolo_id"]),
            }
            calls.append(call)
    return calls


def read_fasta_index(path: Path) -> list[dict]:
    chromosomes = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FASTA index")
            try:
                length = int(fields[1])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-integer sequence length") from error
            chromosomes.append({"chrom": fields[0], "length": length})
    return chromosomes


def read_manifest(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        return {}
    with path.open() as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: input manifest must contain a JSON object")
    return value


def read_mapping_stats(path: Optional[Path]) -> dict[str, str | int | float]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return {}
    stats = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            key = fields[0].rstrip(":")
            raw_value = fields[1]
            number = optional_number(raw_value)
            stats[key] = display_number(number) if number is not None else raw_value
    return stats


def read_vcf_counts(path: Optional[Path]) -> dict[str, int]:
    if path is None or not path.is_file() or path.stat().st_size == 0:
        return {}
    counts: Counter[str] = Counter()
    pattern = re.compile(r"(?:^|;)SVTYPE=([^;]+)")
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"{path}:{line_number}: malformed VCF record")
            match = pattern.search(fields[7])
            variant_type = match.group(1) if match else fields[4].strip("<>")
            counts[variant_type] += 1
    return dict(sorted(counts.items()))


def mean(values: Iterable[Optional[float | int]]) -> Optional[float]:
    present = [float(value) for value in values if value is not None]
    return round(sum(present) / len(present), 6) if present else None


def summarize_families(calls: list[dict]) -> list[dict]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for call in calls:
        grouped[call["family"]].append(call)
    rows = []
    for family, family_calls in grouped.items():
        source_counts = Counter(call["source"] for call in family_calls)
        rows.append(
            {
                "family": family,
                "calls": len(family_calls),
                "insider": source_counts.get("INSIDER", 0),
                "outsider": source_counts.get("OUTSIDER", 0),
                "tsd_confirmed": sum(call["tsd"] is not None for call in family_calls),
                "mean_frequency": mean(call["display_frequency"] for call in family_calls),
            }
        )
    return sorted(rows, key=lambda row: (-row["calls"], row["family"]))


def proximity_groups(calls: list[dict], window: int) -> list[dict]:
    """Create display-only complete-span groups, never biological locus calls."""
    grouped: dict[str, list[dict]] = defaultdict(list)
    for call in calls:
        grouped[call["chrom"]].append(call)

    raw_groups = []
    for chrom in sorted(grouped):
        ordered = sorted(grouped[chrom], key=lambda call: (call["anchor"], call["row"]))
        current: list[dict] = []
        first_anchor: Optional[int] = None
        for call in ordered:
            anchor = call["anchor"]
            if current and first_anchor is not None and anchor - first_anchor > window:
                if len(current) > 1:
                    raw_groups.append((chrom, current))
                current = []
                first_anchor = None
            if not current:
                first_anchor = anchor
            current.append(call)
        if len(current) > 1:
            raw_groups.append((chrom, current))

    rows = []
    for index, (chrom, group_calls) in enumerate(raw_groups, 1):
        anchors = [call["anchor"] for call in group_calls]
        families = sorted({call["family"] for call in group_calls})
        sources = sorted({call["source"] for call in group_calls})
        rows.append(
            {
                "candidate_group": f"P{index:04d}",
                "chrom": chrom,
                "start": min(anchors),
                "end": max(anchors),
                "span": max(anchors) - min(anchors),
                "calls": len(group_calls),
                "families": families,
                "sources": sources,
                "event_ids": [call["event_id"] for call in group_calls],
            }
        )
    return rows


def build_report_data(
    te_infos: Path,
    genome_index: Path,
    manifest_path: Path,
    mapping_stats_path: Optional[Path],
    sv_vcf_path: Optional[Path],
    enabled_sources: list[str],
    title: str,
    author: str,
    work_directory: str,
    locus_window: int,
) -> dict:
    calls = read_te_infos(te_infos)
    chromosomes = read_fasta_index(genome_index)
    chromosome_lengths = {row["chrom"]: row["length"] for row in chromosomes}
    calls_per_chromosome = Counter(call["chrom"] for call in calls)
    chromosome_summary = [
        {
            "chrom": row["chrom"],
            "length": row["length"],
            "calls": calls_per_chromosome.get(row["chrom"], 0),
        }
        for row in chromosomes
    ]
    for chrom in sorted(set(calls_per_chromosome) - set(chromosome_lengths)):
        chromosome_summary.append(
            {"chrom": chrom, "length": None, "calls": calls_per_chromosome[chrom]}
        )

    sources = Counter(call["source"] for call in calls)
    event_types = Counter(call["event_type"] for call in calls)
    event_classes = Counter(call["event_class"] for call in calls)
    strands = Counter(call["strand"] for call in calls)
    tsd_confirmed = sum(call["tsd"] is not None for call in calls)
    frequencies = [call["display_frequency"] for call in calls]
    frequencies_present = [value for value in frequencies if value is not None]

    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "report": {
            "title": title,
            "author": author,
            "work_directory": work_directory,
            "enabled_sources": enabled_sources,
            "locus_window_bp": locus_window,
            "te_infos_sha256": file_sha256(te_infos),
        },
        "summary": {
            "calls": len(calls),
            "families": len({call["family"] for call in calls}),
            "chromosomes_with_calls": len(calls_per_chromosome),
            "tsd_confirmed": tsd_confirmed,
            "tsd_missing": len(calls) - tsd_confirmed,
            "frequency_available": len(frequencies_present),
            "mean_frequency": mean(frequencies_present),
            "sources": dict(sorted(sources.items())),
            "event_types": dict(sorted(event_types.items())),
            "event_classes": dict(sorted(event_classes.items())),
            "strands": dict(sorted(strands.items())),
        },
        "family_summary": summarize_families(calls),
        "chromosome_summary": chromosome_summary,
        "proximity_groups": proximity_groups(calls, locus_window),
        "mapping_stats": read_mapping_stats(mapping_stats_path),
        "sv_counts": read_vcf_counts(sv_vcf_path),
        "input_manifest": read_manifest(manifest_path),
        "calls": calls,
    }


def atomic_write_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent), text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(content)
        os.replace(str(temporary), str(path))
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def render_source(template: str, style: str, script: str, data: dict) -> str:
    for placeholder in (
        PLACEHOLDER_STYLE,
        PLACEHOLDER_DATA,
        PLACEHOLDER_SCRIPT,
        PLACEHOLDER_TITLE,
        PLACEHOLDER_AUTHOR,
    ):
        if template.count(placeholder) != 1:
            raise ValueError(f"report template must contain exactly one {placeholder}")
    compact_json = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    compact_json = compact_json.replace("<", "\\u003c")
    return (
        template.replace(
            PLACEHOLDER_TITLE,
            json.dumps(data["report"]["title"], ensure_ascii=False),
        )
        .replace(
            PLACEHOLDER_AUTHOR,
            json.dumps(data["report"]["author"], ensure_ascii=False),
        )
        .replace(PLACEHOLDER_STYLE, f"<style>\n{style}\n</style>")
        .replace(
            PLACEHOLDER_DATA,
            '<script id="tremolo-report-data" type="application/json">'
            + compact_json
            + "</script>",
        )
        .replace(PLACEHOLDER_SCRIPT, f"<script>\n{script}\n</script>")
    )


def optional_path(value: str) -> Optional[Path]:
    if not value or value == "/dev/null":
        return None
    return Path(value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--te-infos", type=Path, required=True)
    parser.add_argument("--genome-index", type=Path, required=True)
    parser.add_argument("--input-manifest", type=Path, required=True)
    parser.add_argument("--mapping-stats", default="")
    parser.add_argument("--sv-vcf", default="")
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--style", type=Path, required=True)
    parser.add_argument("--script", type=Path, required=True)
    parser.add_argument("--output-qmd", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--enabled-source", action="append", default=[])
    parser.add_argument("--title", default="TrEMOLO analysis report")
    parser.add_argument("--author", default="")
    parser.add_argument("--work-directory", default="")
    parser.add_argument("--locus-window", type=int, default=100)
    args = parser.parse_args()
    if args.locus_window < 0:
        parser.error("--locus-window must be non-negative")
    return args


def main() -> None:
    args = parse_args()
    data = build_report_data(
        te_infos=args.te_infos,
        genome_index=args.genome_index,
        manifest_path=args.input_manifest,
        mapping_stats_path=optional_path(args.mapping_stats),
        sv_vcf_path=optional_path(args.sv_vcf),
        enabled_sources=args.enabled_source,
        title=args.title,
        author=args.author,
        work_directory=args.work_directory,
        locus_window=args.locus_window,
    )
    json_text = json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    source = render_source(
        args.template.read_text(),
        args.style.read_text(),
        args.script.read_text(),
        data,
    )
    atomic_write_text(args.output_json, json_text)
    atomic_write_text(args.output_qmd, source)
    print(
        "Prepared Quarto report data: {} calls, {} families, {} proximity groups.".format(
            data["summary"]["calls"],
            data["summary"]["families"],
            len(data["proximity_groups"]),
        )
    )


if __name__ == "__main__":
    main()
