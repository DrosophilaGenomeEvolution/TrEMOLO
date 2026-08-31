#!/usr/bin/env python3
"""Build a locus-aware population trajectory table and standalone Quarto source."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional


PIPELINE_ROOT = Path(__file__).resolve().parents[3]
LOCUS_LIB = PIPELINE_ROOT / "lib/python/locus"
sys.path.insert(0, str(LOCUS_LIB))
from build_cohort_model import build_cohort  # noqa: E402
from build_locus_model import InputRun  # noqa: E402


SCHEMA_VERSION = "1.0.0"
REQUIRED_MANIFEST_COLUMNS = ("sample_id", "timepoint", "replicate", "te_infos")
OBSERVATION_COLUMNS = (
    "sample_id", "timepoint", "replicate", "locus_id", "allele_id", "chrom",
    "start", "end", "te_family", "strand", "frequency", "observation_status",
    "support_reads", "total_reads", "filter",
)
TRAJECTORY_COLUMNS = (
    "locus_id", "allele_id", "chrom", "start", "end", "te_family", "strand",
    "observed_samples", "unquantified_samples", "missing_samples",
    "observed_timepoints", "slope_per_timepoint", "r_squared", "trend",
)


@dataclass(frozen=True)
class SampleSpec:
    sample_id: str
    timepoint: float
    timepoint_label: str
    replicate: str
    te_infos: Path


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def parse_timepoint(value: str, context: str) -> tuple[float, str]:
    label = value.strip()
    numeric = label[1:] if label.upper().startswith("G") else label
    try:
        timepoint = float(numeric)
    except ValueError as error:
        raise ValueError(f"{context}: timepoint must be numeric or use G<number>") from error
    if not math.isfinite(timepoint):
        raise ValueError(f"{context}: timepoint must be finite")
    return timepoint, label


def resolve_input_path(raw_path: str, base_directory: Path) -> Path:
    path = Path(raw_path).expanduser()
    if not path.is_absolute():
        path = base_directory / path
    if path.is_dir():
        path = path / "TE_INFOS.bed"
    return path.resolve()


def read_tabular_manifest(path: Path) -> list[SampleSpec]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or any(
            column not in reader.fieldnames for column in REQUIRED_MANIFEST_COLUMNS
        ):
            raise ValueError(
                f"{path}: manifest requires columns "
                + ", ".join(REQUIRED_MANIFEST_COLUMNS)
            )
        samples = []
        for line_number, row in enumerate(reader, 2):
            sample_id = row["sample_id"].strip()
            replicate = row["replicate"].strip() or "1"
            if not sample_id:
                raise ValueError(f"{path}:{line_number}: empty sample_id")
            timepoint, label = parse_timepoint(
                row["timepoint"], f"{path}:{line_number}"
            )
            samples.append(
                SampleSpec(
                    sample_id=sample_id,
                    timepoint=timepoint,
                    timepoint_label=label,
                    replicate=replicate,
                    te_infos=resolve_input_path(row["te_infos"].strip(), path.parent),
                )
            )
    return samples


def read_legacy_manifest(path: Path) -> list[SampleSpec]:
    samples = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            value = raw_line.strip()
            if not value or value.startswith("#"):
                continue
            raw_path = value
            label = f"G{len(samples) + 1}"
            prefix, separator, suffix = value.rpartition(":G")
            if separator and suffix.isdigit():
                raw_path = prefix
                label = f"G{suffix}"
            timepoint, label = parse_timepoint(label, f"{path}:{line_number}")
            samples.append(
                SampleSpec(
                    sample_id=label,
                    timepoint=timepoint,
                    timepoint_label=label,
                    replicate="1",
                    te_infos=resolve_input_path(raw_path.rstrip("/"), Path.cwd()),
                )
            )
    return samples


def read_sample_manifest(path: Path) -> list[SampleSpec]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with path.open() as handle:
        first_content = next(
            (line.rstrip("\r\n") for line in handle if line.strip()), ""
        )
    if first_content.split("\t")[:2] == ["sample_id", "timepoint"]:
        samples = read_tabular_manifest(path)
    else:
        samples = read_legacy_manifest(path)
    if len(samples) < 2:
        raise ValueError("population analysis requires at least two samples")
    seen = set()
    for sample in samples:
        if sample.sample_id in seen:
            raise ValueError(f"duplicate sample_id: {sample.sample_id}")
        seen.add(sample.sample_id)
        if not sample.te_infos.is_file():
            raise FileNotFoundError(sample.te_infos)
    return sorted(samples, key=lambda item: (item.timepoint, item.sample_id))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta_lengths(path: Optional[Path]) -> dict[str, int]:
    if path is None:
        return {}
    lengths: dict[str, int] = {}
    name: Optional[str] = None
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
                if not name or name in lengths:
                    raise ValueError(f"{path}:{line_number}: invalid duplicate FASTA id")
                lengths[name] = 0
            elif name is None:
                raise ValueError(f"{path}:{line_number}: sequence before FASTA header")
            else:
                lengths[name] += len(line)
    return lengths


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


def regression(values: list[tuple[float, float]]) -> tuple[Optional[float], Optional[float]]:
    if len(values) < 2 or len({value[0] for value in values}) < 2:
        return None, None
    mean_x = sum(value[0] for value in values) / len(values)
    mean_y = sum(value[1] for value in values) / len(values)
    denominator = sum((value[0] - mean_x) ** 2 for value in values)
    slope = sum(
        (value[0] - mean_x) * (value[1] - mean_y) for value in values
    ) / denominator
    predictions = [mean_y + slope * (value[0] - mean_x) for value in values]
    total = sum((value[1] - mean_y) ** 2 for value in values)
    residual = sum(
        (value[1] - prediction) ** 2
        for value, prediction in zip(values, predictions)
    )
    r_squared = 1.0 if total == 0 else max(0.0, 1.0 - residual / total)
    return slope, r_squared


def classify_trend(
    timepoint_values: list[tuple[float, float]], slope: Optional[float], epsilon: float
) -> str:
    if slope is None:
        return "insufficient_data"
    ordered = sorted(timepoint_values)
    differences = [
        ordered[index + 1][1] - ordered[index][1]
        for index in range(len(ordered) - 1)
    ]
    increasing = any(value > epsilon for value in differences)
    decreasing = any(value < -epsilon for value in differences)
    if increasing and decreasing:
        return "variable"
    if slope > epsilon:
        return "increasing"
    if slope < -epsilon:
        return "decreasing"
    return "stable"


def build_population_tables(
    samples: list[SampleSpec], locus_directory: Path, trend_epsilon: float
) -> tuple[list[dict], list[dict]]:
    loci = {row["locus_id"]: row for row in read_tsv(locus_directory / "loci.tsv")}
    alleles = {row["allele_id"]: row for row in read_tsv(locus_directory / "alleles.tsv")}
    components_by_allele: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(locus_directory / "components.tsv"):
        components_by_allele[row["allele_id"]].append(row)
    observations = {
        (row["sample_id"], row["allele_id"]): row
        for row in read_tsv(locus_directory / "observations.tsv")
    }
    sample_by_id = {sample.sample_id: sample for sample in samples}

    observation_rows = []
    trajectory_rows = []
    for allele_id, allele in sorted(
        alleles.items(),
        key=lambda item: (
            loci[item[1]["locus_id"]]["chrom"],
            int(loci[item[1]["locus_id"]]["start"]),
            item[0],
        ),
    ):
        locus = loci[allele["locus_id"]]
        components = sorted(
            components_by_allele[allele_id], key=lambda row: int(row["rank"])
        )
        te_family = "+".join(row["te_family"] for row in components) or "."
        strand = components[0]["strand"] if len(components) == 1 else "."
        observed = 0
        unquantified = 0
        missing = 0
        by_timepoint: dict[float, list[float]] = defaultdict(list)
        for sample in samples:
            source = observations.get((sample.sample_id, allele_id))
            frequency = "."
            if source is None:
                status = "missing"
                missing += 1
                support_reads = total_reads = filter_value = "."
            else:
                frequency = source["frequency"]
                support_reads = source["support_reads"]
                total_reads = source["total_reads"]
                filter_value = source["filter"]
                if frequency == ".":
                    status = "observed_unquantified"
                    unquantified += 1
                else:
                    status = "observed"
                    observed += 1
                    by_timepoint[sample.timepoint].append(float(frequency))
            observation_rows.append(
                {
                    "sample_id": sample.sample_id,
                    "timepoint": sample.timepoint_label,
                    "replicate": sample.replicate,
                    "locus_id": locus["locus_id"],
                    "allele_id": allele_id,
                    "chrom": locus["chrom"],
                    "start": locus["start"],
                    "end": locus["end"],
                    "te_family": te_family,
                    "strand": strand,
                    "frequency": frequency,
                    "observation_status": status,
                    "support_reads": support_reads,
                    "total_reads": total_reads,
                    "filter": filter_value,
                }
            )
        timepoint_values = [
            (timepoint, sum(values) / len(values))
            for timepoint, values in sorted(by_timepoint.items())
        ]
        slope, r_squared = regression(timepoint_values)
        trajectory_rows.append(
            {
                "locus_id": locus["locus_id"],
                "allele_id": allele_id,
                "chrom": locus["chrom"],
                "start": locus["start"],
                "end": locus["end"],
                "te_family": te_family,
                "strand": strand,
                "observed_samples": observed,
                "unquantified_samples": unquantified,
                "missing_samples": missing,
                "observed_timepoints": len(timepoint_values),
                "slope_per_timepoint": "." if slope is None else f"{slope:.8g}",
                "r_squared": "." if r_squared is None else f"{r_squared:.8g}",
                "trend": classify_trend(timepoint_values, slope, trend_epsilon),
            }
        )
    if set(sample_by_id) != {sample.sample_id for sample in samples}:
        raise AssertionError("sample metadata mismatch")
    return observation_rows, trajectory_rows


def escape_embedded_json(data: dict) -> str:
    return json.dumps(data, ensure_ascii=False, separators=(",", ":")).replace(
        "<", "\\u003c"
    )


def render_qmd(template: str, style: str, script: str, data: dict) -> str:
    replacements = {
        "<!-- POPULATION_REPORT_STYLE -->": f"<style>\n{style}\n</style>",
        "<!-- POPULATION_REPORT_DATA -->": (
            '<script id="population-report-data" type="application/json">'
            + escape_embedded_json(data)
            + "</script>"
        ),
        "<!-- POPULATION_REPORT_SCRIPT -->": f"<script>\n{script}\n</script>",
        '"POPULATION_REPORT_TITLE"': json.dumps(data["report"]["title"]),
    }
    for placeholder in replacements:
        if template.count(placeholder) != 1:
            raise ValueError(f"template requires exactly one {placeholder}")
    for placeholder, value in replacements.items():
        template = template.replace(placeholder, value)
    return template


def build(
    manifest: Path,
    output: Path,
    reference_id: str,
    locus_window: int,
    trend_epsilon: float,
    template: Path,
    style: Path,
    script: Path,
    title: str,
    genome: Optional[Path] = None,
) -> dict:
    if locus_window < 0:
        raise ValueError("locus window must be non-negative")
    if trend_epsilon < 0:
        raise ValueError("trend epsilon must be non-negative")
    samples = read_sample_manifest(manifest)
    chromosome_lengths = read_fasta_lengths(genome)
    if genome is not None:
        expected_reference = f"sha256:{file_sha256(genome)}"
        if reference_id != expected_reference:
            raise ValueError(
                "reference_id must equal sha256:<genome SHA-256> when --genome is used"
            )
    output.mkdir(parents=True, exist_ok=True)
    locus_directory = output / "LOCUS_MODEL"
    counts = build_cohort(
        [InputRun(sample.sample_id, sample.te_infos) for sample in samples],
        locus_directory,
        reference_id,
        locus_window,
    )
    observation_rows, trajectory_rows = build_population_tables(
        samples, locus_directory, trend_epsilon
    )
    write_tsv(output / "population-observations.tsv", OBSERVATION_COLUMNS, observation_rows)
    write_tsv(output / "population-trajectories.tsv", TRAJECTORY_COLUMNS, trajectory_rows)
    status_counts: dict[str, int] = defaultdict(int)
    for row in observation_rows:
        status_counts[row["observation_status"]] += 1
    trend_counts: dict[str, int] = defaultdict(int)
    for row in trajectory_rows:
        trend_counts[row["trend"]] += 1
    data = {
        "schema_version": SCHEMA_VERSION,
        "report": {
            "title": title,
            "reference_id": reference_id,
            "locus_window": locus_window,
            "trend_epsilon": trend_epsilon,
            "manifest": str(manifest.resolve()),
            "manifest_sha256": file_sha256(manifest),
            "genome": str(genome.resolve()) if genome else None,
            "genome_sha256": file_sha256(genome) if genome else None,
        },
        "summary": {
            **counts,
            "samples": len(samples),
            "trajectories": len(trajectory_rows),
            "observation_status": dict(sorted(status_counts.items())),
            "trends": dict(sorted(trend_counts.items())),
        },
        "samples": [
            {
                "sample_id": sample.sample_id,
                "timepoint": sample.timepoint,
                "timepoint_label": sample.timepoint_label,
                "replicate": sample.replicate,
                "te_infos": str(sample.te_infos),
                "sha256": file_sha256(sample.te_infos),
            }
            for sample in samples
        ],
        "chromosome_lengths": chromosome_lengths,
        "observations": observation_rows,
        "trajectories": trajectory_rows,
    }
    atomic_write_text(
        output / "report-data.json",
        json.dumps(data, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
    )
    qmd = render_qmd(
        template.read_text(), style.read_text(), script.read_text(), data
    )
    atomic_write_text(output / "report.qmd", qmd)
    return data


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--genome", type=Path)
    parser.add_argument("--locus-window", type=int, default=20)
    parser.add_argument("--trend-epsilon", type=float, default=0.001)
    parser.add_argument(
        "--template",
        type=Path,
        default=PIPELINE_ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.qmd",
    )
    parser.add_argument(
        "--style",
        type=Path,
        default=PIPELINE_ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.css",
    )
    parser.add_argument(
        "--script",
        type=Path,
        default=PIPELINE_ROOT / "modules/1-FREQUENCY-MULTI-GENERATIONS/report.js",
    )
    parser.add_argument("--title", default="TrEMOLO population trajectories")
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> int:
    args = parse_args(argv)
    data = build(
        manifest=args.input,
        output=args.output,
        reference_id=args.reference_id,
        locus_window=args.locus_window,
        trend_epsilon=args.trend_epsilon,
        template=args.template,
        style=args.style,
        script=args.script,
        title=args.title,
        genome=args.genome,
    )
    print(
        "Population report prepared: "
        f"{data['summary']['samples']} samples, "
        f"{data['summary']['loci']} loci, "
        f"{data['summary']['trajectories']} trajectories."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
