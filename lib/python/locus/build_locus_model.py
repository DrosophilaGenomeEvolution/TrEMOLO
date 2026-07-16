#!/usr/bin/env python3
"""Build the normalized TrEMOLO locus model from legacy TE_INFOS.bed calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


SCHEMA_VERSION = "0.1.0"
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


@dataclass(frozen=True)
class InputRun:
    sample_id: str
    te_infos: Path


def stable_id(prefix: str, *parts: object) -> str:
    payload = "\x1f".join(map(str, parts)).encode()
    return f"{prefix}{hashlib.sha256(payload).hexdigest()[:16]}"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_number(value: str) -> str:
    if value in {"", ".", "NONE", "INSIDER"}:
        return "."
    try:
        return str(float(value))
    except ValueError:
        return "."


def parse_frequency(*values: str) -> str:
    for value in values:
        number = parse_number(value)
        if number != ".":
            frequency = float(number) / 100.0
            if 0.0 <= frequency <= 1.0:
                return f"{frequency:.8g}"
    return "."


def parse_sample(value: str) -> InputRun:
    if "=" not in value:
        raise argparse.ArgumentTypeError("sample must use SAMPLE_ID=/path/TE_INFOS.bed")
    sample_id, path = value.split("=", 1)
    if not sample_id or not path:
        raise argparse.ArgumentTypeError("sample id and TE_INFOS.bed path are required")
    return InputRun(sample_id, Path(path))


def read_calls(run: InputRun) -> Iterable[dict[str, str]]:
    with run.te_infos.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"empty TE_INFOS file: {run.te_infos}") from error
        if len(header) != len(TE_INFO_COLUMNS):
            raise ValueError(
                f"{run.te_infos}: expected {len(TE_INFO_COLUMNS)} columns, got {len(header)}"
            )
        for line_number, fields in enumerate(reader, 2):
            if not fields or not any(fields):
                continue
            if len(fields) != len(TE_INFO_COLUMNS):
                raise ValueError(
                    f"{run.te_infos}:{line_number}: expected {len(TE_INFO_COLUMNS)} "
                    f"columns, got {len(fields)}"
                )
            row = dict(zip(TE_INFO_COLUMNS, fields))
            row["sample_id"] = run.sample_id
            row["line_number"] = str(line_number)
            yield row


def write_tsv(path: Path, columns: tuple[str, ...], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def build(runs: list[InputRun], output: Path, reference_id: str) -> dict[str, int]:
    output.mkdir(parents=True, exist_ok=True)
    loci = []
    alleles = []
    components = []
    observations = []
    evidence = []
    seen_calls: set[tuple[str, str]] = set()

    for run in runs:
        if not run.te_infos.is_file():
            raise FileNotFoundError(run.te_infos)
        for call in read_calls(run):
            call_key = (run.sample_id, call["tremolo_id"])
            if call_key in seen_calls:
                raise ValueError(f"duplicate call id for sample {run.sample_id}: {call['tremolo_id']}")
            seen_calls.add(call_key)

            te_name, separator, original_call_id = call["te_call"].partition("|")
            if not separator:
                original_call_id = call["tremolo_id"]
            locus_id = stable_id(
                "L", reference_id, call["chrom"], call["start"], call["end"], call["tremolo_id"]
            )
            allele_id = stable_id("A", locus_id, call["event_type"], te_name)
            component_id = stable_id("C", allele_id, 1, te_name)
            event_upper = call["event_type"].upper()
            if "DEL" in event_upper or "CONTRACTION" in event_upper:
                allele_type = "deletion"
            else:
                allele_type = "insertion"

            loci.append(
                {
                    "locus_id": locus_id,
                    "chrom": call["chrom"],
                    "start": call["start"],
                    "end": call["end"],
                    "anchor_left": call["new_position"] if call["new_position"] != "NONE" else ".",
                    "anchor_right": ".",
                    "status": "provisional",
                    "confidence": ".",
                    "method": "one-call-per-locus-v0.1",
                }
            )
            alleles.append(
                {
                    "allele_id": allele_id,
                    "locus_id": locus_id,
                    "allele_type": allele_type,
                    "structure": "single",
                    "length": parse_number(call["sv_size"]),
                    "sequence_id": ".",
                    "confidence": ".",
                }
            )
            components.append(
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
                    "strand": call["strand"] if call["strand"] in {"+", "-"} else ".",
                    "identity": parse_number(call["identity"]),
                    "coverage": parse_number(call["coverage"]),
                }
            )
            observations.append(
                {
                    "sample_id": run.sample_id,
                    "locus_id": locus_id,
                    "allele_id": allele_id,
                    "frequency": parse_frequency(
                        call["frequency_with_clipped"], call["frequency"]
                    ),
                    "support_reads": ".",
                    "total_reads": ".",
                    "genotype": ".",
                    "quality": ".",
                    "filter": "PASS",
                }
            )
            source = "INSIDER" if "_INSIDER" in call["tremolo_id"] else "OUTSIDER"
            evidence.append(
                {
                    "call_id": original_call_id,
                    "sample_id": run.sample_id,
                    "locus_id": locus_id,
                    "allele_id": allele_id,
                    "component_id": component_id,
                    "source": source,
                    "source_file": run.te_infos.name,
                    "assignment": "exact",
                    "reason": "legacy_call_projection",
                }
            )

    write_tsv(
        output / "loci.tsv",
        ("locus_id", "chrom", "start", "end", "anchor_left", "anchor_right", "status", "confidence", "method"),
        loci,
    )
    write_tsv(
        output / "alleles.tsv",
        ("allele_id", "locus_id", "allele_type", "structure", "length", "sequence_id", "confidence"),
        alleles,
    )
    write_tsv(
        output / "components.tsv",
        ("component_id", "allele_id", "rank", "parent_component_id", "te_name", "te_family", "te_class", "start_on_allele", "end_on_allele", "strand", "identity", "coverage"),
        components,
    )
    write_tsv(
        output / "observations.tsv",
        ("sample_id", "locus_id", "allele_id", "frequency", "support_reads", "total_reads", "genotype", "quality", "filter"),
        observations,
    )
    write_tsv(
        output / "evidence.tsv",
        ("call_id", "sample_id", "locus_id", "allele_id", "component_id", "source", "source_file", "assignment", "reason"),
        evidence,
    )

    manifest = {
        "schema_version": SCHEMA_VERSION,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "coordinate_system": "0-based-half-open",
        "construction_mode": "one-call-per-locus",
        "reference_id": reference_id,
        "inputs": [
            {
                "sample_id": run.sample_id,
                "te_infos": str(run.te_infos.resolve()),
                "sha256": file_sha256(run.te_infos),
            }
            for run in runs
        ],
        "counts": {
            "loci": len(loci),
            "alleles": len(alleles),
            "components": len(components),
            "observations": len(observations),
            "evidence": len(evidence),
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest["counts"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--sample",
        action="append",
        required=True,
        type=parse_sample,
        help="SAMPLE_ID=/path/to/TE_INFOS.bed; repeat for multiple inputs",
    )
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--reference-id",
        required=True,
        help="stable reference name or checksum shared by all inputs",
    )
    args = parser.parse_args()
    counts = build(args.sample, args.output, args.reference_id)
    print("Locus model written:", ", ".join(f"{key}={value}" for key, value in counts.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
