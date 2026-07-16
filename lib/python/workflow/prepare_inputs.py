#!/usr/bin/env python3
"""Prepare immutable TrEMOLO workflow inputs without shell side effects."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path


FORBIDDEN_HEADER_CHARACTERS = set(":/\\|[]")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fasta_records(path: Path):
    header = None
    sequence = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, 1):
            value = line.strip()
            if not value:
                continue
            if value.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence).upper()
                header = value[1:]
                sequence = []
            elif header is None:
                raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
            else:
                sequence.append(value)
    if header is not None:
        yield header, "".join(sequence).upper()


def relative_symlink(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink() or destination.exists():
        destination.unlink()
    destination.symlink_to(os.path.relpath(source.resolve(), destination.parent.resolve()))


def prepare(
    genome: Path,
    te_database: Path,
    output: Path,
    reference: Path | None = None,
    sample: Path | None = None,
) -> dict:
    inputs = {"genome": genome, "te_database": te_database}
    if reference is not None:
        inputs["reference"] = reference
    if sample is not None:
        inputs["sample"] = sample
    for label, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(f"{label}: {path}")

    output.mkdir(parents=True, exist_ok=True)
    input_directory = output / "INPUT"
    input_directory.mkdir(exist_ok=True)
    relative_symlink(genome, input_directory / "genome.fasta")
    if reference is not None:
        relative_symlink(reference, input_directory / "reference.fasta")
    if sample is not None:
        relative_symlink(sample, input_directory / "reads.fastq")

    records = list(fasta_records(te_database))
    if not records:
        raise ValueError(f"no FASTA record found in {te_database}")
    needs_pseudonyms = any(
        FORBIDDEN_HEADER_CHARACTERS.intersection(header) for header, _ in records
    )
    prepared_database = input_directory / "te_database.fasta"
    mapping = input_directory / "te_headers.tsv"
    with prepared_database.open("w") as fasta, mapping.open("w") as names:
        names.write("original\tprepared\n")
        for index, (header, sequence) in enumerate(records, 1):
            prepared = f"TrEMOLOTE{index}" if needs_pseudonyms else header
            fasta.write(f">{prepared}\n{sequence}\n")
            names.write(f"{header}\t{prepared}\n")

    manifest = {
        "format_version": 1,
        "pseudonyms_used": needs_pseudonyms,
        "te_records": len(records),
        "inputs": {
            label: {"path": str(path.resolve()), "sha256": sha256(path)}
            for label, path in sorted(inputs.items())
        },
        "prepared": {
            "genome": "INPUT/genome.fasta",
            "reference": "INPUT/reference.fasta" if reference is not None else None,
            "sample": "INPUT/reads.fastq" if sample is not None else None,
            "te_database": "INPUT/te_database.fasta",
            "te_headers": "INPUT/te_headers.tsv",
        },
    }
    (output / "input_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (output / ".inputs.prepared").write_text("prepared\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--te-database", required=True, type=Path)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--sample", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    manifest = prepare(args.genome, args.te_database, args.output, args.reference, args.sample)
    print(
        f"Prepared {manifest['te_records']} TE records "
        f"(pseudonyms={manifest['pseudonyms_used']})"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
