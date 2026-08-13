#!/usr/bin/env python3
"""Compare migrated and legacy TrEMOLO variant-calling and TE outputs."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path


INSIDER_FILES = (
    "pm_against_ref.sam.delta",
    "assemblytics_out.Assemblytics.unique_length_filtered_l20000.delta",
    "assemblytics_out.coords.tab",
    "assemblytics_out.coords.csv",
    "assemblytics_out.variants_between_alignments.bed",
    "assemblytics_out.variants_within_alignments.bed",
    "assemblytics_out.Assemblytics_structural_variants.bed",
    "assemblytics_out.Assemblytics_assembly_stats.txt",
)

INSIDER_TE_FILES = (
    "INSERTION_SEQ.fasta",
    "DELETION_SEQ.fasta",
    "INSERTION.bln",
    "DELETION.bln",
    "INSERTION.csv",
    "DELETION.csv",
    "INSERTION_COMBINE_TE.csv",
    "DELETION_COMBINE_TE.csv",
    "INSERTION_TE.bed",
    "DELETION_TE.bed",
    "INSERTION_TE_ON_REF.bed",
    "DELETION_TE_ON_REF.bed",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def normalized_vcf(path: Path) -> list[str]:
    records = []
    with path.open() as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            fields = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                fields[-1] = "<sample>"
            records.append("\t".join(fields))
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("legacy", type=Path)
    parser.add_argument("migrated", type=Path)
    args = parser.parse_args()
    errors = []
    for relative in INSIDER_FILES:
        old = args.legacy / "INSIDER/VARIANT_CALLING" / relative
        new = args.migrated / "INSIDER/VARIANT_CALLING" / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing INSIDER output: {relative}")
        elif digest(old) != digest(new):
            errors.append(f"INSIDER output differs: {relative}")

    for relative in INSIDER_TE_FILES:
        old = args.legacy / "INSIDER/TE_DETECTION" / relative
        new = args.migrated / "INSIDER/TE_DETECTION" / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing INSIDER TE output: {relative}")
        elif digest(old) != digest(new):
            errors.append(f"INSIDER TE output differs: {relative}")

    position_file = Path("POSITION_TE_INSIDER.bed")
    old_position = args.legacy / position_file
    new_position = args.migrated / position_file
    if not old_position.is_file() or not new_position.is_file():
        errors.append(f"missing INSIDER TE output: {position_file}")
    elif digest(old_position) != digest(new_position):
        errors.append(f"INSIDER TE output differs: {position_file}")

    old_vcf = args.legacy / "OUTSIDER/VARIANT_CALLING/SV.vcf"
    new_vcf = args.migrated / "OUTSIDER/VARIANT_CALLING/SV.vcf"
    if not old_vcf.is_file() or not new_vcf.is_file():
        errors.append("missing OUTSIDER VCF")
    elif normalized_vcf(old_vcf) != normalized_vcf(new_vcf):
        errors.append("OUTSIDER VCF records differ")

    if errors:
        print("INSIDER/OUTSIDER migration check: FAILED")
        for error in errors:
            print(f"- {error}")
        return 1
    print("INSIDER/OUTSIDER migration check: PASSED")
    print(f"Exact INSIDER files: {len(INSIDER_FILES)}")
    print(f"Exact INSIDER TE files: {len(INSIDER_TE_FILES) + 1}")
    print("OUTSIDER VCF: identical after volatile header normalization")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
