#!/usr/bin/env python3
"""Compare migrated and legacy TrEMOLO variant-calling and TE outputs."""

from __future__ import annotations

import argparse
from collections import Counter
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

OUTSIDER_TE_EXACT_FILES = (
    "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS.bed",
    "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bed",
    "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.fasta",
    "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln",
    "OUTSIDER/TrEMOLO_SV_TE/INS/RD_NUMBER.txt",
    "OUTSIDER/TrEMOLO_SV_TE/INS/SV_SIZE.tsv",
    "OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.csv",
    "OUTSIDER/TrEMOLO_SV_TE/INS/COMBINE_INS_TREMOLO.csv",
    "OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO_COUNT.csv",
    "OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.bed",
    "OUTSIDER/TrEMOLO_SV_TE/INS/ID.txt",
    "OUTSIDER/TrEMOLO_SV_TE/INS/COUNT_TE_IN_RS.txt",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.vcf",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.vcf.bis",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.bed",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.fasta",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.fasta.bis",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SV_SOFT.bln",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SOFT_TE.csv",
    "OUTSIDER/TrEMOLO_SV_TE/SOFT/SOFT_TE.bed",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/SV_HARD.tr_vcf",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/ID.txt",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/RS_HARD.fastq",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/RS_HARD.fasta",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/RS_HARD.fasta.fa_trml_idx",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD.bed",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD.fasta",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD.bln",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.csv",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/COMBINE_SV_HARD.csv",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE_COUNT.csv",
    "OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.bed",
    "OUTSIDER/VARIANT_CALLING/SEQUENCE_INDEL.fasta",
    "OUTSIDER/TE_DETECTION/TE_VR.bed",
    "OUTSIDER/TE_DETECTION/RD_NUMBER.txt",
    "OUTSIDER/TE_DETECTION/BLAST_SEQUENCE_INDEL_vs_DBTE.bln",
    "OUTSIDER/TE_DETECTION/FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv",
    "OUTSIDER/TE_DETECTION/COMBINE_TE.csv",
    "OUTSIDER/TE_DETECTION/FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE_COUNT.csv",
    "OUTSIDER/TE_DETECTION/POSITION_START_TE.bed",
    "OUTSIDER/TE_DETECTION/ID.txt",
    "OUTSIDER/TE_DETECTION/COUNT_TE_IN_RS.txt",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TrEMOLO_TE_OUTSIDER_CLUSTER_100.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_ID_TrEMOLO.txt",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_TrEMOLO_NOT_FOUND_IN_sniffles.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_SOFT_CLUSTER_100.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_SOFT_NOT_FOUND_IN_sniffles_and_TrEMOLO.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_HARD_CLUSTER_100.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_ID_HARD_commun.txt",
    "OUTSIDER/TE_DETECTION/MERGE_TE/tmp_TE_HARD_NOT_FOUND_IN_SOFT_sniffles_and_TrEMOLO.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL2.bed",
    "OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL_COUNT.csv",
    "POSITION_TE_OUTSIDER.bed",
)

OUTSIDER_CIGAR_CANDIDATE_RECOVERY_FILES = (
    "OUTSIDER/TrEMOLO_SV_TE/INS/INS_FOR_TSD.txt",
)

OUTSIDER_TSD_EXACT_FILES = (
    "OUTSIDER/SIZE_SEQ.tsv",
    "OUTSIDER/FK/ALL_FK_REPORT_FT_READS.bed",
    "OUTSIDER/FK/ALL_FK_REPORT_FT_GENOME.bed",
    "OUTSIDER/FK/ALL_FK_REPORT_FT.bed",
    "OUTSIDER/TSD/TSD_TE.tsv",
)

OUTSIDER_TSD_EXPECTED_COUNTS = {
    "OUTSIDER/FK/ALL_FK_REPORT_FT_READS.bed": 18,
    "OUTSIDER/FK/ALL_FK_REPORT_FT_GENOME.bed": 18,
    "OUTSIDER/FK/ALL_FK_REPORT_FT.bed": 18,
    "OUTSIDER/TSD/TSD_TE.tsv": 9,
}

OUTSIDER_TSD_EXACT_DIRECTORIES = (
    "OUTSIDER/ET_FIND_FA",
    "OUTSIDER/FK/READS",
    "OUTSIDER/FK/GENOME",
)

OUTSIDER_FLANK_ELIGIBILITY_EXPECTED = Counter(
    {
        ("accepted", "accepted"): 18,
        ("rejected", "insufficient_left_flank"): 5,
        ("rejected", "insufficient_right_flank"): 3,
    }
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


def line_counts(path: Path) -> Counter[bytes]:
    return Counter(path.read_bytes().splitlines())


def insertion_keys(path: Path, tsd: bool) -> list[tuple[bytes, ...]]:
    """Project CIGAR and raw flank records onto shared insertion fields."""
    records = []
    for line in path.read_bytes().splitlines():
        fields = line.split(b"\t")
        indexes = (0, 1, 2, 4, 5, 6, 7) if tsd else (0, 1, 3, 4, 5, 6, 7)
        records.append(tuple(fields[index] for index in indexes))
    return sorted(records)


def directory_files(path: Path) -> dict[Path, str]:
    """Return stable hashes for generated directory contents."""

    return {
        item.relative_to(path): digest(item)
        for item in sorted(path.rglob("*"))
        if item.is_file() and item.name != ".snakemake_timestamp"
    }


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

    for relative in OUTSIDER_TE_EXACT_FILES:
        old = args.legacy / relative
        new = args.migrated / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing OUTSIDER TE output: {relative}")
        elif digest(old) != digest(new):
            errors.append(f"OUTSIDER TE output differs: {relative}")

    recovered_cigar_candidates = 0
    for relative in OUTSIDER_CIGAR_CANDIDATE_RECOVERY_FILES:
        old = args.legacy / relative
        new = args.migrated / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing OUTSIDER CIGAR candidate output: {relative}")
        else:
            old_records = line_counts(old)
            new_records = line_counts(new)
            if old_records - new_records:
                errors.append(f"legacy OUTSIDER CIGAR candidates missing: {relative}")
            recovered_cigar_candidates += sum((new_records - old_records).values())

    migrated_insertions = (
        args.migrated / "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS.bed"
    )
    migrated_cigar_candidates = (
        args.migrated / OUTSIDER_CIGAR_CANDIDATE_RECOVERY_FILES[0]
    )
    if migrated_insertions.is_file() and migrated_cigar_candidates.is_file():
        if insertion_keys(migrated_insertions, tsd=False) != insertion_keys(
            migrated_cigar_candidates, tsd=True
        ):
            errors.append(
                "OUTSIDER raw flank candidates do not cover every CIGAR insertion"
            )

    for relative in OUTSIDER_TSD_EXACT_FILES:
        old = args.legacy / relative
        new = args.migrated / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing OUTSIDER flank/TSD output: {relative}")
        elif digest(old) != digest(new):
            errors.append(f"OUTSIDER flank/TSD output differs: {relative}")

    for relative, expected_count in OUTSIDER_TSD_EXPECTED_COUNTS.items():
        path = args.migrated / relative
        if path.is_file():
            actual_count = sum(1 for _ in path.open("rb"))
            if actual_count != expected_count:
                errors.append(
                    f"OUTSIDER flank/TSD line count differs: {relative} "
                    f"({actual_count} != {expected_count})"
                )

    for relative in OUTSIDER_TSD_EXACT_DIRECTORIES:
        old = args.legacy / relative
        new = args.migrated / relative
        if not old.is_dir() or not new.is_dir():
            errors.append(f"missing OUTSIDER flank artifact directory: {relative}")
        elif directory_files(old) != directory_files(new):
            errors.append(f"OUTSIDER flank artifact directory differs: {relative}")

    validation = args.migrated / "OUTSIDER/FK/TSD_FLANK_ELIGIBILITY.tsv"
    if not validation.is_file():
        errors.append("missing OUTSIDER TSD flank eligibility report")
    else:
        rows = validation.read_text().splitlines()
        if not rows:
            errors.append("invalid OUTSIDER TSD flank eligibility report")
        else:
            header = rows[0].split("\t")
            required_columns = {"candidate_index", "event_id", "status", "reason"}
            if not required_columns.issubset(header):
                errors.append("invalid OUTSIDER TSD flank eligibility report")
            else:
                candidate_index = header.index("candidate_index")
                event_index = header.index("event_id")
                status_index = header.index("status")
                reason_index = header.index("reason")
                observed = Counter()
                observed_indexes = []
                accepted_events = Counter()
                for row in rows[1:]:
                    fields = row.split("\t")
                    required_index = max(
                        candidate_index, event_index, status_index, reason_index
                    )
                    if len(fields) <= required_index:
                        errors.append("malformed OUTSIDER TSD flank eligibility row")
                        continue
                    observed[(fields[status_index], fields[reason_index])] += 1
                    observed_indexes.append(fields[candidate_index])
                    if fields[status_index] == "accepted":
                        accepted_events[fields[event_index]] += 1
                if observed != OUTSIDER_FLANK_ELIGIBILITY_EXPECTED:
                    errors.append(
                        "OUTSIDER TSD flank eligibility summary differs: "
                        f"{dict(observed)} != "
                        f"{dict(OUTSIDER_FLANK_ELIGIBILITY_EXPECTED)}"
                    )
                expected_indexes = [str(index) for index in range(1, 27)]
                if observed_indexes != expected_indexes:
                    errors.append("OUTSIDER TSD candidate indexes are not 1..26")
                combined = args.migrated / "OUTSIDER/FK/ALL_FK_REPORT_FT.bed"
                if combined.is_file():
                    flank_events = Counter(
                        line.split("\t")[2]
                        for line in combined.read_text().splitlines()
                        if len(line.split("\t")) == 7
                    )
                    if flank_events != accepted_events:
                        errors.append(
                            "accepted OUTSIDER candidates do not match combined flanks"
                        )

    if errors:
        print("INSIDER/OUTSIDER migration check: FAILED")
        for error in errors:
            print(f"- {error}")
        return 1
    print("INSIDER/OUTSIDER migration check: PASSED")
    print(f"Exact INSIDER files: {len(INSIDER_FILES)}")
    print(f"Exact INSIDER TE files: {len(INSIDER_TE_FILES) + 1}")
    print("OUTSIDER VCF: identical after volatile header normalization")
    print(f"Exact OUTSIDER TE files: {len(OUTSIDER_TE_EXACT_FILES)}")
    print(f"Recovered OUTSIDER CIGAR candidates: {recovered_cigar_candidates}")
    print(f"Exact OUTSIDER flank/TSD files: {len(OUTSIDER_TSD_EXACT_FILES)}")
    print(
        "Exact OUTSIDER flank artifact directories: "
        f"{len(OUTSIDER_TSD_EXACT_DIRECTORIES)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
