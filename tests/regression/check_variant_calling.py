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

OUTSIDER_FREQUENCY_EXACT_FILES = (
    "1-UTILS/TE_SIZE.tsv",
    "OUTSIDER/FREQUENCY/FILTER_BLAST_INS.csv",
    "OUTSIDER/FREQUENCY/FILTER_BLAST_INS.sorted.csv",
    "OUTSIDER/FREQUENCY/SV_SIZE.tsv",
    "OUTSIDER/FREQUENCY/COUNT_READS.txt",
    "OUTSIDER/FREQUENCY/COUNT_TE_IN_RS.txt",
    "OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv",
    "OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv",
)

OUTSIDER_FREQUENCY_EXPECTED_COUNTS = {
    "1-UTILS/TE_SIZE.tsv": 179,
    "OUTSIDER/FREQUENCY/FILTER_BLAST_INS.csv": 70,
    "OUTSIDER/FREQUENCY/FILTER_BLAST_INS.sorted.csv": 70,
    "OUTSIDER/FREQUENCY/SV_SIZE.tsv": 1213,
    "OUTSIDER/FREQUENCY/COUNT_READS.txt": 4295,
    "OUTSIDER/FREQUENCY/COUNT_TE_IN_RS.txt": 60,
    "OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv": 34,
    "OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv": 29,
}

OUTSIDER_FREQUENCY_EXPECTED_CANDIDATE_SOURCES = Counter(
    {"TrEMOLO": 19, "sniffles": 15, "HARD": 35}
)

OUTSIDER_FREQUENCY_EXPECTED_TYPES = Counter({"INS": 29, "DEL": 5})

OUTSIDER_FREQUENCY_REQUIRED_ARTIFACTS = (
    "OUTSIDER/FREQUENCY/MAPPING_POSTION_TE.bam",
    "OUTSIDER/FREQUENCY/MAPPING_POSTION_TE.bam.bai",
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

    for relative in OUTSIDER_FREQUENCY_EXACT_FILES:
        old = args.legacy / relative
        new = args.migrated / relative
        if not old.is_file() or not new.is_file():
            errors.append(f"missing OUTSIDER frequency output: {relative}")
        elif digest(old) != digest(new):
            errors.append(f"OUTSIDER frequency output differs: {relative}")

    for relative, expected_count in OUTSIDER_FREQUENCY_EXPECTED_COUNTS.items():
        path = args.migrated / relative
        if path.is_file():
            actual_count = sum(1 for _ in path.open("rb"))
            if actual_count != expected_count:
                errors.append(
                    f"OUTSIDER frequency line count differs: {relative} "
                    f"({actual_count} != {expected_count})"
                )

    frequency_regions = (
        args.migrated / "OUTSIDER/FREQUENCY/POSITION_START_TE.bed"
    )
    if not frequency_regions.is_file():
        errors.append("missing OUTSIDER frequency region table")
    elif sum(1 for _ in frequency_regions.open("rb")) != 51:
        errors.append("OUTSIDER frequency region count differs from 51")

    for relative in OUTSIDER_FREQUENCY_REQUIRED_ARTIFACTS:
        artifact = args.migrated / relative
        if not artifact.is_file() or artifact.stat().st_size == 0:
            errors.append(f"missing or empty OUTSIDER frequency artifact: {relative}")

    frequency_directory = args.migrated / "OUTSIDER/FREQUENCY"
    candidates = frequency_directory / "FILTER_BLAST_INS.csv"
    frequency = frequency_directory / "FREQUENCY_TE_INS.tsv"
    precise_frequency = frequency_directory / "FREQUENCY_TE_INS_PRECISE.tsv"
    merged_calls = (
        args.migrated
        / "OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed"
    )

    candidate_ids = []
    if candidates.is_file():
        candidate_rows = candidates.read_text().splitlines()[1:]
        candidate_sources = Counter()
        for row in candidate_rows:
            fields = row.split("\t")
            query_parts = fields[1].split(":") if len(fields) > 1 else []
            if len(query_parts) <= 4:
                errors.append("malformed OUTSIDER frequency candidate row")
                continue
            event_id = query_parts[4]
            candidate_ids.append(event_id)
            candidate_sources[event_id.split(".", 1)[0]] += 1
        if candidate_sources != OUTSIDER_FREQUENCY_EXPECTED_CANDIDATE_SOURCES:
            errors.append(
                "OUTSIDER frequency candidate sources differ: "
                f"{dict(candidate_sources)} != "
                f"{dict(OUTSIDER_FREQUENCY_EXPECTED_CANDIDATE_SOURCES)}"
            )
        if len(candidate_ids) != 69 or len(set(candidate_ids)) != 69:
            errors.append(
                "OUTSIDER frequency candidates are not 69 unique events"
            )

    frequency_ids = []
    frequency_types = Counter()
    frequency_sources = Counter()
    if frequency.is_file():
        for row in frequency.read_text().splitlines():
            fields = row.split("\t")
            if len(fields) < 11:
                errors.append("malformed OUTSIDER frequency row")
                continue
            event_id = fields[1]
            frequency_ids.append(event_id)
            frequency_sources[event_id.split(".", 1)[0]] += 1
            frequency_types[fields[10]] += 1
        if frequency_types != OUTSIDER_FREQUENCY_EXPECTED_TYPES:
            errors.append(
                "OUTSIDER frequency event types differ: "
                f"{dict(frequency_types)} != "
                f"{dict(OUTSIDER_FREQUENCY_EXPECTED_TYPES)}"
            )
        expected_sources = Counter({"TrEMOLO": 19, "sniffles": 15})
        if frequency_sources != expected_sources:
            errors.append(
                "OUTSIDER frequency event sources differ: "
                f"{dict(frequency_sources)} != {dict(expected_sources)}"
            )
        if len(frequency_ids) != 34 or len(set(frequency_ids)) != 34:
            errors.append("OUTSIDER frequencies are not 34 unique events")
        retained_candidate_ids = {
            event_id
            for event_id in candidate_ids
            if not event_id.startswith(("HARD.", "SOFT."))
        }
        if candidate_ids and set(frequency_ids) != retained_candidate_ids:
            errors.append(
                "OUTSIDER frequency output does not retain exactly the "
                "non-HARD/SOFT compatibility candidates"
            )

    precise_ids = []
    precise_types = Counter()
    if precise_frequency.is_file():
        for row in precise_frequency.read_text().splitlines():
            fields = row.split("\t")
            if len(fields) < 11:
                errors.append("malformed OUTSIDER precise-frequency row")
                continue
            precise_ids.append(fields[1])
            precise_types[fields[10]] += 1
        if precise_types != Counter({"INS": 29}):
            errors.append(
                "OUTSIDER precise-frequency event types differ: "
                f"{dict(precise_types)} != {{'INS': 29}}"
            )
        frequency_ins_ids = {
            row.split("\t")[1]
            for row in frequency.read_text().splitlines()
            if len(row.split("\t")) >= 11 and row.split("\t")[10] == "INS"
        } if frequency.is_file() else set()
        if frequency_ins_ids and set(precise_ids) != frequency_ins_ids:
            errors.append(
                "OUTSIDER precise frequencies do not cover all frequency INS events"
            )

    if frequency_ids and merged_calls.is_file():
        merged_ids = set()
        for row in merged_calls.read_text().splitlines():
            fields = row.split("\t")
            if len(fields) < 4 or "|" not in fields[3]:
                errors.append("malformed merged OUTSIDER TE row")
                continue
            merged_ids.add(fields[3].split("|", 1)[1])
        in_merged_calls = len(set(frequency_ids) & merged_ids)
        outside_merged_calls = len(set(frequency_ids) - merged_ids)
        if (in_merged_calls, outside_merged_calls) != (25, 9):
            errors.append(
                "OUTSIDER frequency merged-call membership differs: "
                f"{in_merged_calls} merged and {outside_merged_calls} extra "
                "events != 25 and 9"
            )

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
    print(
        "Exact OUTSIDER frequency files: "
        f"{len(OUTSIDER_FREQUENCY_EXACT_FILES)}"
    )
    print(
        "OUTSIDER frequency semantics: 69 candidates, 34 frequencies "
        "(25 merged + 9 extra), 29 precise INS"
    )
    print("OUTSIDER frequency regional BAM/index: present and non-empty")
    print(f"Recovered OUTSIDER CIGAR candidates: {recovered_cigar_candidates}")
    print(f"Exact OUTSIDER flank/TSD files: {len(OUTSIDER_TSD_EXACT_FILES)}")
    print(
        "Exact OUTSIDER flank artifact directories: "
        f"{len(OUTSIDER_TSD_EXACT_DIRECTORIES)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
