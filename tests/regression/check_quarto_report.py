#!/usr/bin/env python3
"""Check that a rendered Quarto report is a faithful TE_INFOS projection."""

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source(tremolo_id):
    if "_INSIDER" in tremolo_id:
        return "INSIDER"
    if "_OUTSIDER" in tremolo_id:
        return "OUTSIDER"
    return "UNKNOWN"


def optional_number(value):
    if value in {"", ".", "NONE", "NONE-FREQA", "NONE-FREQB", "INSIDER"}:
        return None
    return float(value.replace(",", "."))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work_directory", type=Path)
    args = parser.parse_args()
    te_infos = args.work_directory / "TE_INFOS.bed"
    report_data = args.work_directory / "REPORT/report-data.json"
    report_html = args.work_directory / "REPORT/report.html"
    for path in (te_infos, report_data, report_html):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing or empty report artifact: {path}")

    with te_infos.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if len(rows[0]) != 15 or len(rows) < 1:
        raise SystemExit("invalid TE_INFOS schema")
    calls = rows[1:]
    data = json.loads(report_data.read_text())
    summary = data["summary"]
    expected_sources = Counter(source(row[13]) for row in calls)
    expected_events = Counter(row[14] for row in calls)
    expected_strands = Counter(row[4] for row in calls)
    expected_tsd = sum(row[5] not in {"", ".", "NONE"} for row in calls)
    expected_frequencies = []
    for row in calls:
        raw = optional_number(row[10])
        clipped = optional_number(row[11])
        expected_frequencies.append(clipped if clipped is not None else raw)

    expected = {
        "calls": len(calls),
        "families": len({row[3].partition("|")[0] for row in calls}),
        "chromosomes_with_calls": len({row[0] for row in calls}),
        "tsd_confirmed": expected_tsd,
        "frequency_available": sum(value is not None for value in expected_frequencies),
        "sources": dict(sorted(expected_sources.items())),
        "event_types": dict(sorted(expected_events.items())),
        "strands": dict(sorted(expected_strands.items())),
    }
    for key, value in expected.items():
        if summary.get(key) != value:
            raise SystemExit(f"report summary mismatch for {key}: {summary.get(key)!r} != {value!r}")
    if data["report"]["te_infos_sha256"] != sha256(te_infos):
        raise SystemExit("report TE_INFOS checksum mismatch")
    if len(data["calls"]) != len(calls):
        raise SystemExit("report call projection is incomplete")
    for line_number, (row, projected, frequency) in enumerate(
        zip(calls, data["calls"], expected_frequencies), 2
    ):
        expected_projection = {
            "chrom": row[0],
            "start": int(row[1]),
            "event_type": row[14],
            "display_frequency": frequency,
        }
        observed_projection = {
            key: projected.get(key) for key in expected_projection
        }
        if observed_projection != expected_projection:
            raise SystemExit(
                f"report frequency-position projection mismatch at TE_INFOS line {line_number}: "
                f"{observed_projection!r} != {expected_projection!r}"
            )

    html = report_html.read_text(errors="replace")
    required_ids = (
        "tremolo-report-data",
        "trm-summary-cards",
        "trm-family-chart",
        "trm-frequency-position-chart",
        "trm-frequency-position-chrom",
        "trm-frequency-position-legend",
        "trm-call-table-body",
        "trm-proximity-table-body",
    )
    for identifier in required_ids:
        if f'id="{identifier}"' not in html:
            raise SystemExit(f"rendered HTML lacks {identifier}")

    candidates_path = args.work_directory / "TE_CALL_CANDIDATES.tsv"
    ambiguous = data.get("ambiguous_calls", {})
    if not candidates_path.is_file() or candidates_path.stat().st_size == 0:
        raise SystemExit("missing normalized TE_CALL_CANDIDATES.tsv")
    with candidates_path.open(newline="") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    expected_groups = {row["candidate_group_id"] for row in candidate_rows}
    if ambiguous.get("summary", {}).get("ambiguous_calls") != len(expected_groups):
        raise SystemExit("ambiguous-call group count mismatch")
    if len(ambiguous.get("candidates", [])) != len(candidate_rows):
        raise SystemExit("ambiguous-call candidate projection is incomplete")
    if any(int(row["candidate_count"]) < 2 for row in candidate_rows):
        raise SystemExit("TE_CALL_CANDIDATES contains a non-ambiguous call")
    if 'id="trm-ambiguous-call-table-body"' not in html:
        raise SystemExit("rendered HTML lacks ambiguous-call table")

    resident = data.get("resident_te", {})
    if resident.get("available"):
        resident_rows = []
        targets = resident.get("summary", {}).get("targets", {})
        for target in targets:
            copies_path = args.work_directory / f"TE_GENOME/{target}/ALL_TE_COPIES.tsv"
            if not copies_path.is_file():
                raise SystemExit(
                    f"report exposes resident TE data but its {target} source table is missing"
                )
            with copies_path.open(newline="") as handle:
                resident_rows.extend(csv.DictReader(handle, delimiter="\t"))
        projected = resident.get("copies", [])
        if len(projected) != len(resident_rows):
            raise SystemExit("resident TE report projection is incomplete")
        expected_ambiguous = sum(int(row["candidate_count"]) > 1 for row in resident_rows)
        expected_multi_te = sum(len(row["te_candidates"].split(";")) > 1 for row in resident_rows)
        if resident["summary"].get("ambiguous_copies") != expected_ambiguous:
            raise SystemExit("resident TE ambiguous-copy count mismatch")
        if resident["summary"].get("multi_te_copies") != expected_multi_te:
            raise SystemExit("resident TE multi-label copy count mismatch")
        if any(
            copy["candidate_count"] > 1
            and len(copy.get("candidates", [])) != copy["candidate_count"]
            for copy in projected
        ):
            raise SystemExit("resident TE alternative assignments are incomplete")
        for identifier in (
            "trm-resident-cards",
            "trm-resident-calibration-chart",
            "trm-resident-table-body",
        ):
            if f'id="{identifier}"' not in html:
                raise SystemExit(f"rendered HTML lacks {identifier}")
    if re.search(r"<(?:script|link)[^>]+(?:src|href)=[\"']https?://", html, re.I):
        raise SystemExit("rendered report loads an external script or stylesheet")
    print(
        "Quarto report regression: PASSED "
        f"({len(calls)} calls, {summary['families']} families, "
        f"{len(data['proximity_groups'])} proximity groups)"
    )


if __name__ == "__main__":
    main()
