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

    expected = {
        "calls": len(calls),
        "families": len({row[3].partition("|")[0] for row in calls}),
        "chromosomes_with_calls": len({row[0] for row in calls}),
        "tsd_confirmed": expected_tsd,
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

    html = report_html.read_text(errors="replace")
    required_ids = (
        "tremolo-report-data",
        "trm-summary-cards",
        "trm-family-chart",
        "trm-call-table-body",
        "trm-proximity-table-body",
    )
    for identifier in required_ids:
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

