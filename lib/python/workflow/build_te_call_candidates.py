#!/usr/bin/env python3
"""Build a normalized table containing only ambiguous variable TE calls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional


HEADER = (
    "candidate_group_id",
    "source",
    "chrom",
    "start",
    "end",
    "event_id",
    "tremolo_id",
    "event_type",
    "reported_te",
    "candidate_te",
    "assignment",
    "candidate_rank",
    "candidate_count",
    "evidence_count",
    "evidence_fraction",
    "evidence_channels",
    "ambiguity_type",
)


def classify_source(tremolo_id: str) -> str:
    if "_INSIDER" in tremolo_id:
        return "INSIDER"
    if "_OUTSIDER" in tremolo_id:
        return "OUTSIDER"
    return "UNKNOWN"


def stable_group_id(source: str, chrom: str, start: int, event_id: str) -> str:
    value = "\t".join((source, chrom, str(start), event_id)).encode("utf-8")
    return "TCA" + hashlib.sha1(value).hexdigest()[:13]


def read_calls(path: Path) -> dict[str, dict]:
    calls = {}
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"{path}: empty TE_INFOS table") from error
        if len(header) != 15 or header[:4] != ["#chrom", "start", "end", "TE|ID"]:
            raise ValueError(f"{path}: unexpected TE_INFOS header")
        for line_number, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) != 15:
                raise ValueError(f"{path}:{line_number}: expected 15 columns")
            family, separator, event_id = row[3].partition("|")
            if not separator or not family or not event_id:
                raise ValueError(f"{path}:{line_number}: malformed TE|ID")
            try:
                start = int(row[1])
                end = int(row[2])
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: non-integer coordinates") from error
            if start < 0 or end < start:
                raise ValueError(f"{path}:{line_number}: invalid interval")
            call = {
                "source": classify_source(row[13]),
                "chrom": row[0],
                "start": start,
                "end": end,
                "event_id": event_id,
                "tremolo_id": row[13],
                "event_type": row[14],
                "reported_te": family,
            }
            previous = calls.get(event_id)
            if previous is not None and previous != call:
                raise ValueError(f"{path}:{line_number}: conflicting duplicate event {event_id}")
            calls[event_id] = call
    return calls


def parse_evidence_argument(value: str) -> tuple[str, Path]:
    channel, separator, raw_path = value.partition("=")
    if not separator or not channel or not raw_path:
        raise argparse.ArgumentTypeError("evidence must use CHANNEL=PATH")
    return channel, Path(raw_path)


def read_evidence(
    inputs: list[tuple[str, Path]],
) -> dict[str, dict[str, dict[str, object]]]:
    evidence: dict[str, dict[str, dict[str, object]]] = defaultdict(
        lambda: defaultdict(lambda: {"count": 0, "channels": set()})
    )
    for channel, path in inputs:
        if str(path) == "/dev/null":
            continue
        with path.open() as handle:
            for line_number, raw_line in enumerate(handle, 1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                fields = line.split("\t")
                if len(fields) != 2:
                    raise ValueError(f"{path}:{line_number}: expected two columns")
                event_and_te, count_text = fields
                event_id, separator, candidate_te = event_and_te.rpartition(":")
                if not separator or not event_id or not candidate_te:
                    raise ValueError(f"{path}:{line_number}: malformed event:TE key")
                try:
                    count = int(count_text)
                except ValueError as error:
                    raise ValueError(f"{path}:{line_number}: non-integer evidence count") from error
                if count < 1:
                    raise ValueError(f"{path}:{line_number}: evidence count must be positive")
                item = evidence[event_id][candidate_te]
                item["count"] = int(item["count"]) + count
                channels = item["channels"]
                assert isinstance(channels, set)
                channels.add(channel)
    return evidence


def format_fraction(value: float) -> str:
    return format(value, ".6g")


def build_rows(calls: dict[str, dict], evidence: dict) -> list[tuple[str, ...]]:
    rows = []
    for event_id, candidates in evidence.items():
        call = calls.get(event_id)
        if call is None or len(candidates) < 2:
            continue
        if call["reported_te"] not in candidates:
            raise ValueError(
                f"{event_id}: reported TE {call['reported_te']} is absent from candidates"
            )
        ordered = sorted(
            candidates.items(),
            key=lambda item: (
                -int(item[1]["count"]),
                item[0] != call["reported_te"],
                item[0],
            ),
        )
        total = sum(int(item["count"]) for _, item in ordered)
        group_id = stable_group_id(
            call["source"], call["chrom"], call["start"], event_id
        )
        for rank, (candidate_te, item) in enumerate(ordered, 1):
            channels = item["channels"]
            assert isinstance(channels, set)
            count = int(item["count"])
            rows.append(
                (
                    group_id,
                    call["source"],
                    call["chrom"],
                    str(call["start"]),
                    str(call["end"]),
                    event_id,
                    call["tremolo_id"],
                    call["event_type"],
                    call["reported_te"],
                    candidate_te,
                    "reported_primary" if candidate_te == call["reported_te"] else "alternative",
                    str(rank),
                    str(len(ordered)),
                    str(count),
                    format_fraction(count / total),
                    ";".join(sorted(channels)),
                    "unresolved_family_candidates",
                )
            )
    return sorted(rows, key=lambda row: (row[1], row[2], int(row[3]), row[5], int(row[11])))


def atomic_write(path: Path, rows: list[tuple[str, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=str(path.parent), text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
        os.replace(str(temporary), str(path))
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def build(te_infos: Path, evidence_inputs: list[tuple[str, Path]], output: Path) -> tuple[int, int]:
    calls = read_calls(te_infos)
    evidence = read_evidence(evidence_inputs)
    rows = build_rows(calls, evidence)
    atomic_write(output, rows)
    return len({row[0] for row in rows}), len(rows)


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--te-infos", required=True, type=Path)
    parser.add_argument(
        "--evidence",
        action="append",
        default=[],
        type=parse_evidence_argument,
        metavar="CHANNEL=PATH",
    )
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = parse_args(argv)
    group_count, row_count = build(args.te_infos, args.evidence, args.output)
    print(
        f"Wrote {row_count} candidates for {group_count} ambiguous variable TE calls."
    )


if __name__ == "__main__":
    main()
