#!/usr/bin/env python3
"""Convert lifted flank pairs into OUTSIDER insertion coordinates."""

from __future__ import print_function

import argparse
import csv
from collections import OrderedDict
from pathlib import Path


AUDIT_HEADER = (
    "event_id",
    "te_family",
    "status",
    "reason",
    "left_chromosome",
    "left_end",
    "left_coverage",
    "right_chromosome",
    "right_start",
    "right_coverage",
    "projected_chromosome",
    "projected_start",
    "projected_end",
    "projected_gap",
)


def attributes(text):
    values = {}
    for field in text.strip().strip(";").split(";"):
        if "=" in field:
            key, value = field.split("=", 1)
            values[key] = value
    return values


def read_lifted_features(path):
    grouped = OrderedDict()
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    "{}:{}: expected nine GFF columns".format(path, line_number)
                )
            values = attributes(fields[8])
            identifier = values.get("ID", "")
            side = values.get("SIDE", "")
            if not identifier or side not in ("L", "R"):
                continue
            grouped.setdefault(identifier, []).append(
                {
                    "chromosome": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "family": values.get("NAME", "UNKNOWN"),
                    "side": side,
                    "coverage": float(values.get("coverage", "0") or 0),
                }
            )
    return grouped


def summarize(grouped, max_gap):
    good = []
    bad = []
    mapped_ids = []
    audit = []
    for identifier, features in grouped.items():
        lefts = [feature for feature in features if feature["side"] == "L"]
        rights = [feature for feature in features if feature["side"] == "R"]
        row = {
            "event_id": identifier,
            "te_family": (
                features[0]["family"] if features else "UNKNOWN"
            ),
            "status": "rejected",
            "reason": "",
            "left_chromosome": lefts[0]["chromosome"] if len(lefts) == 1 else "",
            "left_end": str(lefts[0]["end"]) if len(lefts) == 1 else "",
            "left_coverage": str(lefts[0]["coverage"]) if len(lefts) == 1 else "",
            "right_chromosome": rights[0]["chromosome"] if len(rights) == 1 else "",
            "right_start": str(rights[0]["start"]) if len(rights) == 1 else "",
            "right_coverage": str(rights[0]["coverage"]) if len(rights) == 1 else "",
            "projected_chromosome": "",
            "projected_start": "",
            "projected_end": "",
            "projected_gap": "",
        }
        if len(lefts) != 1 or len(rights) != 1:
            row["reason"] = "missing_or_ambiguous_flank"
            audit.append(row)
            continue
        mapped_ids.append(identifier)
        left = lefts[0]
        right = rights[0]
        family = right["family"]
        name = "{}|{}".format(family, identifier)
        gap = right["start"] - left["end"]
        row["projected_gap"] = str(gap)
        if left["chromosome"] == right["chromosome"]:
            if abs(gap) <= max_gap:
                start = min(left["end"], right["start"])
                end = max(left["end"], right["start"])
                good.append((left["chromosome"], start, end, name, gap, "INSIDER"))
                row["status"] = "projected"
                row["reason"] = "concordant_flanks"
                row["projected_chromosome"] = left["chromosome"]
                row["projected_start"] = str(start)
                row["projected_end"] = str(end)
            else:
                row["reason"] = "gap_exceeds_limit"
                bad.append(
                    (
                        "{}:{}".format(right["chromosome"], right["start"]),
                        "{}:{}".format(left["chromosome"], left["end"]),
                        left["end"],
                        name,
                        gap,
                        "INSIDER",
                    )
                )
        else:
            # A pair projected to two chromosomes is not a defensible locus.
            # The historical shell selected one flank by coverage and emitted
            # it as a regular call, sometimes combining the selected left
            # coordinate with the right chromosome.  Keep the diagnostic in
            # BAD_POS_TE_LIFT.bed, but do not promote it to the public BED.
            selected = left if left["coverage"] >= right["coverage"] else right
            bad.append(
                (
                    "{}:{}".format(right["chromosome"], right["start"]),
                    "{}:{}".format(left["chromosome"], left["end"]),
                    left["end"],
                    name,
                    gap,
                    "INSIDER",
                    selected["side"],
                )
            )
            row["reason"] = "discordant_chromosomes"
        audit.append(row)
    return good, bad, mapped_ids, audit


def write_rows(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as handle:
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def write_audit(path, rows):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build(args):
    grouped = read_lifted_features(args.lifted_gff)
    good, bad, mapped_ids, audit = summarize(grouped, args.max_gap)
    write_rows(args.good_bed, good)
    write_rows(args.bad_bed, bad)
    write_rows(args.mapped_ids, ((identifier,) for identifier in mapped_ids))
    write_rows(args.public_bed, (row[:4] for row in good))
    write_audit(args.audit, audit)

    combined = []
    if args.insider_bed and Path(args.insider_bed).is_file():
        with Path(args.insider_bed).open() as handle:
            combined.extend(line.rstrip("\r\n") for line in handle if line.strip())
    combined.extend("\t".join(str(value) for value in row[:4]) for row in good)
    Path(args.combined_bed).parent.mkdir(parents=True, exist_ok=True)
    with Path(args.combined_bed).open("w") as handle:
        for line in combined:
            handle.write(line + "\n")
    return len(good), len(bad)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lifted-gff", type=Path, required=True)
    parser.add_argument("--insider-bed", type=Path)
    parser.add_argument("--good-bed", type=Path, required=True)
    parser.add_argument("--bad-bed", type=Path, required=True)
    parser.add_argument("--mapped-ids", type=Path, required=True)
    parser.add_argument("--public-bed", type=Path, required=True)
    parser.add_argument("--combined-bed", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--max-gap", type=int, default=20000)
    args = parser.parse_args(argv)
    if args.max_gap < 0:
        parser.error("--max-gap must be zero or greater")
    return args


def main(argv=None):
    args = parse_args(argv)
    good, bad = build(args)
    print("Projected {} insertions; {} mappings require review.".format(good, bad))


if __name__ == "__main__":
    main()
