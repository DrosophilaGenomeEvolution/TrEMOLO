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
    "left_strand",
    "right_strand",
    "projected_strand",
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
                    "strand": fields[6],
                }
            )
    return grouped


def read_expected_insertions(path):
    """Keep every integrated event, including those with no source flanks."""

    expected = OrderedDict()
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 4 or ":" not in fields[3]:
                raise ValueError("{}:{}: expected integrated TE:event BED record".format(path, line_number))
            family, identifier = fields[3].rsplit(":", 1)
            if not family or not identifier:
                raise ValueError("{}:{}: empty TE family or event ID".format(path, line_number))
            if identifier in expected:
                raise ValueError("duplicate integrated event {}".format(identifier))
            expected[identifier] = family
    return expected


def summarize(grouped, max_gap, expected=None):
    if expected is None:
        expected = OrderedDict((identifier, features[0]["family"] if features else "UNKNOWN")
                               for identifier, features in grouped.items())
    unexpected = set(grouped).difference(expected)
    if unexpected:
        raise ValueError("Lifted GFF contains unexpected events: {}".format(
            ", ".join(sorted(unexpected))))
    good = []
    bad = []
    mapped_ids = []
    audit = []
    for identifier, family in expected.items():
        features = grouped.get(identifier, [])
        lefts = [feature for feature in features if feature["side"] == "L"]
        rights = [feature for feature in features if feature["side"] == "R"]
        row = {
            "event_id": identifier,
            "te_family": family,
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
            "left_strand": lefts[0]["strand"] if len(lefts) == 1 else "",
            "right_strand": rights[0]["strand"] if len(rights) == 1 else "",
            "projected_strand": "",
        }
        if not features:
            row["reason"] = "no_mapped_flanks"
            audit.append(row)
            continue
        if len(lefts) != 1 or len(rights) != 1:
            row["reason"] = "missing_or_ambiguous_flank"
            audit.append(row)
            continue
        mapped_ids.append(identifier)
        left = lefts[0]
        right = rights[0]
        name = "{}|{}".format(right["family"], identifier)
        if left["chromosome"] != right["chromosome"]:
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
                    right["start"] - 1 - left["end"],
                    "INSIDER",
                    selected["side"],
                )
            )
            row["reason"] = "discordant_chromosomes"
        elif left["family"] != right["family"]:
            row["reason"] = "discordant_families"
        elif left["strand"] not in ("+", "-") or right["strand"] not in ("+", "-"):
            row["reason"] = "unknown_flank_strand"
        elif left["strand"] != right["strand"]:
            row["reason"] = "discordant_flank_strands"
        elif min(left["start"], right["start"]) < 1 or left["end"] < left["start"] or right["end"] < right["start"]:
            row["reason"] = "invalid_flank_interval"
        else:
            # Inner boundaries in zero-based, half-open coordinates. A pair
            # mapped on '-' reverses both the order and the relevant ends.
            if left["strand"] == "+":
                ordered = left["start"] <= right["start"] and left["end"] <= right["end"]
                lower, upper = left["end"], right["start"] - 1
            else:
                ordered = right["start"] <= left["start"] and right["end"] <= left["end"]
                lower, upper = right["end"], left["start"] - 1
            gap = upper - lower
            row["projected_gap"] = str(gap)
            if not ordered:
                row["reason"] = "discordant_flank_order"
            elif abs(gap) > max_gap:
                row["reason"] = "gap_exceeds_limit"
                bad.append(("{}:{}".format(left["chromosome"], upper),
                            "{}:{}".format(left["chromosome"], lower),
                            lower, name, gap, "INSIDER"))
            else:
                # Small overlaps can reflect a target-site duplication;
                # retain the signed gap in the audit rather than rejecting it.
                start, end = min(lower, upper), max(lower, upper)
                good.append((left["chromosome"], start, end, name, gap, "INSIDER"))
                row["status"] = "projected"
                row["reason"] = "concordant_flanks"
                row["projected_chromosome"] = left["chromosome"]
                row["projected_start"] = str(start)
                row["projected_end"] = str(end)
                row["projected_strand"] = left["strand"]
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
    expected = read_expected_insertions(args.positions)
    good, bad, mapped_ids, audit = summarize(grouped, args.max_gap, expected)
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
    return len(good), sum(row["status"] == "rejected" for row in audit)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lifted-gff", type=Path, required=True)
    parser.add_argument("--positions", type=Path, required=True,
                        help="Canonical integrated BED defining every expected insertion")
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
    good, rejected = build(args)
    print("Projected {} insertions; {} mappings rejected.".format(good, rejected))


if __name__ == "__main__":
    main()
