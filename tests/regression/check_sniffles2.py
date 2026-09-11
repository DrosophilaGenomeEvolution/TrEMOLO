#!/usr/bin/env python3
"""Check a completed bundled Sniffles 2 run without asserting legacy counts."""

import argparse
import csv
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("--expect-empty", action="store_true")
    args = parser.parse_args()
    variant = args.workdir / "OUTSIDER/VARIANT_CALLING"
    metadata = json.loads((variant / "sniffles-run.json").read_text())
    assert metadata["caller"] == "sniffles2"
    assert metadata["version"] == "2.8.1"
    assert "--reference" in metadata["command"] and "--output-rnames" in metadata["command"]
    expected = {}
    count = 0
    for line in (variant / "SV.vcf").read_text().splitlines():
        if line.startswith("#"):
            continue
        count += 1
        fields = line.split("\t")
        info = dict(item.split("=", 1) if "=" in item else (item, True)
                    for item in fields[7].split(";"))
        # The bundled configuration selects insertions >1000 bp. Sniffles 2
        # ran with a reference, so the ALT includes exactly one anchor base.
        if info["SVTYPE"] != "INS" or fields[4].startswith("<"):
            continue
        length = abs(int(info["SVLEN"]))
        assert len(fields[4]) == length + 1
        assert fields[4][0] == fields[3]
        sequence = fields[4][1:]
        if length <= 1000 or set(sequence) <= {"N"}:
            continue
        identifier = "sniffles.INS." + fields[2]
        expected[identifier] = (fields[0], int(fields[1]), int(info["SUPPORT"]), sequence)
    observed = {}
    lines = (variant / "SEQUENCE_INDEL.fasta").read_text().splitlines()
    for index, line in enumerate(lines):
        if not line.startswith(">"):
            continue
        parts = line[1:].split(":")
        if parts[-1] == "0":
            observed[parts[4]] = (parts[0], int(parts[2]), int(parts[5]), lines[index + 1])
    assert observed == expected, "Normalized consensus sequences, IDs, coordinates or support differ"
    if args.expect_empty:
        assert count == 0 and not observed
    else:
        assert count > 0 and observed, "Fixture did not exercise Sniffles 2 sequence extraction"
    calls = [line.split("\t") for line in (args.workdir / "TE_INFOS.bed").read_text().splitlines()[1:]]
    assert calls and all(len(row) == 15 for row in calls)
    with (args.workdir / "OUTSIDER/TE_TOWARD_GENOME/INTEGRATION_TE.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            assert row["status"] == "integrated", row
    print("Sniffles 2 regression: PASSED ({} SVs, {} consensus sequences, {} final calls)".format(
        count, len(observed), len(calls)))


if __name__ == "__main__":
    main()
