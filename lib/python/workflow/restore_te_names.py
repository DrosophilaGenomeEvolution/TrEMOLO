#!/usr/bin/env python3
"""Publish original TE names in one streaming pass per declared output."""
import argparse
import csv
import json
import os
import shutil
import tempfile
from pathlib import Path


def read_names(path):
    names = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["original", "prepared"]:
            raise ValueError("invalid TE-header mapping: {}".format(path))
        for row in reader:
            original, prepared = row["original"], row["prepared"]
            if not original or not prepared or prepared in names:
                raise ValueError("empty or duplicate prepared TE identifier")
            if any(c in original for c in "\t\r\n"):
                raise ValueError("TE names cannot contain tab or newline characters")
            names[prepared] = original
    return names


def publish(source, destination, kind, names):
    source, destination = Path(source), Path(destination)
    if source.resolve() == destination.resolve():
        raise ValueError("public output must differ from internal input")
    destination.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix="." + destination.name, dir=str(destination.parent))
    os.close(descriptor)
    changed = 0
    try:
        if kind == "copy" or all(key == value for key, value in names.items()):
            shutil.copyfile(str(source), temporary)
        else:
            with source.open(newline="") as reader_handle, open(temporary, "w", newline="") as writer_handle:
                reader = csv.reader(reader_handle, delimiter="\t")
                writer = csv.writer(writer_handle, delimiter="\t", lineterminator="\n")
                if kind == "table":
                    header = next(reader, None)
                    indices = [] if header is None else [i for i, value in enumerate(header)
                        if value in {"sseqid", "te_family", "reported_te", "candidate_te", "te_name"}]
                    if header is not None:
                        writer.writerow(header)
                    compound = False
                else:
                    mode, index = kind.split(":")
                    if mode not in ("compound", "column"):
                        raise ValueError("unknown name restoration kind: " + kind)
                    indices, compound = [int(index)], mode == "compound"
                for row in reader:
                    if not row:
                        writer.writerow(row)
                        continue
                    for index in indices:
                        if index >= len(row):
                            raise ValueError("{}: missing TE column {}".format(source, index))
                        value = row[index]
                        if compound:
                            family, separator, identifier = value.partition("|")
                            restored = names.get(family, family) + separator + identifier
                        else:
                            restored = names.get(value, value)
                        changed += restored != value
                        row[index] = restored
                    writer.writerow(row)
        os.replace(temporary, str(destination))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return changed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--te-headers", type=Path, required=True)
    parser.add_argument("--exports", required=True, help="JSON list of source/output/kind records")
    parser.add_argument("--manifest", type=Path, required=True)
    args = parser.parse_args()
    names = read_names(args.te_headers)
    exports = json.loads(args.exports)
    outputs = [Path(item["output"]).resolve() for item in exports]
    if len(outputs) != len(set(outputs)):
        raise ValueError("duplicate public output")
    results = []
    for item in exports:
        count = publish(item["source"], item["output"], item["kind"], names)
        results.append(dict(item, restored_cells=count))
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps({"te_records": len(names), "exports": results}, indent=2) + "\n")
    print("Published {} files; {} TE-name cells restored.".format(len(results), sum(r["restored_cells"] for r in results)))


if __name__ == "__main__":
    main()
