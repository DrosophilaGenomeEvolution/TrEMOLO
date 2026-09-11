#!/usr/bin/env python3
"""Describe both flanks of integrated OUTSIDER insertions as GFF3 features."""

from __future__ import print_function

import argparse
from collections import OrderedDict
from pathlib import Path


def read_fai(path):
    lengths = OrderedDict()
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 2:
                raise ValueError("{}:{}: malformed FASTA index".format(path, line_number))
            lengths[fields[0]] = int(fields[1])
    return lengths


def prepare(position_bed, genome_index, output_gff, feature_file, flank_size):
    lengths = read_fai(genome_index)
    features = []
    seen = set()
    with Path(position_bed).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    "{}:{}: expected four BED columns".format(position_bed, line_number)
                )
            chromosome, start_text, end_text, name = fields[:4]
            if ".INS." not in name:
                continue
            if ":" not in name:
                raise ValueError(
                    "{}:{}: expected TE:event name".format(position_bed, line_number)
                )
            family, identifier = name.rsplit(":", 1)
            if identifier in seen:
                raise ValueError("duplicate integrated event {}".format(identifier))
            seen.add(identifier)
            if chromosome not in lengths:
                raise ValueError("{} is absent from {}".format(chromosome, genome_index))
            start = int(start_text)
            end = int(end_text)
            if start < 0 or end <= start or end > lengths[chromosome]:
                raise ValueError(
                    "invalid integrated interval {}:{}-{}".format(
                        chromosome, start, end
                    )
                )
            # Keep the historical coordinates here: the left feature ends at
            # BED start and the right feature starts at BED end.  GFF itself is
            # one-based, so the outer boundary is clamped to one.
            attributes_left = "ID={};NAME={};SIDE=L".format(identifier, family)
            attributes_right = "ID={};NAME={};SIDE=R".format(identifier, family)
            features.append(
                (
                    chromosome,
                    "Liftoff",
                    "repeat_element",
                    max(1, start - flank_size),
                    start,
                    ".",
                    "+",
                    ".",
                    attributes_left,
                )
            )
            features.append(
                (
                    chromosome,
                    "Liftoff",
                    "repeat_element",
                    end,
                    min(lengths[chromosome], end + flank_size),
                    ".",
                    "+",
                    ".",
                    attributes_right,
                )
            )

    Path(output_gff).parent.mkdir(parents=True, exist_ok=True)
    with Path(output_gff).open("w") as handle:
        for feature in features:
            handle.write("\t".join(str(value) for value in feature) + "\n")
    Path(feature_file).write_text("repeat_element\n")
    return len(features) // 2


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--positions", type=Path, required=True)
    parser.add_argument("--genome-index", type=Path, required=True)
    parser.add_argument("--output-gff", type=Path, required=True)
    parser.add_argument("--feature-file", type=Path, required=True)
    parser.add_argument("--flank-size", type=int, default=100000)
    args = parser.parse_args(argv)
    if args.flank_size <= 0:
        parser.error("--flank-size must be greater than zero")
    return args


def main(argv=None):
    args = parse_args(argv)
    count = prepare(
        args.positions,
        args.genome_index,
        args.output_gff,
        args.feature_file,
        args.flank_size,
    )
    print("Prepared Liftoff flanks for {} integrated insertions.".format(count))


if __name__ == "__main__":
    main()
