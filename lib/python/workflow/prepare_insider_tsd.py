#!/usr/bin/env python3
"""Prepare INSIDER genome flanks for legacy-compatible TSD calling."""

from __future__ import print_function

import argparse
from collections import OrderedDict
from pathlib import Path

try:
    import pysam
except ImportError:  # pragma: no cover - exercised only outside the container
    pysam = None


SUPPORTED_BASES = frozenset("ACGTN")


def positive_integer(value):
    try:
        integer = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if integer <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return integer


def read_candidates(path):
    """Read the first four columns of an INSIDER TE BED in input order."""

    candidates = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    "{}:{}: expected at least 4 tab-separated columns".format(
                        path, line_number
                    )
                )
            chromosome, start_text, end_text, name = fields[:4]
            try:
                start = int(start_text)
                end = int(end_text)
            except ValueError:
                raise ValueError(
                    "{}:{}: start and end must be integers".format(path, line_number)
                )
            name_fields = name.split("|", 1)
            if (
                not chromosome
                or start < 0
                or end < start
                or len(name_fields) != 2
                or not all(name_fields)
            ):
                raise ValueError(
                    "{}:{}: invalid INSIDER TE interval".format(path, line_number)
                )
            candidates.append(
                {
                    "chromosome": chromosome,
                    "start": start,
                    "end": end,
                    "start_text": start_text,
                    "end_text": end_text,
                    "name": name,
                    "family": name_fields[0],
                    "identifier": name_fields[1],
                }
            )
    return candidates


def open_indexed_genome(fasta_path, index_path):
    if pysam is None:
        raise RuntimeError("pysam is required to read the indexed genome FASTA")
    return pysam.FastaFile(str(fasta_path), filepath_index=str(index_path))


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def unsupported_bases(*sequences):
    bases = set("".join(sequences).upper())
    return sorted(bases.difference(SUPPORTED_BASES))


def build_outputs(candidates, genome, flank_size):
    """Return the three historical intermediate tables in candidate order."""

    genome_lengths = OrderedDict(zip(genome.references, genome.lengths))
    merged_rows = []
    flank_rows = []
    formatted_rows = []

    for index, candidate in enumerate(candidates, 1):
        chromosome = candidate["chromosome"]
        start = candidate["start"]
        end = candidate["end"]
        name = candidate["name"]
        if chromosome not in genome_lengths:
            raise ValueError(
                "candidate {}: chromosome not found in genome: {}".format(
                    index, chromosome
                )
            )

        left_start = start - flank_size
        right_end = end + flank_size
        if left_start < 0 or right_end > genome_lengths[chromosome]:
            raise ValueError(
                "candidate {}: {}:{}-{} lacks {} bp genome flanks".format(
                    index, chromosome, start, end, flank_size
                )
            )

        left_sequence = genome.fetch(chromosome, left_start, start)
        right_sequence = genome.fetch(chromosome, end, right_end)
        if len(left_sequence) != flank_size or len(right_sequence) != flank_size:
            raise ValueError(
                "candidate {}: indexed FASTA returned incomplete flanks".format(index)
            )
        ambiguous = unsupported_bases(left_sequence, right_sequence)
        if ambiguous:
            raise ValueError(
                "candidate {}: unsupported ambiguous base(s): {}".format(
                    index, ",".join(ambiguous)
                )
            )

        size = end - start
        merged_rows.append(
            (
                chromosome,
                candidate["start_text"],
                candidate["end_text"],
                name,
            )
        )
        flank_rows.extend(
            (
                (chromosome, left_start, start, "{}:{}:FK_L".format(size, name)),
                (chromosome, end, right_end, "{}:{}:FK_R".format(size, name)),
            )
        )
        infos = "{}:{}:{}-{}:{}".format(
            size, chromosome, left_start, start, name
        )
        formatted_rows.append(
            (
                infos,
                candidate["family"],
                candidate["identifier"],
                left_sequence[-flank_size:],
                right_sequence[:flank_size],
                ".",
                ".",
            )
        )

    return merged_rows, flank_rows, formatted_rows


def prepare(args):
    for path in (args.genome, args.genome_index, args.insertion_bed, args.deletion_bed):
        if not path.is_file():
            raise FileNotFoundError(path)

    candidates = read_candidates(args.insertion_bed)
    candidates.extend(read_candidates(args.deletion_bed))
    genome = open_indexed_genome(args.genome, args.genome_index)
    try:
        merged_rows, flank_rows, formatted_rows = build_outputs(
            candidates, genome, args.flank_size
        )
    finally:
        genome.close()

    write_rows(args.merged_bed, merged_rows)
    write_rows(args.flank_bed, flank_rows)
    write_rows(args.formatted_flanks, formatted_rows)
    return len(candidates)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--genome-index", required=True, type=Path)
    parser.add_argument("--insertion-bed", required=True, type=Path)
    parser.add_argument("--deletion-bed", required=True, type=Path)
    parser.add_argument("--flank-size", required=True, type=positive_integer)
    parser.add_argument("--merged-bed", required=True, type=Path)
    parser.add_argument("--flank-bed", required=True, type=Path)
    parser.add_argument("--formatted-flanks", required=True, type=Path)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    count = prepare(args)
    print("Prepared {} INSIDER TSD candidates.".format(count))


if __name__ == "__main__":
    main()
