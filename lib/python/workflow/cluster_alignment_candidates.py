#!/usr/bin/env python3
"""Format bedtools clusters without dropping boundary records or retaining old clusters."""
import argparse
from itertools import groupby
from pathlib import Path


def format_insertions(rows, fasta, positions, sizes, insertion_fasta=None, flank_size=0):
    for cluster_id, cluster_rows in groupby(rows, key=lambda row: row[-1]):
        records = []
        seen = set()
        for row in cluster_rows:
            if row[3] not in seen:
                seen.add(row[3])
                records.append(row)
        first = records[0]
        prefix = f"{first[0]}:<INS>:{first[1]}:{first[2]}:TrEMOLO.INS.{cluster_id}"
        support = len(records)
        for number, row in enumerate(records, 1):
            header = f"{prefix}:{support}:IMPRECISE:{number}"
            fasta.write(f">{header}\n{row[4]}\n")
            if insertion_fasta is not None:
                # count_read is the query position immediately after the CIGAR I.
                # A read can start within flank_size bases of the insertion.
                length = int(row[7])
                left = min(flank_size, int(row[5]) - length)
                sequence = row[4][left:left + length]
                if left < 0 or length < 1 or len(sequence) != length:
                    raise ValueError(f"invalid insertion fragment for {header}")
                insertion_fasta.write(f">{header}\n{sequence}\n")
            sizes.write(f"{header}\t{len(row[4])}\n")
            positions.write(f"{prefix}:{number}:{row[3]}:{row[6]}:{row[7]}:{int(row[6])+int(row[7])}\n")


def format_hard(rows, fasta):
    for cluster_id, cluster_rows in groupby(rows, key=lambda row: row[-1]):
        records = list(cluster_rows)
        first = records[0]
        prefix = f"{first[0]}:<HARD>:{first[1]}:{first[2]}:HARD.{cluster_id}.{first[5]}"
        for number, row in enumerate(records):
            fasta.write(f">{prefix}:{len(records)}:IMPRECISE:{number}\n{row[4]}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("kind", choices=("ins", "hard"))
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--positions", type=Path)
    parser.add_argument("--sizes", type=Path)
    parser.add_argument("--insertion-fasta", type=Path)
    parser.add_argument("--flank-size", type=int, default=0)
    args = parser.parse_args()
    if args.flank_size < 0:
        parser.error("--flank-size must be non-negative")
    with args.input.open() as source, args.fasta.open("w") as fasta:
        rows = (line.rstrip("\n").split("\t") for line in source if line.strip())
        if args.kind == "hard":
            format_hard(rows, fasta)
        else:
            if args.positions is None or args.sizes is None:
                parser.error("ins requires --positions and --sizes")
            with args.positions.open("w") as positions, args.sizes.open("w") as sizes:
                if args.insertion_fasta is None:
                    format_insertions(rows, fasta, positions, sizes)
                else:
                    with args.insertion_fasta.open("w") as insertion_fasta:
                        format_insertions(rows, fasta, positions, sizes, insertion_fasta, args.flank_size)


if __name__ == "__main__":
    main()
