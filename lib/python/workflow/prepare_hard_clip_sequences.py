#!/usr/bin/env python3
"""Extract HARD query fragments in one FASTQ pass, retaining legacy audit artifacts."""
import argparse
from pathlib import Path
import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=Path, required=True)
    parser.add_argument('--fastq', type=Path, required=True)
    parser.add_argument('--reads-fasta', type=Path, required=True)
    parser.add_argument('--reads-index', type=Path, required=True)
    parser.add_argument('--variants', type=Path, required=True)
    parser.add_argument('--max-size', type=int, required=True)
    args = parser.parse_args()
    candidates = []
    by_read = {}
    with args.candidates.open() as source:
        for line in source:
            if line.startswith('#') or not line.strip():
                continue
            row = line.rstrip('\n').split('\t')
            candidates.append(row)
            by_read.setdefault(row[4], []).append(len(candidates)-1)
    fragments = {}
    with args.reads_fasta.open('w') as fasta, args.reads_index.open('w') as index:
        if args.fastq.stat().st_size:
            with pysam.FastxFile(str(args.fastq)) as reads:
                for number, read in enumerate(reads):
                    fasta.write(f'>{read.name}\n{read.sequence}\n')
                    index.write(f'{2*number+1}:{read.name}\n')
                    for candidate_index in by_read.get(read.name, ()):
                        row = candidates[candidate_index]
                        size = int(row[5])
                        side = row[7]
                        # Keep historical slice and reported-length semantics.
                        sequence = read.sequence[:size] if side == 'L' else read.sequence[-size:]
                        if size > args.max_size:
                            sequence = sequence[-args.max_size:] if side == 'L' else sequence[:args.max_size]
                        fragments[candidate_index] = (sequence, min(size, args.max_size))
    with args.variants.open('w') as output:
        for number, row in enumerate(candidates):
            if number not in fragments:
                raise ValueError(f'Hard-clipped read {row[4]} is missing from extracted FASTQ')
            sequence, size = fragments[number]
            output.write('\t'.join([row[0], row[1], row[2], row[4], sequence, row[7], str(size)]) + '\n')


if __name__ == '__main__':
    main()
