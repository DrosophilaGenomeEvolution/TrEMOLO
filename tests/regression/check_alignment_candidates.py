#!/usr/bin/env python3
"""Validate corrected CIGAR evidence and complete cluster emission on any work directory."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path
import pysam


def fasta_sequences(path):
    records = {}
    header = None
    with path.open() as source:
        for line in source:
            if line.startswith('>'):
                header = line[1:].strip()
                if header in records:
                    raise ValueError(f'Duplicate FASTA ID: {header}')
                records[header] = ''
            elif header is not None:
                records[header] += line.strip()
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work_directory', type=Path)
    parser.add_argument('--no-clipped', action='store_true')
    parser.add_argument('--min-size', type=int, default=30)
    args = parser.parse_args()
    root = args.work_directory/'OUTSIDER/TrEMOLO_SV_TE'
    clustered = defaultdict(dict)
    with (root/'INS/SV_INS_CLUST.bed').open() as source:
        for line in source:
            row = line.rstrip('\n').split('\t')
            clustered[row[-1]].setdefault(row[3], row)
    expected_fasta = {}
    expected_positions = []
    for cluster, unique_reads in clustered.items():
        rows = list(unique_reads.values())
        first = rows[0]
        prefix = f'{first[0]}:<INS>:{first[1]}:{first[2]}:TrEMOLO.INS.{cluster}'
        for number, row in enumerate(rows, 1):
            expected_fasta[f'{prefix}:{len(rows)}:IMPRECISE:{number}'] = row[4]
            expected_positions.append(f'{prefix}:{number}:{row[3]}:{row[6]}:{row[7]}:{int(row[6])+int(row[7])}')
    assert fasta_sequences(root/'INS/SV_INS_CLUST.fasta') == expected_fasta, 'INS cluster emission lost or changed candidates'
    assert (root/'INS/RD_NUMBER.txt').read_text().splitlines() == expected_positions, 'Incomplete INS read-position provenance'
    expected_soft = Counter()
    expected_hard = Counter()
    expected_ins = 0
    with pysam.AlignmentFile(str(args.work_directory/'OUTSIDER/MAPPING/SAMPLE_mapping_GENOME_MD.sorted.bam'), 'rb') as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or not read.query_name:
                continue
            ref = query = real = 0
            seen_hard = set()
            for op, length in read.cigartuples or ():
                if op in (0, 2, 7): ref += length
                if op in (0, 1, 7, 4): query += length
                if op in (0, 1, 7, 4, 5): real += length
                if op == 1 and length >= args.min_size and read.query_sequence: expected_ins += 1
                if args.no_clipped or op not in (4, 5) or length <= args.min_size or not read.query_sequence: continue
                pos = read.reference_start+ref
                if op == 4:
                    side = 'L' if query-length == 0 else 'R'
                    real_pos = real if side == 'L' else real-length
                    expected_soft[(read.reference_name, pos, read.query_name, real_pos, read.flag, side, read.query_sequence[query-length:query])] += 1
                else:
                    identity = (pos, read.query_name, real, len(read.query_sequence))
                    if identity in seen_hard: continue
                    seen_hard.add(identity)
                    side = 'L' if real-length == 0 else 'R'
                    real_pos = real if side == 'L' else real-length
                    expected_hard[(read.reference_name, pos, read.query_name, real_pos, read.flag, side, read.query_sequence)] += 1
    assert len((root/'INS/SV_INS.bed').read_text().splitlines()) == expected_ins, 'CIGAR insertion count mismatch'
    for kind, path, expected in [('SOFT', root/'SOFT/SV_SOFT.vcf.bis', expected_soft), ('HARD', root/'HARD/SV_HARD.tr_vcf', expected_hard)]:
        actual = Counter()
        with path.open() as source:
            for line in source:
                if line.startswith('#') or not line.strip(): continue
                row = line.rstrip('\n').split('\t')
                actual[(row[0], int(row[1]), row[4], int(row[5]), int(row[6]), row[7], row[9])] += 1
        assert actual == expected, f'{kind} evidence differs from BAM provenance'
    hard_sequences = Counter()
    with (root/'HARD/HARD.bed').open() as source:
        for line in source:
            if line.strip(): hard_sequences[line.rstrip('\n').split('\t')[4]] += 1
    assert Counter(fasta_sequences(root/'HARD/HARD.fasta').values()) == hard_sequences, 'HARD cluster emission lost candidates'
    print(f'Corrected alignment candidates: PASSED ({expected_ins} INS records, {sum(expected_soft.values())} SOFT records, {sum(expected_hard.values())} HARD records; all INS/HARD clusters retained)')


if __name__ == '__main__':
    main()
