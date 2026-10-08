"""Extract CIGAR insertions and clipped evidence with deterministic, balanced BAM jobs."""
import argparse
import logging
import multiprocessing
import os
from pathlib import Path
import shutil
import sys
import tempfile
import time

import pysam


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('bam_file')
    parser.add_argument('--output-soft', default='SOFT.txt')
    parser.add_argument('--output-ins', default='INS.txt')
    parser.add_argument('--output-hard', default='HARD.txt')
    parser.add_argument('--output-seq-tsd')
    parser.add_argument('--no-clipped', action='store_true')
    parser.add_argument('--time-limit', type=int, default=0, help='Maximum elapsed hours; 0 disables the limit')
    parser.add_argument('-t', '--threads', type=int, default=4)
    parser.add_argument('-f', '--flank-size', type=int, default=30)
    parser.add_argument('-w', '--window', type=int, default=50)
    parser.add_argument('-s', '--min-size', type=int, default=30)
    parser.add_argument('--chunk-size', type=int, default=10000000, help='Reference bases per extraction job; 0 uses whole chromosomes')
    args = parser.parse_args()
    if args.threads < 1 or min(args.window, args.time_limit, args.flank_size, args.chunk_size) < 0 or args.min_size < 1:
        parser.error('threads and minimum size must be positive; window, time, flank and chunk size must be non-negative')
    return args


def initialize_worker(options, stop_event, directory):
    global worker_options, worker_stop_event, worker_directory
    worker_options, worker_stop_event, worker_directory = options, stop_event, Path(directory)


def write_row(output, fields):
    output.write('\t'.join(map(str, fields)) + '\n')


def process_region(task):
    task_id, chrom, start, end = task
    args = worker_options
    counts = [0, 0, 0]
    with pysam.AlignmentFile(args.bam_file, 'rb') as bam, \
            (worker_directory / f'ins.{task_id}').open('w') as ins, \
            (worker_directory / f'tsd.{task_id}').open('w') as tsd, \
            (worker_directory / f'soft.{task_id}').open('w') as soft, \
            (worker_directory / f'hard.{task_id}').open('w') as hard:
        for number, read in enumerate(bam.fetch(chrom, start, end)):
            # fetch also returns alignments overlapping from a preceding block.
            # Each alignment belongs exclusively to the block containing its start.
            if not start <= read.reference_start < end:
                continue
            if args.time_limit and number % 1024 == 0 and worker_stop_event.is_set():
                raise TimeoutError('Alignment candidate extraction exceeded TIME_LIMIT')
            sequence = read.query_sequence
            name = read.query_name
            if not name:
                continue
            count_ref = count_read = count_real = 0
            hard_seen = set()
            for operation, length in read.cigartuples or ():
                # Preserve the historical CIGAR-coordinate contract in this change.
                if operation in (0, 2, 7):
                    count_ref += length
                if operation in (0, 1, 7, 4):
                    count_read += length
                if operation in (0, 1, 7, 4, 5):
                    count_real += length
                position = read.reference_start + count_ref
                if operation == 1 and length >= args.min_size and sequence:
                    left = max(count_read-length-args.flank_size, 0)
                    right = min(count_read+args.flank_size, len(sequence))
                    fragment = sequence[left:right]
                    write_row(ins, [chrom, position, position+1, name, fragment, count_read, count_real, length])
                    if args.output_seq_tsd:
                        flanks = sequence[left:left+args.flank_size] + '|' + fragment + '|' + sequence[right-args.flank_size:right]
                        write_row(tsd, [chrom, position, name, flanks, fragment, count_read, count_real, length])
                    counts[0] += 1
                if args.no_clipped or operation not in (4, 5) or length <= args.min_size:
                    continue
                if operation == 5 and sequence is None:
                    raise ValueError(f'Hard-clipped alignment {name} has no sequence')
                if not sequence:
                    continue
                if operation == 5:
                    identity = (position, name, count_real, len(sequence))
                    if identity in hard_seen:
                        continue
                    hard_seen.add(identity)
                    side = 'L' if count_real-length == 0 else 'R'
                    real_position = count_real if side == 'L' else count_real-length
                    write_row(hard, [chrom, position, position+1, '.', name, real_position, read.flag, side, len(sequence), sequence])
                    counts[2] += 1
                else:
                    side = 'L' if count_read-length == 0 else 'R'
                    real_position = count_real if side == 'L' else count_real-length
                    clipped = sequence[count_read-length:count_read]
                    write_row(soft, [chrom, position, position+1, '.', name, real_position, read.flag, side, len(clipped), clipped])
                    counts[1] += 1
    logging.info('%s:%s-%s: INS=%s SOFT=%s HARD=%s', chrom, start, end, *counts)
    return counts


def merge_soft(parts, output, all_output, chrom_index, window):
    # Group after ordered block merging, so boundary reads and cluster IDs are
    # identical for any chunk size or worker count. Buckets avoid global scans.
    clusters = []
    buckets = {}
    for part in parts:
        with part.open() as source:
            for line in source:
                row = line.rstrip('\n').split('\t')
                position = int(row[1])
                bucket = position // max(window, 1)
                matches = [index for key in (bucket-1, bucket, bucket+1)
                           for index in buckets.get(key, ())
                           if abs(clusters[index]['position'] - position) <= window]
                if matches:
                    index = min(matches)
                    cluster = clusters[index]
                else:
                    index = len(clusters)
                    cluster = {'chrom': row[0], 'position': position, 'id': f'SOFT.{chrom_index}_{index}',
                               'L': ['NONE', '0', 'NONE'], 'R': ['NONE', '0', 'NONE'],
                               'reads_L': {}, 'reads_R': {}, 'support': 0}
                    clusters.append(cluster)
                    buckets.setdefault(bucket, []).append(index)
                side = row[7]
                if int(row[8]) > int(cluster[side][1]):
                    cluster[side] = [row[4], row[8], row[9]]
                cluster['reads_' + side][row[4] + ':' + row[9]] = None
                cluster['support'] += 1
                row[3] = cluster['id']
                write_row(all_output, row)
    for cluster in clusters:
        left, right = cluster['L'], cluster['R']
        write_row(output, [cluster['chrom'], cluster['position'], cluster['id'],
                          f'BEST_L_RS={left[0]};BEST_L_SIZE={left[1]};BEST_L_SEQ={left[2]}',
                          f'BEST_R_RS={right[0]};BEST_R_SIZE={right[1]};BEST_R_SEQ={right[2]}',
                          'RS_LEFT=' + ','.join(cluster['reads_L']),
                          'RS_RIGHT=' + ','.join(cluster['reads_R']),
                          'NB_RS=' + str(cluster['support'])])


def main():
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    with pysam.AlignmentFile(args.bam_file, 'rb') as bam:
        chromosomes = list(zip(bam.references, bam.lengths))
        mapped = {entry.contig: entry.mapped for entry in bam.get_index_statistics()}
    tasks = []
    chromosome_tasks = []
    for chrom, length in chromosomes:
        current = []
        if mapped.get(chrom, 0):
            width = args.chunk_size or length
            for start in range(0, length, width):
                task = (len(tasks), chrom, start, min(start+width, length))
                tasks.append(task)
                current.append(task)
        chromosome_tasks.append(current)
    logging.info('chromosomes: %s; jobs: %s; --no-clipped: %s', len(chromosomes), len(tasks), args.no_clipped)
    stop_event = multiprocessing.Event()
    # Keep final-order IDs separate from scheduling. Larger expected jobs start first.
    lengths = dict(chromosomes)
    scheduled = sorted(tasks, key=lambda task: -mapped[task[1]] * (task[3]-task[2]) / lengths[task[1]])
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='.alignment-candidates-', dir=str(Path(args.output_ins).resolve().parent)) as directory:
        scratch = Path(directory)
        if tasks:
            with multiprocessing.Pool(processes=min(args.threads, len(tasks)),
                                      initializer=initialize_worker, initargs=(args, stop_event, directory)) as pool:
                result = pool.map_async(process_region, scheduled, chunksize=1)
                try:
                    result.get(timeout=args.time_limit*3600 if args.time_limit else None)
                except multiprocessing.TimeoutError:
                    stop_event.set()
                    raise TimeoutError('Alignment candidate extraction exceeded TIME_LIMIT')
        output_paths = {'ins': args.output_ins, 'soft': args.output_soft,
                        'soft_all': args.output_soft + '.bis', 'hard': args.output_hard}
        if args.output_seq_tsd:
            output_paths['tsd'] = args.output_seq_tsd
        from contextlib import ExitStack
        with ExitStack() as stack:
            outputs = {key: stack.enter_context((scratch / ('final.' + key)).open('w')) for key in output_paths}
            if not args.no_clipped:
                write_row(outputs['soft_all'], ['#REF', 'START', 'END', 'ID', 'ID_READ', 'POS_REAL_READ', 'FLAG', 'SIDE', 'SIZE_SEQ', 'SEQ'])
                write_row(outputs['hard'], ['#REF', 'START', 'END', 'ID', 'ID_READ', 'POS_REAL_READ', 'FLAG', 'SIDE', 'SIZE_SEQ', 'SEQ'])
            for chrom_index, current in enumerate(chromosome_tasks):
                for task in current:
                    for key in ('ins', 'tsd'):
                        if key in outputs:
                            with (scratch / f'{key}.{task[0]}').open() as source:
                                shutil.copyfileobj(source, outputs[key])
                if not args.no_clipped:
                    merge_soft([scratch / f'soft.{task[0]}' for task in current], outputs['soft'], outputs['soft_all'], chrom_index, args.window)
                    hard_index = 0
                    for task in current:
                        with (scratch / f'hard.{task[0]}').open() as source:
                            for line in source:
                                row = line.rstrip('\n').split('\t')
                                row[3] = f'HARD.{chrom_index}_{hard_index}'
                                write_row(outputs['hard'], row)
                                hard_index += 1
                if args.time_limit and time.monotonic()-started >= args.time_limit*3600:
                    raise TimeoutError('Alignment candidate extraction exceeded TIME_LIMIT')
        # Outputs are published only after every worker and merge has succeeded.
        for key, destination in output_paths.items():
            shutil.move(str(scratch / ('final.' + key)), destination)
    logging.info('Candidate extraction completed in %.3f seconds', time.monotonic()-started)


if __name__ == '__main__':
    main()
