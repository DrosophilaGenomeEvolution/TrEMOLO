import argparse
import sys
from multiprocessing import Pool
from typing import NamedTuple, Optional, Sequence

from utils import BM


MIN_TSD_SIZE = 4
SUPPORTED_BASES = frozenset("ACGTN")


class Candidate(NamedTuple):
    infos: str
    te: str
    identifier: str
    flank_left: str
    flank_right: str
    te_size: int
    genome_position: int


def positive_integer(value):
    try:
        integer = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be an integer")
    if integer <= 0:
        raise argparse.ArgumentTypeError("must be greater than zero")
    return integer


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "bed_seq", metavar="<bed-seq>", type=argparse.FileType("r"),
        help="BED-like file containing paired genome flanks",
    )
    parser.add_argument(
        "output_file", metavar="<output>", type=argparse.FileType("w"),
        help="Name of the tabular output file",
    )
    parser.add_argument(
        "-t", "--threads", dest="threads", type=positive_integer, default=4,
        help="Number of worker processes (must be greater than zero)",
    )
    return parser


def parse_candidate(line, line_number=None):
    fields = line.strip().split("\t")
    location = "line {}".format(line_number) if line_number else "record"
    if len(fields) != 7:
        return None, "{}: expected 7 tab-separated columns, got {}".format(
            location, len(fields)
        )

    infos, te, identifier, flank_left, flank_right, _genome_left, _genome_right = fields
    info_fields = infos.split(":")
    if len(info_fields) < 4 or "-" not in info_fields[2]:
        return None, "{}: malformed INSIDER flank metadata".format(location)
    try:
        te_size = int(info_fields[0])
        genome_position = int(info_fields[2].split("-", 1)[0])
    except ValueError:
        return None, "{}: TE size and genome position must be integers".format(
            location
        )

    bases = set((flank_left + flank_right).upper())
    unsupported = sorted(bases.difference(SUPPORTED_BASES))
    if unsupported:
        return None, "{}: unsupported ambiguous base(s): {}".format(
            location, ",".join(unsupported)
        )

    return Candidate(
        infos=infos,
        te=te,
        identifier=identifier,
        flank_left=flank_left,
        flank_right=flank_right,
        te_size=te_size,
        genome_position=genome_position,
    ), None


def find_tsd(candidate):
    flank_left = candidate.flank_left
    flank_right = candidate.flank_right
    shorter_flank = flank_right
    other_flank = flank_left
    mode = "R"
    if len(flank_right) > len(flank_left):
        mode = "L"
        shorter_flank = flank_left
        other_flank = flank_right

    motif_size = len(shorter_flank)
    while motif_size >= MIN_TSD_SIZE:
        for offset in range(len(shorter_flank) - motif_size + 1):
            tsd = shorter_flank[offset : offset + motif_size]
            motif_positions = BM(tsd, other_flank)
            if not motif_positions:
                continue
            if mode == "R":
                shift_left = len(flank_left) - motif_size - max(motif_positions)
                shift_right = offset
            else:
                shift_left = len(flank_left) - motif_size - offset
                shift_right = min(motif_positions)

            return [
                "{}|{}".format(candidate.te, candidate.identifier),
                tsd,
                str(candidate.genome_position - shift_left),
                str(candidate.te_size + shift_left + shift_right),
                candidate.infos,
            ]
        motif_size -= 1
    return None


def read_candidates(lines, diagnostics=sys.stderr):
    for line_number, line in enumerate(lines, 1):
        candidate, error = parse_candidate(line, line_number)
        if error:
            print("Skipping candidate: {}".format(error), file=diagnostics)
            continue
        yield candidate


def call_tsds(lines, threads):
    candidates = list(read_candidates(lines))
    if not candidates:
        return []
    with Pool(processes=threads) as pool:
        return pool.map(find_tsd, candidates)


def write_results(results, output_file):
    for result in results:
        if result is not None:
            output_file.write("\t".join(result) + "\n")


def main(argv: Optional[Sequence[str]] = None):
    args = build_parser().parse_args(argv)
    write_results(call_tsds(args.bed_seq, args.threads), args.output_file)


if __name__ == "__main__":
    main()
