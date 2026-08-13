import argparse
import sys
from multiprocessing import Pool
from typing import NamedTuple, Optional, Sequence

from utils import BM, find_svi_position


MIN_TSD_SIZE = 4
SUPPORTED_BASES = frozenset("ACGTN")


class Candidate(NamedTuple):
    """Validated input needed to find the TSD of one insertion."""

    infos: str
    te: str
    identifier: str
    flank_left: str
    flank_right: str
    genome_flank_left: str
    genome_flank_right: str
    te_size: int
    genome_position: int


def positive_integer(value):
    """Argparse converter accepting strictly positive integers."""

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
        "bed_seq",
        metavar="<bed-seq>",
        type=argparse.FileType("r"),
        help="BED-like file containing the read and genome flanks",
    )
    parser.add_argument(
        "output_file",
        metavar="<output>",
        type=argparse.FileType("w"),
        help="Name of the tabular output file",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=positive_integer,
        default=4,
        help="Number of worker processes (must be greater than zero)",
    )
    return parser


def _unsupported_bases(sequences):
    bases = set("".join(sequences).upper())
    return sorted(bases.difference(SUPPORTED_BASES))


def parse_candidate(line, line_number=None):
    """Parse and validate one seven-column input record.

    The Boyer-Moore implementation used by TrEMOLO supports A, C, G, T and N.
    Other IUPAC ambiguity codes are rejected here so that one such candidate
    cannot terminate the complete multiprocessing run.
    """

    fields = line.strip().split("\t")
    location = "line {}".format(line_number) if line_number is not None else "record"

    if len(fields) != 7:
        return None, "{}: expected 7 tab-separated columns, got {}".format(
            location, len(fields)
        )

    infos, te, identifier, flank_left, flank_right, genome_left, genome_right = fields
    info_fields = infos.split(":")
    if len(info_fields) < 4:
        return None, "{}: malformed insertion metadata in column 1".format(location)

    try:
        te_size = int(info_fields[0])
        genome_position = int(info_fields[3])
    except ValueError:
        return None, (
            "{}: insertion size and genome position in column 1 must be integers".format(
                location
            )
        )

    unsupported = _unsupported_bases(
        (flank_left, flank_right, genome_left, genome_right)
    )
    if unsupported:
        return None, "{}: unsupported ambiguous base(s): {}".format(
            location, ",".join(unsupported)
        )

    candidate = Candidate(
        infos=infos,
        te=te,
        identifier=identifier,
        flank_left=flank_left,
        flank_right=flank_right,
        genome_flank_left=genome_left,
        genome_flank_right=genome_right,
        te_size=te_size,
        genome_position=genome_position,
    )
    return candidate, None


def find_tsd(candidate):
    """Return the historical five-column TSD result for one candidate."""

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

            te_size = candidate.te_size + shift_left + shift_right
            empty_site = candidate.genome_flank_left + candidate.genome_flank_right
            genome_position = find_svi_position(
                tsd,
                empty_site,
                len(candidate.genome_flank_left),
                candidate.genome_position,
            )
            return [
                "{}|{}".format(candidate.te, candidate.identifier),
                tsd,
                str(genome_position),
                str(te_size),
                candidate.infos,
            ]

        motif_size -= 1

    return None


def read_candidates(lines, diagnostics=sys.stderr):
    """Yield valid candidates and report invalid records in input order."""

    for line_number, line in enumerate(lines, 1):
        candidate, error = parse_candidate(line, line_number)
        if error:
            print("Skipping candidate: {}".format(error), file=diagnostics)
            continue
        yield candidate


def call_tsds(lines, threads):
    """Call TSDs while retaining the input order guaranteed by Pool.map."""

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
    results = call_tsds(args.bed_seq, args.threads)
    write_results(results, args.output_file)


if __name__ == "__main__":
    main()
