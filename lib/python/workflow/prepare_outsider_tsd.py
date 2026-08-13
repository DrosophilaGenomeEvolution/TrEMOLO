#!/usr/bin/env python3
"""Prepare OUTSIDER TE sequences and flanks for TSD detection.

This module replaces the historical ``GET_SEQ_TE`` shell block with explicit
inputs and outputs.  The final ``ALL_FK_REPORT_FT*.bed`` files intentionally
retain the legacy layout (including its 11-base left read flank when the
configured flank size is 10), so that the existing TSD caller can be reused.

Two historical selection quirks are also retained for regression
compatibility: the first direct-insertion row matching ``MERGE_TE_ALL.bed`` is
discarded, followed by the first row of the combined Sniffles/direct candidate
list.  They are isolated in :func:`select_legacy_candidates` and can therefore
be removed deliberately in a later data-model migration.
"""

from __future__ import print_function

import argparse
import re
import shutil
from collections import OrderedDict, defaultdict
from pathlib import Path
from urllib.parse import quote

try:
    import pysam
except ImportError:  # pragma: no cover - exercised only outside the container
    pysam = None


CANDIDATE_HEADER = (
    "sseqid",
    "qseqid",
    "pident",
    "size_per",
    "size_el",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "source",
)

VALIDATION_HEADER = (
    "candidate_index",
    "source",
    "family",
    "qseqid",
    "event_id",
    "status",
    "reason",
    "qstart",
    "qend",
    "sequence_length",
    "chromosome",
    "breakpoint",
    "genome_length",
)


def read_fasta(path):
    """Return FASTA records as an ordered ``name -> sequence`` mapping.

    Sequence lines may be wrapped.  Only the first whitespace-delimited token
    of a header is used, matching samtools/bedtools sequence-name semantics.
    """

    records = OrderedDict()
    name = None
    sequence = []
    with path.open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(sequence)
                header = line[1:].strip()
                if not header:
                    raise ValueError("{}:{}: empty FASTA header".format(path, line_number))
                name = header.split()[0]
                sequence = []
            elif name is None:
                raise ValueError(
                    "{}:{}: sequence before first FASTA header".format(path, line_number)
                )
            else:
                sequence.append(line)
    if name is not None:
        records[name] = "".join(sequence)
    return records


def read_rows(path, skip_header=False):
    """Read non-empty tab-delimited rows, optionally dropping one header."""

    rows = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line:
                rows.append(line.split("\t"))
    return rows[1:] if skip_header and rows else rows


def whole_word_match(text, pattern):
    """Implement the relevant ``grep -w`` behaviour used by the old rule."""

    expression = r"(?<![A-Za-z0-9_]){}(?![A-Za-z0-9_])".format(re.escape(pattern))
    return re.search(expression, text) is not None


def event_id(qseqid):
    """Extract a Sniffles/TrEMOLO event identifier from a query name."""

    fields = qseqid.split(":")
    return fields[4] if len(fields) > 4 else ""


def strip_orientation(qseqid):
    """Reproduce legacy ``sed -e 's/:[-+]//g'`` normalization."""

    return re.sub(r":[+-]", "", qseqid)


def select_legacy_candidates(
    sniffles_filtered,
    sniffles_combine,
    merged_bed,
    direct_calls,
    direct_combine,
):
    """Select candidate rows while preserving the two legacy row drops."""

    sniffles_ids = [
        row[1]
        for row in read_rows(sniffles_filtered, skip_header=True)
        if len(row) > 1
    ]

    merged_ids = []
    for row in read_rows(merged_bed):
        if len(row) > 3:
            name_fields = row[3].split("|", 1)
            if len(name_fields) == 2:
                merged_ids.append(name_fields[1])

    # ``grep -w -f merged_ids INS_TREMOLO.csv | awk 'NR>1 {print $2}'``:
    # grep does not select the CSV header, hence NR>1 drops a real candidate.
    matched_direct = []
    for row in read_rows(direct_calls, skip_header=True):
        raw_line = "\t".join(row)
        if any(whole_word_match(raw_line, identifier) for identifier in merged_ids):
            matched_direct.append(row)
    direct_ids = [row[1] for row in matched_direct[1:] if len(row) > 1]

    selected_ids = set(sniffles_ids + direct_ids)
    combined = []
    for source, path in (
        ("sniffles", sniffles_combine),
        ("direct", direct_combine),
    ):
        for row in read_rows(path, skip_header=True):
            if len(row) > 1 and row[1] in selected_ids:
                combined.append((source, row))

    # The generated tmp_combine.csv had no header, but the old awk used NR>1.
    combined = combined[1:]
    # The following shell guard required more than one surviving line before
    # creating tmp_TE_all.csv.  This odd minimum is part of legacy mode too.
    if len(combined) <= 1:
        combined = []
    candidates = []
    for source, row in combined:
        if len(row) < 9:
            # Keep malformed rows visible to validation by padding their fields.
            row = row + [""] * (9 - len(row))
        candidates.append(
            {
                "family": row[0],
                "qseqid": strip_orientation(row[1]),
                "pident": row[2],
                "size_per": row[3],
                "size_el": row[4],
                "qstart": row[5],
                "qend": row[6],
                "sstart": row[7],
                "send": row[8],
                "source": source,
            }
        )
    return candidates


def read_declared_sizes(path):
    """Read the legacy SV query-length table."""

    sizes = OrderedDict()
    for row in read_rows(path):
        if len(row) < 2:
            continue
        try:
            sizes[row[0]] = int(row[1])
        except ValueError:
            continue
    return sizes


def candidate_columns(candidate):
    """Return the historical 13 columns plus an explicit source column."""

    return (
        candidate["family"],
        candidate["qseqid"],
        candidate["pident"],
        candidate["size_per"],
        candidate["size_el"],
        "mismatch",
        "gapopen",
        candidate["qstart"],
        candidate["qend"],
        candidate["sstart"],
        candidate["send"],
        "evalue",
        "bitscore",
        candidate["source"],
    )


def validate_candidate(
    candidate, query_sequences, declared_sizes, genome_lengths, flank_size=10
):
    """Validate one candidate and return fields used by flank extraction."""

    qseqid = candidate["qseqid"]
    sequence = query_sequences.get(qseqid)
    chromosome = ""
    breakpoint_text = ""
    fields = qseqid.split(":")
    if fields:
        chromosome = fields[0]
    if len(fields) > 2:
        breakpoint_text = fields[2]

    context = {
        "sequence": sequence,
        "sequence_length": "" if sequence is None else len(sequence),
        "chromosome": chromosome,
        "breakpoint": breakpoint_text,
        "genome_length": (
            "" if chromosome not in genome_lengths else genome_lengths[chromosome]
        ),
        "qstart": None,
        "qend": None,
        "te_start": None,
        "te_end": None,
    }

    if sequence is None:
        return "rejected", "missing_sequence", context

    try:
        qstart = int(candidate["qstart"])
        qend = int(candidate["qend"])
    except (TypeError, ValueError):
        return "rejected", "invalid_coordinates", context
    context["qstart"] = qstart
    context["qend"] = qend
    if qstart <= 0 or qend <= 0 or qstart == qend:
        return "rejected", "invalid_coordinates", context

    te_start = min(qstart, qend) - 1
    te_end = max(qstart, qend)
    context["te_start"] = te_start
    context["te_end"] = te_end

    # The historical acceptance test used SV_SIZE.tsv.  Prefer it to the
    # physical FASTA length to retain that behaviour, while still preventing
    # extraction outside the actual sequence.
    declared_length = declared_sizes.get(qseqid)
    if declared_length is None:
        declared_length = len(sequence)
    if te_start < 0 or te_end > declared_length or te_end > len(sequence):
        return "rejected", "breakpoint_out_of_bounds", context
    if te_start < 4:
        return "rejected", "insufficient_left_flank", context
    if declared_length - te_end < 4 or len(sequence) - te_end < 4:
        return "rejected", "insufficient_right_flank", context

    if chromosome not in genome_lengths:
        return "rejected", "missing_genome_chromosome", context
    try:
        breakpoint = int(breakpoint_text)
    except (TypeError, ValueError):
        return "rejected", "breakpoint_out_of_bounds", context
    context["breakpoint"] = breakpoint
    genome_length = genome_lengths[chromosome]
    if breakpoint < 1 or breakpoint > genome_length:
        return "rejected", "breakpoint_out_of_bounds", context
    # Compatibility mode accepts truncated genome flanks at contig ends, as
    # bedtools did, but both extracted intervals must contain at least one base.
    genome_left_start = max(0, breakpoint - flank_size - 1)
    genome_left_end = breakpoint - 1
    genome_right_start = breakpoint - 1
    genome_right_end = min(breakpoint + flank_size - 1, genome_length)
    if genome_left_end <= genome_left_start:
        return "rejected", "insufficient_genome_left_flank", context
    if genome_right_end <= genome_right_start:
        return "rejected", "insufficient_genome_right_flank", context

    return "accepted", "accepted", context


def open_indexed_genome(fasta_path, index_path):
    """Open a FASTA through its declared faidx without loading its sequence."""

    if pysam is None:
        raise ImportError(
            "pysam is required to read the indexed genome; run this command "
            "inside the TrEMOLO container"
        )
    return pysam.FastaFile(str(fasta_path), filepath_index=str(index_path))


def family_filename_component(family):
    """Encode an arbitrary family name as one safe path component."""

    return quote(family, safe="._-")


def reset_directory(path):
    """Create an empty output directory without touching neighbouring files."""

    if path.exists():
        if not path.is_dir():
            raise ValueError("output directory path is not a directory: {}".format(path))
        shutil.rmtree(str(path))
    path.mkdir(parents=True)


def write_tsv(path, header, rows):
    """Write a deterministic tab-separated file."""

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        if header:
            handle.write("\t".join(str(value) for value in header) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")


def bedtools_fasta_header(name, sequence_name, start, end):
    """Return the header generated by ``bedtools getfasta -name+``."""

    return "{}::{}:{}-{}".format(name, sequence_name, start, end)


def prepare(args):
    """Generate all explicit GET_SEQ_TE replacement outputs."""

    compatibility_mode = getattr(args, "compatibility_mode", "legacy")
    if compatibility_mode != "legacy":
        raise ValueError("unsupported compatibility mode: {}".format(compatibility_mode))

    input_paths = (
        args.genome,
        args.genome_index,
        args.sniffles_calls,
        args.sniffles_combined,
        args.sniffles_fasta,
        args.merged_bed,
        args.direct_calls,
        args.direct_combined,
        args.direct_fasta,
        args.sequence_sizes,
    )
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.flank_size < 1:
        raise ValueError("--flank-size must be positive")

    genome_fasta = open_indexed_genome(args.genome, args.genome_index)
    try:
        return prepare_with_genome(args, genome_fasta)
    finally:
        genome_fasta.close()


def prepare_with_genome(args, genome_fasta):
    """Generate outputs using an already opened indexed genome handle."""

    genome_lengths = OrderedDict(zip(genome_fasta.references, genome_fasta.lengths))
    query_sequences = read_fasta(args.sniffles_fasta)
    # ``cat sniffles.fasta direct.fasta`` made later duplicate names resolve to
    # the direct record in indexed access; OrderedDict assignment mirrors that.
    for name, sequence in read_fasta(args.direct_fasta).items():
        query_sequences[name] = sequence
    declared_sizes = read_declared_sizes(args.sequence_sizes)

    candidates = select_legacy_candidates(
        args.sniffles_calls,
        args.sniffles_combined,
        args.merged_bed,
        args.direct_calls,
        args.direct_combined,
    )

    write_tsv(
        args.candidates,
        CANDIDATE_HEADER,
        [candidate_columns(candidate) for candidate in candidates],
    )

    size_rows = list(genome_lengths.items())
    size_rows.extend((name, length) for name, length in declared_sizes.items())
    write_tsv(args.all_sequence_sizes, None, size_rows)

    reset_directory(args.te_fasta_dir)
    reset_directory(args.read_artifact_dir)
    reset_directory(args.genome_artifact_dir)

    family_et_records = defaultdict(list)
    family_read_records = defaultdict(list)
    family_genome_records = defaultdict(list)
    validation_rows = []

    for index, candidate in enumerate(candidates, 1):
        status, reason, context = validate_candidate(
            candidate,
            query_sequences,
            declared_sizes,
            genome_lengths,
            args.flank_size,
        )
        qseqid = candidate["qseqid"]
        sequence = context["sequence"]

        # GET_SEQ_TE created TE FASTAs even for rows that later failed the
        # flank checks.  Preserve that useful diagnostic output when the TE
        # interval itself can be represented.
        try:
            raw_qstart = int(candidate["qstart"])
            raw_qend = int(candidate["qend"])
            te_start = min(raw_qstart, raw_qend) - 1
            te_end = max(raw_qstart, raw_qend)
        except (TypeError, ValueError):
            te_start = None
            te_end = None
        if (
            sequence is not None
            and te_start is not None
            and te_end is not None
            and 0 <= te_start < te_end <= len(sequence)
        ):
            family_et_records[candidate["family"]].append(
                (qseqid, te_start, te_end, sequence[te_start:te_end])
            )

        validation_rows.append(
            (
                index,
                candidate["source"],
                candidate["family"],
                qseqid,
                event_id(qseqid),
                status,
                reason,
                candidate["qstart"],
                candidate["qend"],
                context["sequence_length"],
                context["chromosome"],
                context["breakpoint"],
                context["genome_length"],
            )
        )
        if status != "accepted":
            continue

        qstart = context["qstart"]
        qend = context["qend"]
        te_start = context["te_start"]
        te_end = context["te_end"]
        declared_length = declared_sizes.get(qseqid, len(sequence))
        # Legacy BED intervals used the declared query length.  The validation
        # above guarantees that slicing these coordinates is safe.
        read_left_start = 0
        read_left_end = te_start
        read_right_start = te_end
        read_right_end = declared_length
        read_left = sequence[read_left_start:read_left_end]
        read_right = sequence[read_right_start:read_right_end]

        fields = qseqid.split(":")
        chromosome = fields[0]
        breakpoint = int(fields[2])
        genome_left_start = max(0, breakpoint - args.flank_size - 1)
        genome_left_end = breakpoint - 1
        genome_right_start = breakpoint - 1
        genome_right_end = min(
            breakpoint + args.flank_size - 1, genome_lengths[chromosome]
        )
        genome_left = genome_fasta.fetch(
            chromosome, genome_left_start, genome_left_end
        )
        genome_right = genome_fasta.fetch(
            chromosome, genome_right_start, genome_right_end
        )

        family_read_records[candidate["family"]].append(
            {
                "candidate": candidate,
                "candidate_key": (index, candidate["source"], candidate["family"], qseqid),
                "left_start": read_left_start,
                "left_end": read_left_end,
                "right_start": read_right_start,
                "right_end": read_right_end,
                "left": read_left,
                "right": read_right,
            }
        )
        family_genome_records[candidate["family"]].append(
            {
                "candidate": candidate,
                "candidate_key": (index, candidate["source"], candidate["family"], qseqid),
                "chromosome": chromosome,
                "left_start": genome_left_start,
                "left_end": genome_left_end,
                "right_start": genome_right_start,
                "right_end": genome_right_end,
                "left": genome_left,
                "right": genome_right,
            }
        )

    write_tsv(args.validation, VALIDATION_HEADER, validation_rows)

    for family in sorted(family_et_records):
        path = args.te_fasta_dir / "TE_REPORT_FOUND_{}.fasta".format(
            family_filename_component(family)
        )
        with path.open("w") as handle:
            for qseqid, start, end, sequence in family_et_records[family]:
                handle.write(">{}:{}-{}\n{}\n".format(qseqid, start, end, sequence))

    all_read_entries = []
    all_genome_entries = []
    for family in sorted(family_read_records):
        read_bed_rows = []
        read_report_rows = []
        genome_bed_rows = []
        genome_report_rows = []

        for record in family_read_records[family]:
            candidate = record["candidate"]
            qseqid = candidate["qseqid"]
            left_name = "{}:{}:FKL".format(candidate["size_el"], family)
            right_name = "{}:{}:FKR".format(candidate["size_el"], family)
            read_bed_rows.extend(
                (
                    (qseqid, record["left_start"], record["left_end"], left_name),
                    (qseqid, record["right_start"], record["right_end"], right_name),
                )
            )
            read_report_rows.extend(
                (
                    (
                        bedtools_fasta_header(
                            left_name, qseqid, record["left_start"], record["left_end"]
                        ),
                        record["left"],
                    ),
                    (
                        bedtools_fasta_header(
                            right_name,
                            qseqid,
                            record["right_start"],
                            record["right_end"],
                        ),
                        record["right"],
                    ),
                )
            )
            # Historical awk started one character earlier than a conventional
            # N-base suffix, producing flank_size + 1 bases when available.
            read_left = record["left"][-(args.flank_size + 1) :]
            read_right = record["right"][: args.flank_size]
            all_read_entries.append(
                (
                    record["candidate_key"],
                    (
                    "{}:{}".format(candidate["size_el"], qseqid),
                    family,
                    event_id(qseqid),
                    read_left,
                    read_right,
                    ),
                )
            )

        for record in family_genome_records[family]:
            candidate = record["candidate"]
            qseqid = candidate["qseqid"]
            chromosome = record["chromosome"]
            genome_bed_rows.extend(
                (
                    (
                        chromosome,
                        record["left_start"],
                        record["left_end"],
                        qseqid,
                        family + ":FKL",
                    ),
                    (
                        chromosome,
                        record["right_start"],
                        record["right_end"],
                        qseqid,
                        family + ":FKR",
                    ),
                )
            )
            genome_report_rows.extend(
                (
                    (
                        bedtools_fasta_header(
                            qseqid,
                            chromosome,
                            record["left_start"],
                            record["left_end"],
                        ),
                        record["left"],
                    ),
                    (
                        bedtools_fasta_header(
                            qseqid,
                            chromosome,
                            record["right_start"],
                            record["right_end"],
                        ),
                        record["right"],
                    ),
                )
            )
            # The odd leading colon and embedded left-region coordinates are
            # part of ALL_FK_REPORT_FT_GENOME.bed's historical first column.
            genome_info = ":{}::{}:{}-{}".format(
                qseqid, chromosome, record["left_start"], record["left_end"]
            )
            all_genome_entries.append(
                (
                    record["candidate_key"],
                    (
                    genome_info,
                    chromosome,
                    event_id(qseqid),
                    record["left"][-args.flank_size :],
                    record["right"][: args.flank_size],
                    ),
                )
            )

        family_component = family_filename_component(family)
        write_tsv(
            args.read_artifact_dir
            / "TE_REPORT_flanking_{}.bed".format(family_component),
            None,
            read_bed_rows,
        )
        write_tsv(
            args.read_artifact_dir
            / "FK_REPORT_FOUND_{}.bed".format(family_component),
            None,
            read_report_rows,
        )
        write_tsv(
            args.genome_artifact_dir
            / "TE_REPORT_SVI_{}.bed".format(family_component),
            None,
            genome_bed_rows,
        )
        write_tsv(
            args.genome_artifact_dir
            / "FK_REPORT_FOUND_SVI_{}.bed".format(family_component),
            None,
            genome_report_rows,
        )

    write_tsv(args.read_flanks, None, [row for _, row in all_read_entries])
    write_tsv(args.genome_flanks, None, [row for _, row in all_genome_entries])
    genome_by_candidate = {
        key: (row[3], row[4])
        for key, row in all_genome_entries
        if row[3] and row[4]
    }
    all_rows = []
    for key, row in all_read_entries:
        flanks = genome_by_candidate.get(key)
        if flanks:
            all_rows.append(row + flanks)
    write_tsv(args.combined_flanks, None, all_rows)

    return len(candidates), len(all_rows)


def build_parser():
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", required=True, type=Path)
    parser.add_argument("--genome-index", required=True, type=Path)
    parser.add_argument("--sniffles-calls", required=True, type=Path)
    parser.add_argument("--sniffles-combined", required=True, type=Path)
    parser.add_argument("--sniffles-fasta", required=True, type=Path)
    parser.add_argument("--merged-bed", required=True, type=Path)
    parser.add_argument("--direct-calls", required=True, type=Path)
    parser.add_argument("--direct-combined", required=True, type=Path)
    parser.add_argument("--direct-fasta", required=True, type=Path)
    parser.add_argument("--sequence-sizes", required=True, type=Path)
    parser.add_argument("--flank-size", required=True, type=int)
    parser.add_argument(
        "--compatibility-mode", choices=("legacy",), default="legacy"
    )
    parser.add_argument("--te-fasta-dir", required=True, type=Path)
    parser.add_argument("--read-artifact-dir", required=True, type=Path)
    parser.add_argument("--genome-artifact-dir", required=True, type=Path)
    parser.add_argument("--candidates", required=True, type=Path)
    parser.add_argument("--validation", required=True, type=Path)
    parser.add_argument("--all-sequence-sizes", required=True, type=Path)
    parser.add_argument("--read-flanks", required=True, type=Path)
    parser.add_argument("--genome-flanks", required=True, type=Path)
    parser.add_argument("--combined-flanks", required=True, type=Path)
    return parser


def main():
    """Run the standalone preparation command."""

    args = build_parser().parse_args()
    candidate_count, accepted_count = prepare(args)
    print(
        "Prepared {} OUTSIDER TSD candidates ({} accepted)".format(
            candidate_count, accepted_count
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
