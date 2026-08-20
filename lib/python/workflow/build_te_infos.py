#!/usr/bin/env python3
"""Build the legacy 15-column TrEMOLO ``TE_INFOS.bed`` table."""

from __future__ import print_function

import argparse
import ctypes
import os
import re
from collections import OrderedDict, defaultdict
from pathlib import Path


HEADER = (
    "#chrom",
    "start",
    "end",
    "TE|ID",
    "strand",
    "TSD",
    "pident",
    "psize_TE",
    "SIZE_TE",
    "NEW_POS",
    "FREQ",
    "FREQ_WITH_CLIPPED",
    "SV_SIZE",
    "ID_TrEMOLO",
    "TYPE",
)

CONTRACTION_TYPES = frozenset(
    ("Deletion", "Repeat_contraction", "Tandem_contraction")
)
EXPANSION_TYPES = frozenset(
    ("Insertion", "Repeat_expansion", "Tandem_expansion")
)


def read_rows(path, skip_header=False):
    rows = []
    with path.open() as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\r\n")
            if line:
                rows.append(line.split("\t"))
    return rows[1:] if skip_header and rows else rows


def require_width(row, minimum, source):
    if len(row) < minimum:
        raise ValueError(
            "{}: expected at least {} columns, got {}".format(
                source, minimum, len(row)
            )
        )


def outsider_event_id(qseqid):
    fields = qseqid.split(":")
    return fields[4] if len(fields) > 4 else ""


def insider_event_id(qseqid):
    return qseqid.split(":", 1)[0]


def first_by_event(rows, event_parser, source):
    indexed = OrderedDict()
    for row in rows:
        require_width(row, 2, source)
        identifier = event_parser(row[1])
        if not identifier:
            raise ValueError("{}: malformed qseqid: {}".format(source, row[1]))
        indexed.setdefault(identifier, row)
    return indexed


def tsd_by_event(rows, source):
    indexed = OrderedDict()
    for row in rows:
        require_width(row, 5, source)
        name = row[0].split("|", 1)
        if len(name) != 2 or not all(name):
            raise ValueError("{}: malformed TE|ID: {}".format(source, row[0]))
        indexed.setdefault(name[1], row)
    return indexed


def legacy_chain_id(chain):
    """Reproduce the float32 arithmetic of ``chain_to_id.cpp`` exactly."""

    value = ctypes.c_float(0.0).value
    for index, character in enumerate(chain):
        code = ord(character)
        numerator = ctypes.c_float(float(7 + index)).value
        denominator = ctypes.c_float(float(5 - code)).value
        fraction = ctypes.c_float(numerator / denominator).value
        term = ctypes.c_float(ctypes.c_float(float(code)).value + fraction).value
        value = ctypes.c_float(value + term).value

    integer = int(value)
    fraction = ctypes.c_float(value - integer).value
    value = ctypes.c_float(
        ctypes.c_float(fraction * ctypes.c_float(100000.0).value).value + integer
    ).value
    return format(value, ".0f")


def parse_position(row, source, minimum_columns):
    require_width(row, minimum_columns, source)
    name = row[3].split("|", 1)
    if len(name) != 2 or not all(name):
        raise ValueError("{}: malformed TE|ID: {}".format(source, row[3]))
    try:
        start = int(row[1])
        end = int(row[2])
    except ValueError:
        raise ValueError("{}: start and end must be integers".format(source))
    strand = row[5] if minimum_columns >= 6 else row[4]
    if not row[0] or start < 0 or end < start or strand not in ("+", "-"):
        raise ValueError("{}: invalid genomic interval".format(source))
    return {
        "chromosome": row[0],
        "start": row[1],
        "end": row[2],
        "start_value": start,
        "end_value": end,
        "family": name[0],
        "identifier": name[1],
        "strand": strand,
    }


def build_deletion_windows(rows):
    by_chromosome = defaultdict(list)
    for row in rows:
        position = parse_position(row, "INSIDER deletion BED", 5)
        by_chromosome[position["chromosome"]].append(position)
    return by_chromosome


def intersects_same_family_deletion(position, deletions, window=30):
    for deletion in deletions.get(position["chromosome"], ()):
        if deletion["family"] != position["family"]:
            continue
        if (
            deletion["start_value"] < position["end_value"] + window
            and deletion["end_value"] > position["start_value"] - window
        ):
            return True
    return False


def index_frequency(rows, source):
    indexed = OrderedDict()
    for row in rows:
        require_width(row, 10, source)
        indexed.setdefault(row[1], row)
    return indexed


def index_sizes(rows, source):
    indexed = defaultdict(list)
    for row in rows:
        require_width(row, 2, source)
        identifier = outsider_event_id(row[0])
        if identifier:
            indexed[identifier].append(row[1])
    return indexed


def maximum_size(values):
    maximum = None
    maximum_text = "NONE"
    for value in values:
        try:
            numeric = float(value)
        except ValueError:
            continue
        if maximum is None or numeric > maximum:
            maximum = numeric
            maximum_text = value
    return maximum_text


def valid_percentage(value):
    try:
        return float(value) <= 100
    except ValueError:
        return False


def outsider_frequency_values(identifier, precise, raw):
    precise_row = precise.get(identifier)
    if (
        precise_row is not None
        and valid_percentage(precise_row[8])
        and valid_percentage(precise_row[9])
    ):
        return precise_row[8], precise_row[9]
    raw_row = raw.get(identifier)
    if raw_row is not None:
        return raw_row[8], raw_row[9]
    return "NONE", "NONE"


def outsider_variant_type(identifier, sniffles_calls, direct_calls):
    row = sniffles_calls.get(identifier) or direct_calls.get(identifier)
    if row is None:
        return ""
    fields = row[1].split(":")
    if len(fields) < 2:
        return ""
    return fields[1].replace("<", "").replace(">", "")


def render_outsider_rows(args):
    positions = read_rows(args.outsider_positions)
    source_tables = [
        first_by_event(
            read_rows(args.outsider_sniffles_combine, skip_header=True),
            outsider_event_id,
            "OUTSIDER Sniffles combined calls",
        ),
        first_by_event(
            read_rows(args.outsider_direct_combine, skip_header=True),
            outsider_event_id,
            "OUTSIDER direct combined calls",
        ),
        first_by_event(
            read_rows(args.outsider_soft_calls, skip_header=True),
            outsider_event_id,
            "OUTSIDER soft-clipped calls",
        ),
        first_by_event(
            read_rows(args.outsider_hard_calls, skip_header=True),
            outsider_event_id,
            "OUTSIDER hard-clipped calls",
        ),
    ]
    tsds = tsd_by_event(read_rows(args.outsider_tsd), "OUTSIDER TSD")
    precise = index_frequency(
        read_rows(args.outsider_frequency_precise), "OUTSIDER precise frequency"
    )
    raw_frequency = index_frequency(
        read_rows(args.outsider_frequency), "OUTSIDER frequency"
    )
    sizes = index_sizes(read_rows(args.outsider_sv_sizes), "OUTSIDER SV sizes")
    sniffles_calls = first_by_event(
        read_rows(args.outsider_sniffles_calls, skip_header=True),
        outsider_event_id,
        "OUTSIDER Sniffles calls",
    )
    direct_calls = first_by_event(
        read_rows(args.outsider_direct_calls, skip_header=True),
        outsider_event_id,
        "OUTSIDER direct calls",
    )
    deletions = build_deletion_windows(read_rows(args.insider_deletion_bed))

    rendered = []
    for raw_position in positions:
        position = parse_position(raw_position, "OUTSIDER positions", 6)
        identifier = position["identifier"]
        source_row = None
        for table in source_tables:
            if identifier in table:
                source_row = table[identifier]
                break
        pident = source_row[2] if source_row is not None else "NONE"
        size_percent = source_row[3] if source_row is not None else "NONE"

        event_type_match = re.search(r"SOFT|HARD|DEL|INS", identifier)
        event_type = event_type_match.group(0) if event_type_match else "UNDEFINED"
        if intersects_same_family_deletion(position, deletions):
            event_type += "_DEL"

        tsd_row = tsds.get(identifier)
        if tsd_row is None:
            tsd = "NONE"
            sniffles_row = source_tables[0].get(identifier)
            te_size = sniffles_row[6] if sniffles_row is not None else "NONE"
            new_position = position["start"]
        else:
            tsd = tsd_row[1]
            te_size = tsd_row[3].split(":", 1)[0] or "NONE"
            new_position = tsd_row[2]

        frequency_before, frequency_after = outsider_frequency_values(
            identifier, precise, raw_frequency
        )
        sv_size = maximum_size(sizes.get(identifier, ()))

        if re.search(r"SOFT|HARD", identifier):
            tremolo_id = "TE_ID_OUTSIDER.{}".format(identifier)
        else:
            call_type = outsider_variant_type(
                identifier, sniffles_calls, direct_calls
            )
            support_row = raw_frequency.get(identifier)
            read_support = support_row[3] if support_row is not None else ""
            tremolo_id = "TE_ID_OUTSIDER.{}.{}.{}".format(
                legacy_chain_id(
                    "{}:{}".format(position["chromosome"], position["start"])
                ),
                call_type,
                legacy_chain_id(read_support),
            )

        rendered.append(
            (
                position["chromosome"],
                position["start"],
                position["end"],
                "{}|{}".format(position["family"], identifier),
                position["strand"],
                tsd,
                pident,
                size_percent,
                te_size,
                new_position,
                frequency_before,
                frequency_after,
                sv_size,
                tremolo_id,
                event_type,
            )
        )
    return rendered


def index_insider_depth(rows):
    indexed = OrderedDict()
    for row in rows:
        require_width(row, 8, "INSIDER frequency")
        indexed.setdefault(row[7], row)
    return indexed


def index_variants(rows):
    indexed = OrderedDict()
    for row in rows:
        require_width(row, 10, "INSIDER structural variants")
        indexed.setdefault(row[3], row)
    return indexed


def insider_sv_size(identifier, variants):
    row = variants.get(identifier)
    if row is None:
        return "NONE"
    variant_type = row[6]
    if variant_type in CONTRACTION_TYPES:
        return row[4] or "NONE"
    if variant_type in EXPANSION_TYPES:
        return row[8] or "NONE"
    return "NONE"


def render_insider_source(
    position_rows, combine, tsds, depth, variants, source_name, minimum_columns
):
    rendered = []
    for raw_position in position_rows:
        position = parse_position(raw_position, source_name, minimum_columns)
        identifier = position["identifier"]
        call = combine.get(identifier)
        if call is None:
            raise ValueError(
                "{}: no classified TE call for {}".format(source_name, identifier)
            )
        require_width(call, 5, source_name)
        qseqid_fields = call[1].split(":")
        if len(qseqid_fields) < 3:
            raise ValueError("{}: malformed qseqid: {}".format(source_name, call[1]))
        event_type = qseqid_fields[2]

        tsd_row = tsds.get(identifier)
        if tsd_row is None:
            tsd = "NONE"
            te_size = call[4] or "NONE"
            new_position = position["start"]
        else:
            tsd = tsd_row[1]
            te_size = tsd_row[3] or "NONE"
            new_position = tsd_row[2]

        name = "{}|{}".format(position["family"], identifier)
        depth_row = depth.get(name)
        frequency = depth_row[5] if depth_row is not None else "NONE"
        tremolo_id = "TE_ID_INSIDER.{}.{}".format(
            legacy_chain_id(
                "{}:{}".format(position["chromosome"], position["start"])
            ),
            event_type,
        )
        rendered.append(
            (
                position["chromosome"],
                position["start"],
                position["end"],
                name,
                position["strand"],
                tsd,
                call[2],
                call[3],
                te_size,
                new_position,
                frequency,
                "INSIDER",
                insider_sv_size(identifier, variants),
                tremolo_id,
                event_type,
            )
        )
    return rendered


def render_insider_rows(args):
    insertion_combine = first_by_event(
        read_rows(args.insider_insertion_combine, skip_header=True),
        insider_event_id,
        "INSIDER insertion calls",
    )
    deletion_combine = first_by_event(
        read_rows(args.insider_deletion_combine, skip_header=True),
        insider_event_id,
        "INSIDER deletion calls",
    )
    tsds = tsd_by_event(read_rows(args.insider_tsd), "INSIDER TSD")
    depth = index_insider_depth(
        read_rows(args.insider_frequency, skip_header=True)
    )
    variants = index_variants(read_rows(args.insider_variants))

    rendered = render_insider_source(
        read_rows(args.insider_positions),
        insertion_combine,
        tsds,
        depth,
        variants,
        "INSIDER insertion positions",
        6,
    )
    rendered.extend(
        render_insider_source(
            read_rows(args.insider_deletion_bed),
            deletion_combine,
            tsds,
            depth,
            variants,
            "INSIDER deletion positions",
            5,
        )
    )
    return rendered


def render(args):
    rows = render_outsider_rows(args)
    rows.extend(render_insider_rows(args))
    return rows


def write_output(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    try:
        with temporary.open("w") as handle:
            handle.write("\t".join(HEADER) + "\n")
            for row in rows:
                if len(row) != len(HEADER):
                    raise ValueError("TE_INFOS row does not have 15 columns")
                handle.write("\t".join(str(value) for value in row) + "\n")
        os.replace(str(temporary), str(path))
    finally:
        if temporary.exists():
            temporary.unlink()


def build(args):
    rows = render(args)
    write_output(args.output, rows)
    return len(rows)


def build_parser():
    parser = argparse.ArgumentParser()
    for option in (
        "outsider-positions",
        "outsider-sniffles-combine",
        "outsider-direct-combine",
        "outsider-soft-calls",
        "outsider-hard-calls",
        "outsider-tsd",
        "outsider-frequency-precise",
        "outsider-frequency",
        "outsider-sv-sizes",
        "outsider-sniffles-calls",
        "outsider-direct-calls",
        "insider-deletion-bed",
        "insider-positions",
        "insider-insertion-combine",
        "insider-deletion-combine",
        "insider-tsd",
        "insider-frequency",
        "insider-variants",
        "output",
    ):
        parser.add_argument("--" + option, required=True, type=Path)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    count = build(args)
    print("Wrote {} TE_INFOS records.".format(count))


if __name__ == "__main__":
    main()
