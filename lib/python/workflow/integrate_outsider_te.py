#!/usr/bin/env python3
"""Build genomes containing the non-clipped OUTSIDER insertion calls.

The historical rule implemented this transformation through nested ``grep``
loops and inserted the event identifier itself before every canonical TE.  The
identifier made ``PSEUDO_GENOME_TE_DB_ID.fasta`` an invalid DNA FASTA.  This
implementation keeps the public filenames and coordinate semantics, but stores
event identity in BED/audit records and inserts nucleotide sequence only.
"""

from __future__ import print_function

import argparse
import csv
import sys
from collections import OrderedDict
from pathlib import Path


DNA_BASES = frozenset("ACGTURYKMSWBDHVN-")
AUDIT_HEADER = (
    "event_id",
    "te_family",
    "source",
    "qseqid",
    "chromosome",
    "breakpoint",
    "strand",
    "canonical_length",
    "observed_length",
    "status",
    "reason",
)


def read_fasta(path):
    """Read a possibly wrapped FASTA into an ordered mapping."""

    records = OrderedDict()
    name = None
    sequence = []
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    if name in records:
                        raise ValueError("{}: duplicate FASTA record {}".format(path, name))
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
        if name in records:
            raise ValueError("{}: duplicate FASTA record {}".format(path, name))
        records[name] = "".join(sequence)
    return records


def reverse_complement(sequence):
    translation = str.maketrans(
        "ACGTURYKMSWBDHVNacgturykmswbdhvn-",
        "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn-",
    )
    return sequence.translate(translation)[::-1]


def event_id(qseqid):
    fields = qseqid.split(":")
    return fields[4] if len(fields) > 4 else ""


def sequence_record_id(qseqid):
    """Return the FASTA identifier after removing the BLAST-only strand tag."""

    return qseqid[:-2] if qseqid.endswith((":+", ":-")) else qseqid


def read_call_table(path, source):
    """Index classified BLAST calls by their stable event identifier."""

    calls = OrderedDict()
    with Path(path).open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for line_number, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) < 2:
                raise ValueError(
                    "{}:{}: expected at least two columns".format(path, line_number)
                )
            identifier = event_id(row[1])
            if not identifier:
                raise ValueError(
                    "{}:{}: cannot extract event ID from {}".format(
                        path, line_number, row[1]
                    )
                )
            calls.setdefault(
                identifier,
                {"family": row[0], "qseqid": row[1], "source": source},
            )
    return calls


def read_merged_calls(path):
    """Return non-clipped merged calls in their declared order."""

    records = []
    with Path(path).open() as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(
                    "{}:{}: expected at least four BED columns".format(
                        path, line_number
                    )
                )
            name = fields[3]
            if "HARD" in name or "SOFT" in name:
                continue
            if "|" not in name:
                raise ValueError(
                    "{}:{}: malformed TE|event name {}".format(path, line_number, name)
                )
            family, identifier = name.rsplit("|", 1)
            records.append(
                {
                    "family": family,
                    "event_id": identifier,
                    "merged_chromosome": fields[0],
                    "merged_start": fields[1],
                }
            )
    return records


def validate_dna(sequence, label):
    unsupported = sorted(set(sequence.upper()).difference(DNA_BASES))
    if unsupported:
        raise ValueError(
            "{} contains unsupported nucleotide symbols: {}".format(
                label, ",".join(unsupported)
            )
        )


def resolve_insertions(
    merged_calls,
    sniffles_calls,
    direct_calls,
    sniffles_sequences,
    direct_sequences,
    te_sequences,
    genome_sequences,
):
    """Resolve each merged call to canonical and observed sequence records."""

    insertions = []
    audit = []
    errors = []
    seen_events = set()
    for merged in merged_calls:
        identifier = merged["event_id"]
        row = {
            "event_id": identifier,
            "te_family": merged["family"],
            "source": "",
            "qseqid": "",
            "chromosome": merged["merged_chromosome"],
            "breakpoint": merged["merged_start"],
            "strand": "",
            "canonical_length": "",
            "observed_length": "",
            "status": "rejected",
            "reason": "",
        }
        if identifier in seen_events:
            row["reason"] = "duplicate_merged_event"
            audit.append(row)
            errors.append("duplicate merged event {}".format(identifier))
            continue
        seen_events.add(identifier)

        if identifier.startswith("sniffles."):
            call = sniffles_calls.get(identifier)
            source_sequences = sniffles_sequences
        elif identifier.startswith("TrEMOLO."):
            call = direct_calls.get(identifier)
            source_sequences = direct_sequences
        else:
            call = sniffles_calls.get(identifier) or direct_calls.get(identifier)
            source_sequences = (
                sniffles_sequences if identifier in sniffles_calls else direct_sequences
            )
        if call is None:
            row["reason"] = "missing_classified_call"
            audit.append(row)
            errors.append("{} has no classified call".format(identifier))
            continue

        qseqid = call["qseqid"]
        fields = qseqid.split(":")
        row["source"] = call["source"]
        row["qseqid"] = qseqid
        if len(fields) < 6:
            row["reason"] = "malformed_qseqid"
            audit.append(row)
            errors.append("{} has malformed qseqid {}".format(identifier, qseqid))
            continue
        chromosome = fields[0]
        strand = fields[-1] if fields[-1] in ("+", "-") else "+"
        row["chromosome"] = chromosome
        row["breakpoint"] = fields[2]
        row["strand"] = strand
        try:
            breakpoint = int(fields[2])
        except ValueError:
            row["reason"] = "invalid_breakpoint"
            audit.append(row)
            errors.append("{} has invalid breakpoint {}".format(identifier, fields[2]))
            continue
        if chromosome not in genome_sequences:
            row["reason"] = "missing_genome_chromosome"
            audit.append(row)
            errors.append("{} refers to missing chromosome {}".format(identifier, chromosome))
            continue
        if breakpoint < 0 or breakpoint > len(genome_sequences[chromosome]):
            row["reason"] = "breakpoint_out_of_bounds"
            audit.append(row)
            errors.append("{} breakpoint is outside {}".format(identifier, chromosome))
            continue
        if call["family"] != merged["family"]:
            row["reason"] = "family_mismatch"
            audit.append(row)
            errors.append(
                "{} family differs between merged BED ({}) and call table ({})".format(
                    identifier, merged["family"], call["family"]
                )
            )
            continue
        canonical = te_sequences.get(merged["family"])
        observed = source_sequences.get(sequence_record_id(qseqid))
        if canonical is None:
            row["reason"] = "missing_canonical_sequence"
            audit.append(row)
            errors.append("{} has no canonical TE sequence".format(identifier))
            continue
        if observed is None:
            row["reason"] = "missing_observed_sequence"
            audit.append(row)
            errors.append("{} has no observed insertion sequence".format(identifier))
            continue
        canonical = canonical.upper()
        observed = observed.upper()
        if strand == "-":
            canonical = reverse_complement(canonical)
        try:
            validate_dna(canonical, "canonical sequence for {}".format(identifier))
            validate_dna(observed, "observed sequence for {}".format(identifier))
        except ValueError as error:
            row["reason"] = "invalid_sequence"
            audit.append(row)
            errors.append(str(error))
            continue

        row["canonical_length"] = str(len(canonical))
        row["observed_length"] = str(len(observed))
        row["status"] = "integrated"
        row["reason"] = "integrated"
        audit.append(row)
        insertions.append(
            {
                "event_id": identifier,
                "family": merged["family"],
                "chromosome": chromosome,
                "breakpoint": breakpoint,
                "strand": strand,
                "canonical": canonical,
                "observed": observed,
            }
        )
    return insertions, audit, errors


def integrate(genome_sequences, insertions, sequence_key):
    """Insert sequences and return the new genome plus shifted BED rows."""

    by_chromosome = OrderedDict((name, []) for name in genome_sequences)
    for order, insertion in enumerate(insertions):
        record = dict(insertion)
        record["order"] = order
        by_chromosome[insertion["chromosome"]].append(record)

    output = OrderedDict()
    positions = []
    for chromosome, genome in genome_sequences.items():
        records = sorted(
            by_chromosome[chromosome],
            key=lambda record: (record["breakpoint"], record["order"]),
        )
        chunks = []
        cursor = 0
        shift = 0
        for record in records:
            breakpoint = record["breakpoint"]
            sequence = record[sequence_key]
            chunks.append(genome[cursor:breakpoint])
            chunks.append(sequence)
            start = breakpoint + shift
            end = start + len(sequence)
            positions.append(
                (
                    chromosome,
                    start,
                    end,
                    "{}:{}".format(record["family"], record["event_id"]),
                )
            )
            cursor = breakpoint
            shift += len(sequence)
        chunks.append(genome[cursor:])
        output[chromosome] = "".join(chunks)
    return output, positions


def write_fasta(records, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as handle:
        for name, sequence in records.items():
            handle.write(">{}\n{}\n".format(name, sequence))


def write_bed(rows, path, public_names=False):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w") as handle:
        for chromosome, start, end, name in rows:
            if public_names:
                name = name.replace(":", "|")
            handle.write("{}\t{}\t{}\t{}\n".format(chromosome, start, end, name))


def write_audit(rows, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def build(args):
    genome = read_fasta(args.genome)
    te_sequences = read_fasta(args.te_database)
    sniffles_sequences = read_fasta(args.sniffles_fasta)
    direct_sequences = read_fasta(args.direct_fasta)
    merged = read_merged_calls(args.merged_bed)
    sniffles_calls = read_call_table(args.sniffles_calls, "sniffles")
    direct_calls = read_call_table(args.direct_calls, "direct")

    insertions, audit, errors = resolve_insertions(
        merged,
        sniffles_calls,
        direct_calls,
        sniffles_sequences,
        direct_sequences,
        te_sequences,
        genome,
    )
    write_audit(audit, args.audit)
    if errors:
        for error in errors:
            print("Rejected OUTSIDER integration: {}".format(error), file=sys.stderr)

    canonical_genome, canonical_positions = integrate(genome, insertions, "canonical")
    observed_genome, observed_positions = integrate(genome, insertions, "observed")
    write_fasta(canonical_genome, args.canonical_genome)
    write_fasta(observed_genome, args.observed_genome)
    write_bed(canonical_positions, args.canonical_bed)
    write_bed(observed_positions, args.observed_bed)
    write_bed(canonical_positions, args.canonical_public_bed, public_names=True)
    write_bed(observed_positions, args.observed_public_bed, public_names=True)
    return len(insertions)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--te-database", type=Path, required=True)
    parser.add_argument("--merged-bed", type=Path, required=True)
    parser.add_argument("--sniffles-calls", type=Path, required=True)
    parser.add_argument("--direct-calls", type=Path, required=True)
    parser.add_argument("--sniffles-fasta", type=Path, required=True)
    parser.add_argument("--direct-fasta", type=Path, required=True)
    parser.add_argument("--canonical-genome", type=Path, required=True)
    parser.add_argument("--observed-genome", type=Path, required=True)
    parser.add_argument("--canonical-bed", type=Path, required=True)
    parser.add_argument("--observed-bed", type=Path, required=True)
    parser.add_argument("--canonical-public-bed", type=Path, required=True)
    parser.add_argument("--observed-public-bed", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    count = build(args)
    print("Integrated {} OUTSIDER calls into both genomes.".format(count))


if __name__ == "__main__":
    main()
