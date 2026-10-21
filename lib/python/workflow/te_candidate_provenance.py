"""Describe raw sequence evidence without inferring biological alleles."""
import csv
import hashlib
import math
from collections import defaultdict
from pathlib import Path
from genomic_alignment import ALIGNMENT_COLUMNS, alignment_metadata

SUMMARY_COLUMNS = ("query_count", "sequence_count", "read_count", "vcf_query_count", "provenance_status", "relation_types")
EVIDENCE_COLUMNS = ("candidate_group_id", "event_id", "evidence_channel", "query_id", "read_id",
                    "sequence_origin", "insertion_sha256", "insertion_length", "candidate_te",
                    "query_start", "query_end", "te_start", "te_end", "strand", "identity",
                    "aligned_bp", "evalue", "bitscore", "relation_types") + ALIGNMENT_COLUMNS
PROVENANCE_COLUMNS = ("query_id", "event_id", "read_id", "sequence_origin")


def read_provenance(path, events=None):
    result = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        header = tuple(reader.fieldnames or ())
        if header not in (PROVENANCE_COLUMNS, PROVENANCE_COLUMNS + ALIGNMENT_COLUMNS):
            raise ValueError("Invalid sequence provenance header: " + str(path))
        for row in reader:
            if events is not None and row["event_id"] not in events:
                continue
            if row["query_id"] in result:
                raise ValueError("Duplicate provenance query: " + row["query_id"])
            if row["sequence_origin"] not in {"read_sequence", "vcf_reported"}:
                raise ValueError("Unknown sequence origin")
            if row["sequence_origin"] == "vcf_reported" and row["read_id"] != ".":
                raise ValueError("A VCF allele cannot be assigned to a scanned read")
            if row["sequence_origin"] == "read_sequence" and row["read_id"] in ("", "."):
                raise ValueError("A read sequence must have its source read identifier")
            metadata = alignment_metadata(row.get("genomic_alignment_flag"), row["sequence_origin"])
            if "genomic_alignment_status" in row and row["genomic_alignment_status"] != metadata["genomic_alignment_status"]:
                raise ValueError("Genomic alignment status disagrees with SAM flag")
            row.update(metadata)
            result[row["query_id"]] = row
    return result


def read_sequences(path, events=None):
    result, name, pieces = {}, None, []
    def save():
        if name is not None:
            if events is not None and name.split(":")[4] not in events:
                return
            if name in result:
                raise ValueError("Duplicate insertion query: " + name)
            sequence = "".join(pieces).upper()
            result[name] = (hashlib.sha256(sequence.encode("ascii")).hexdigest(), len(sequence))
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                name, pieces = line[1:], []
            else:
                if name is None:
                    raise ValueError("Sequence before FASTA header")
                pieces.append(line)
    save()
    return result


def load_evidence(alignments, provenances, fastas, events=None):
    provenance = {c: read_provenance(p, events) for c, p in provenances if str(p) != "/dev/null"}
    sequences = {c: read_sequences(p, events) for c, p in fastas if str(p) != "/dev/null"}
    by_event = defaultdict(list)
    unknown_alignment = alignment_metadata()
    for channel, path in alignments:
        if str(path) == "/dev/null":
            continue
        with Path(path).open() as handle:
            for number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) != 12:
                    raise ValueError(f"{path}:{number}: expected 12 BLAST columns")
                query, family = fields[:2]
                parts = query.split(":")
                if len(parts) < 5:
                    raise ValueError("Missing event in query: " + query)
                if events is not None and parts[4] not in events:
                    continue
                q1, q2, s1, s2 = map(int, fields[6:10])
                identity, evalue, score = float(fields[2]), float(fields[10]), float(fields[11])
                if min(q1, q2, s1, s2) < 1 or not all(map(math.isfinite, (identity, evalue, score))) or not 0 <= identity <= 100 or evalue < 0 or int(fields[3]) < 1:
                    raise ValueError("Invalid BLAST alignment: " + query)
                origin = provenance.get(channel, {}).get(query)
                if origin is not None and origin["event_id"] != parts[4]:
                    raise ValueError("Inconsistent event provenance: " + query)
                digest, length = sequences.get(channel, {}).get(query, (".", "."))
                if channel in provenance and origin is None:
                    raise ValueError("Missing declared provenance: " + query)
                if channel in sequences and digest == ".":
                    raise ValueError("Missing declared insertion sequence: " + query)
                by_event[parts[4]].append(dict(event_id=parts[4], evidence_channel=channel, query_id=query,
                    read_id=origin["read_id"] if origin else ".", sequence_origin=origin["sequence_origin"] if origin else "unknown",
                    **{key: (origin or unknown_alignment)[key] for key in ALIGNMENT_COLUMNS},
                    insertion_sha256=digest, insertion_length=length, candidate_te=family,
                    query_start=min(q1,q2)-1, query_end=max(q1,q2), te_start=min(s1,s2)-1, te_end=max(s1,s2),
                    strand="+" if (q2-q1)*(s2-s1) >= 0 else "-", identity=identity,
                    aligned_bp=int(fields[3]), evalue=evalue, bitscore=score))
    return by_event


def summarize(rows, families):
    rows = [row for row in rows if row["candidate_te"] in families]
    by_query = defaultdict(list)
    for row in rows:
        by_query[(row["evidence_channel"], row["query_id"])].append(row)
    relations = set()
    for values in by_query.values():
        labels = sorted({row["candidate_te"] for row in values})
        for index, family in enumerate(labels):
            for other in labels[index+1:]:
                overlaps = [max(0, min(a["query_end"], b["query_end"])-max(a["query_start"], b["query_start"]))
                            for a in values if a["candidate_te"] == family for b in values if b["candidate_te"] == other]
                relations.add("same_query_overlapping_families" if any(overlaps) else "same_query_disjoint_families")
    digests = {row["insertion_sha256"] for row in rows if row["insertion_sha256"] != "."}
    reads = {row["read_id"] for row in rows if row["read_id"] != "."}
    if len(digests) > 1:
        relations.add("different_insertion_sequences")
    if len(reads) > 1:
        relations.add("different_reads")
    if len(digests) == 1 and len(by_query) > 1:
        relations.add("identical_insertion_multiple_queries")
    complete = bool(rows) and all(row["sequence_origin"] != "unknown" and row["insertion_sha256"] != "." for row in rows)
    annotation = "same_query_overlapping_families" in relations
    components = "same_query_disjoint_families" in relations
    different = len(digests) > 1
    if sum((annotation, components, different)) > 1:
        kind = "mixed_sequence_evidence"
    elif annotation:
        kind = "same_sequence_annotation_ambiguity"
    elif components:
        kind = "same_sequence_multiple_regions"
    elif different:
        kind = "different_sequence_evidence"
    else:
        kind = "unresolved_family_candidates"
    summary = dict(query_count=len(by_query), sequence_count=len(digests), read_count=len(reads),
                   vcf_query_count=sum(any(row["sequence_origin"] == "vcf_reported" for row in values) for values in by_query.values()),
                   provenance_status="complete" if complete else "partial" if rows else "unavailable",
                   relation_types=";".join(sorted(relations)) or ".", ambiguity_type=kind)
    return summary, rows
