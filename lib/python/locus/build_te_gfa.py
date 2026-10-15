#!/usr/bin/env python3
"""Build a reference-backed GFA with ranked, synthetic TE insertion walks."""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

from build_cohort_model import anchor, cluster_calls, event_class
from build_locus_model import InputRun, file_sha256, read_calls, stable_id, write_tsv
from gfa_summary import summarize, write_summary

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow"))
from te_candidate_provenance import load_evidence
from restore_te_names import read_names


def read_fasta(path, full_headers=False, allow_empty=False):
    sequences, name, pieces = {}, None, []
    def save():
        if name is not None:
            sequence = "".join(pieces).upper()
            if name in sequences or not sequence or re.fullmatch("[ACGTRYSWKMBDHVN]+", sequence) is None:
                raise ValueError(f"{path}: duplicate, empty or invalid DNA sequence: {name}")
            sequences[name] = sequence
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                save()
                header = line[1:]
                if not header:
                    raise ValueError(f"{path}: empty FASTA header")
                name, pieces = (header if full_headers else header.split()[0]), []
            elif name is None:
                raise ValueError(f"{path}: sequence before FASTA header")
            else:
                pieces.append(line)
    save()
    if not sequences and not allow_empty:
        raise ValueError(f"{path}: empty FASTA")
    return sequences


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def best_chain(rows):
    """Maximum bitscore collinear chain, non-overlapping on query AND TE."""
    best = []
    for strand in ("+", "-"):
        ordered = sorted((r for r in rows if r["strand"] == strand), key=lambda r: (
            r["query_end"], r["query_start"], r["te_start"], r["te_end"], -r["bitscore"]))
        chains = []
        for i, row in enumerate(ordered):
            chain = [row]
            for j in range(i):
                previous = ordered[j]
                subject_ok = (previous["te_end"] <= row["te_start"] if strand == "+"
                              else previous["te_start"] >= row["te_end"])
                if previous["query_end"] <= row["query_start"] and subject_ok:
                    candidate = chains[j] + [row]
                    if sum(r["bitscore"] for r in candidate) > sum(r["bitscore"] for r in chain):
                        chain = candidate
            chains.append(chain)
            if sum(r["bitscore"] for r in chain) > sum(r["bitscore"] for r in best):
                best = chain
    return best


def rank(row):
    return (-row["bitscore"], -row["coverage"], -row["identity"],
            -row["support_reads"], row["sample_id"], row["channel"], row["query_id"], row["family"])


def run_input(value):
    sample, separator, directory = value.partition("=")
    if not separator or not sample or not directory:
        raise ValueError("Run must use SAMPLE=/path/to/analysis")
    root = Path(directory).resolve()
    spec = dict(sample_id=sample, te_infos=str(root / "TE_INFOS.bed"),
                reference=str(root / "INPUT/genome.fasta"), channels=[])
    mapping = root / "INPUT/te_headers.tsv"
    if mapping.is_file():
        spec["te_headers"] = str(mapping)
    for name, fasta, provenance, blast in (
        ("outsider_direct_alignment", "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.insertions.fasta",
         "OUTSIDER/TrEMOLO_SV_TE/INS/SEQUENCE_PROVENANCE.tsv", "OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln"),
        ("outsider_sniffles", "OUTSIDER/VARIANT_CALLING/SEQUENCE_INDEL.insertions.fasta",
         "OUTSIDER/VARIANT_CALLING/SEQUENCE_PROVENANCE.tsv", "OUTSIDER/TE_DETECTION/BLAST_SEQUENCE_INDEL_vs_DBTE.bln")):
        paths = [root / p for p in (fasta, provenance, blast)]
        if any(p.exists() for p in paths):
            for path in paths:
                if not path.is_file():
                    raise ValueError(f"Incomplete channel {name}: missing {path}; rerun sequence extraction")
            spec["channels"].append(dict(name=name, insertion_fasta=str(paths[0]),
                                         provenance=str(paths[1]), alignment=str(paths[2])))
    return spec


def collect(inputs, database, min_identity, min_coverage, min_aligned_bp, chromosomes=None):
    calls, candidates, exclusions, used_paths = [], [], [], []
    seen_samples = set()
    for spec in inputs:
        sample = spec["sample_id"]
        if not sample or sample in seen_samples:
            raise ValueError("Empty or duplicate sample_id")
        seen_samples.add(sample)
        path = Path(spec["te_infos"])
        used_paths.append(path)
        accepted = {}
        for call in read_calls(InputRun(sample, path)):
            event = call["te_call"].rpartition("|")[2]
            if not event or event in accepted:
                raise ValueError(f"{path}: missing or duplicate event ID: {event}")
            accepted[event] = call
            if chromosomes is not None and call["chrom"] not in chromosomes:
                exclusions.append(dict(sample_id=sample, event_id=event, reason="excluded_by_chrom_filter"))
                continue
            if event_class(call["event_type"]) == "insertion":
                calls.append(call)
            else:
                exclusions.append(dict(sample_id=sample, event_id=event, reason="unsupported_event_type"))
        names = read_names(spec["te_headers"]) if spec.get("te_headers") else {}
        if spec.get("te_headers"):
            used_paths.append(Path(spec["te_headers"]))
        seen_channels = set()
        for channel in spec.get("channels", []):
            label = channel["name"]
            if label in seen_channels:
                raise ValueError("Duplicate evidence channel: " + label)
            seen_channels.add(label)
            paths = {key: Path(channel[key]) for key in ("alignment", "provenance", "insertion_fasta")}
            used_paths.extend(paths.values())
            sequences = read_fasta(paths["insertion_fasta"], allow_empty=True)
            events = {e for e, c in accepted.items() if event_class(c["event_type"]) == "insertion"
                      and (chromosomes is None or c["chrom"] in chromosomes)}
            evidence = load_evidence([(label, paths["alignment"])], [(label, paths["provenance"])],
                                     [(label, paths["insertion_fasta"])], events)
            grouped = defaultdict(list)
            for rows in evidence.values():
                for row in rows:
                    row["candidate_te"] = names.get(row["candidate_te"], row["candidate_te"])
                    if row["candidate_te"] not in database:
                        raise ValueError("TE absent from database: " + row["candidate_te"])
                    if row["te_end"] > len(database[row["candidate_te"]]) or row["bitscore"] < 0:
                        raise ValueError("Invalid TE coordinates or negative bitscore")
                    if row["identity"] >= min_identity:
                        grouped[(row["event_id"], row["query_id"], row["candidate_te"])].append(row)
            for (event, query, family), rows in sorted(grouped.items()):
                chain = best_chain(rows)
                if not chain:
                    continue
                weight = sum(r["aligned_bp"] for r in chain)
                coverage = 100 * sum(r["te_end"] - r["te_start"] for r in chain) / len(database[family])
                first = chain[0]
                sequence = sequences[query]
                candidates.append(dict(sample_id=sample, event_id=event, channel=label, query_id=query,
                    read_id=first["read_id"], sequence_origin=first["sequence_origin"],
                    insertion_sha256=first["insertion_sha256"], family=family, strand=first["strand"],
                    bitscore=sum(r["bitscore"] for r in chain), coverage=coverage,
                    identity=sum(r["identity"] * r["aligned_bp"] for r in chain) / weight,
                    support_reads=0, sequence=sequence, locus_id=".", node_id=".", level=".",
                    status="eligible" if coverage >= min_coverage and weight >= min_aligned_bp else "below_threshold"))
    return calls, candidates, exclusions, used_paths


def build(reference, te_database, inputs, output, reference_name="reference", locus_window=0,
          min_identity=80, min_coverage=80, min_aligned_bp=100, selected_tes=None, chrom_regex=None,
          only_called_chromosomes=True):
    try:
        chrom_pattern = re.compile(chrom_regex) if chrom_regex is not None else None
    except re.error as error:
        raise ValueError(f"Invalid chromosome regex: {error}") from error
    if selected_tes is not None:
        selected_tes = set(selected_tes)
        if not selected_tes:
            raise ValueError("TE selection must contain at least one name")
    if locus_window < 0 or min_aligned_bp < 1 or not all(
        math.isfinite(v) and 0 <= v <= 100 for v in (min_identity, min_coverage)):
        raise ValueError("Invalid window or alignment thresholds")
    if re.fullmatch(r"[A-Za-z0-9_.-]+", reference_name) is None or reference_name.startswith("level_"):
        raise ValueError("Invalid reference sample name (level_ is reserved)")
    genome, database = read_fasta(reference), read_fasta(te_database, full_headers=True)
    for spec in inputs:
        if spec.get("reference") and read_fasta(spec["reference"]) != genome:
            raise ValueError("Input run uses a different reference: " + spec["sample_id"])
    if chrom_pattern is not None:
        genome = {chrom: sequence for chrom, sequence in genome.items() if chrom_pattern.search(chrom)}
        if not genome:
            raise ValueError("Chromosome regex matches no reference sequences")
    if only_called_chromosomes:
        called_chromosomes = {
            call["chrom"] for spec in inputs
            for call in read_calls(InputRun(spec["sample_id"], Path(spec["te_infos"])))
        }
        genome = {chrom: sequence for chrom, sequence in genome.items() if chrom in called_chromosomes}
        if not genome:
            raise ValueError("No reference chromosomes have calls in TE_INFOS.bed after chromosome filtering")
    for chrom in genome:
        if re.fullmatch(r"[A-Za-z0-9_.:-]+", chrom) is None:
            raise ValueError("Invalid chromosome name for GFA: " + chrom)
    calls, candidates, exclusions, paths = collect(inputs, database, min_identity, min_coverage, min_aligned_bp,
                                                  set(genome) if chrom_pattern is not None or only_called_chromosomes else None)
    paths.extend(Path(spec["reference"]) for spec in inputs if spec.get("reference"))
    destinations = {Path(output).resolve() / name for name in (
        "te-graph.gfa", "sequence-alternatives.tsv", "loci.tsv", "excluded-events.tsv", "walks.tsv", "manifest.json",
        "summary.json", "summary.html", "call-representation.tsv")}
    if any(Path(p).resolve() in destinations for p in [reference, te_database] + paths):
        raise ValueError("Output would overwrite an input file")
    by_event = defaultdict(list)
    for row in candidates:
        by_event[(row["sample_id"], row["event_id"])].append(row)
    by_chrom = defaultdict(list)
    locus_rows = []
    reference_hash = file_sha256(reference)
    for members in cluster_calls(calls, locus_window):
        chrom = members[0]["chrom"]
        position = min(anchor(c) for c in members)
        if chrom not in genome or not 0 <= position <= len(genome[chrom]):
            raise ValueError(f"Invalid reference anchor: {chrom}:{position}")
        locus = stable_id("L", reference_hash, chrom, position, max(anchor(c) for c in members), locus_window)
        rows = []
        for call in members:
            event = call["te_call"].rpartition("|")[2]
            values = by_event[(call["sample_id"], event)]
            if not values:
                exclusions.append(dict(sample_id=call["sample_id"], event_id=event, reason="no_usable_sequence_evidence"))
            for row in values:
                row["locus_id"] = locus
            rows.extend(values)
        # One sequence with several family labels cannot become several variants.
        families_by_digest = defaultdict(set)
        for row in rows:
            if row["status"] == "eligible":
                families_by_digest[row["insertion_sha256"]].add(row["family"])
        for row in rows:
            if row["status"] == "eligible" and len(families_by_digest[row["insertion_sha256"]]) > 1:
                row["status"] = "ambiguous_or_composite_annotation"
        locus_families = {call["te_call"].rpartition("|")[0] for call in members}
        locus_families.update(row["family"] for row in rows if row["status"] != "below_threshold")
        selected = selected_tes is None or bool(locus_families & selected_tes)
        for row in rows:
            row["locus_selected"] = "yes" if selected else "no"
            if not selected and row["status"] == "eligible":
                row["status"] = "excluded_by_te_filter"
        grouped = defaultdict(list)
        for row in rows:
            if row["status"] == "eligible":
                grouped[row["family"]].append(row)
        representatives = []
        for family, values in sorted(grouped.items()):
            for row in values:
                row["support_reads"] = len({(r["sample_id"], r["read_id"]) for r in values
                    if r["read_id"] != "." and r["insertion_sha256"] == row["insertion_sha256"]})
            winner = min(values, key=rank)
            winner["status"] = "representative"
            for row in values:
                if row is not winner:
                    row["status"] = "same_family_alternative"
            representatives.append(winner)
        representatives.sort(key=rank)
        for level, row in enumerate(representatives, 1):
            row["level"] = level
        if selected:
            by_chrom[chrom].append((position, representatives))
        locus_rows.append(dict(locus_id=locus, chrom=chrom, position=position,
            anchor_min=position, anchor_max=max(anchor(c) for c in members),
            status="exact_anchor" if len({anchor(c) for c in members}) == 1 else "provisional_anchor",
            alternatives=len(representatives), selected="yes" if selected else "no"))
    lines = [f"H\tVN:Z:1.1\tRS:Z:{reference_name}"]
    node = 0
    link_count = 0
    walk_rows = []
    for chrom, sequence in genome.items():
        segments, backbone, choices = [], [], []
        previous = 0
        for position, variants in sorted(by_chrom[chrom], key=lambda item: item[0]):
            if position > previous:
                node += 1
                segments.append((node, sequence[previous:position]))
                backbone.append((node, "+"))
                choices.append([(node, "+")])
            branch = []
            for variant in variants:
                node += 1
                variant["node_id"] = node
                dna = variant["sequence"]
                segments.append((node, dna if variant["strand"] == "+" else reverse_complement(dna)))
                branch.append((node, variant["strand"]))
            if branch:
                choices.append(branch)
            previous = position
        if previous < len(sequence):
            node += 1
            segments.append((node, sequence[previous:]))
            backbone.append((node, "+"))
            choices.append([(node, "+")])
        lines.extend(f"S\t{i}\t{dna}" for i, dna in segments)
        lengths = {i: len(dna) for i, dna in segments}
        max_level = max((len(v) for _, v in by_chrom[chrom]), default=0)
        walks = [(reference_name, backbone)]
        for level in range(max_level):
            walk = [alternatives[min(level, len(alternatives)-1)] for alternatives in choices]
            walks.append((f"level_{level+1}", walk))
        links = set()
        for sample, walk in walks:
            dna_length = sum(lengths[i] for i, _ in walk)
            traversal = "".join((">" if orient == "+" else "<") + str(i) for i, orient in walk)
            lines.append(f"W\t{sample}\t0\t{chrom}\t0\t{dna_length}\t{traversal}")
            links.update((a, ao, b, bo) for (a, ao), (b, bo) in zip(walk, walk[1:]))
            walk_rows.append(dict(sample=sample, chrom=chrom, length=dna_length,
                                  interpretation="reference" if sample == reference_name else "synthetic_unphased_level"))
        lines.extend(f"L\t{a}\t{ao}\t{b}\t{bo}\t0M" for a, ao, b, bo in sorted(links))
        link_count += len(links)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "te-graph.gfa").write_text("\n".join(lines) + "\n", encoding="ascii")
    columns = ("locus_id", "sample_id", "event_id", "channel", "query_id", "read_id", "sequence_origin",
               "insertion_sha256", "family", "strand", "bitscore", "coverage", "identity", "support_reads",
               "status", "node_id", "level", "locus_selected")
    write_tsv(output / "sequence-alternatives.tsv", columns,
              ({k: r[k] for k in columns} for r in sorted(candidates, key=lambda r: (r["locus_id"], rank(r)))))
    write_tsv(output / "loci.tsv", ("locus_id", "chrom", "position", "anchor_min", "anchor_max", "status", "alternatives", "selected"), locus_rows)
    write_tsv(output / "excluded-events.tsv", ("sample_id", "event_id", "reason"), exclusions)
    write_tsv(output / "walks.tsv", ("sample", "chrom", "length", "interpretation"), walk_rows)
    summary, call_rows = summarize(inputs, candidates, locus_rows, exclusions, walk_rows, node, link_count)
    write_summary(output, summary, call_rows)
    manifest = dict(schema_version="0.4.0", reference_name=reference_name,
        chrom_regex=chrom_regex, chromosomes=list(genome),
        only_called_chromosomes=only_called_chromosomes,
        selected_tes=sorted(selected_tes) if selected_tes is not None else None,
        construction="family_representatives_synthetic_levels", locus_window=locus_window,
        ranking=["collinear_nonoverlapping_bitscore_sum", "coverage", "identity", "family_read_support", "identifier"],
        thresholds=dict(identity=min_identity, coverage=min_coverage, aligned_bp=min_aligned_bp),
        inputs=[dict(path=str(Path(p).resolve()), sha256=file_sha256(p)) for p in [reference, te_database] + paths],
        counts=dict(nodes=node, loci=len(locus_rows), selected_loci=sum(r["selected"] == "yes" for r in locus_rows),
                    representatives=sum(r["status"] == "representative" for r in candidates),
                    input_calls=summary["calls"]["total"], represented_calls=summary["calls"]["represented"],
                    excluded_calls=summary["calls"]["excluded"]), representation_summary="summary.json")
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--te-database", required=True, type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--inputs", type=Path, help="JSON list of sample/channel input paths")
    source.add_argument("--run", action="append", help="SAMPLE=/path/to/analysis; repeat for a cohort")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--reference-name", default="reference")
    parser.add_argument("--chrom-regex", help="Additional Python regex searched on reference sequence names")
    chromosome_scope = parser.add_mutually_exclusive_group()
    chromosome_scope.add_argument("--all-reference-chromosomes", dest="only_called_chromosomes", action="store_false",
                                 help="Include reference chromosomes absent from TE_INFOS.bed; still applies --chrom-regex")
    chromosome_scope.add_argument("--only-called-chromosomes", dest="only_called_chromosomes", action="store_true",
                                 help=argparse.SUPPRESS)
    parser.set_defaults(only_called_chromosomes=True)
    parser.add_argument("--locus-window", type=int, default=0)
    parser.add_argument("--min-identity", type=float, default=80)
    parser.add_argument("--min-coverage", type=float, default=80)
    parser.add_argument("--min-aligned-bp", type=int, default=100)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--te-names", help="Comma-separated original TE names; retains entire matching loci")
    selection.add_argument("--te-names-file", type=Path, help="File with one original TE name per line")
    args = parser.parse_args()
    inputs = [run_input(value) for value in args.run] if args.run else json.loads(args.inputs.read_text())
    selected_tes = None
    if args.te_names is not None:
        selected_tes = {name.strip() for name in args.te_names.split(",") if name.strip()}
    elif args.te_names_file is not None:
        selected_tes = {name.strip() for name in args.te_names_file.read_text().splitlines() if name.strip()}
    manifest = build(args.reference, args.te_database, inputs, args.output,
                     args.reference_name, args.locus_window, args.min_identity, args.min_coverage,
                     args.min_aligned_bp, selected_tes, args.chrom_regex, args.only_called_chromosomes)
    print(json.dumps(manifest["counts"]))
    counts = manifest["counts"]
    print(f"BED calls: {counts['represented_calls']}/{counts['input_calls']} represented; "
          f"{counts['excluded_calls']} excluded. Family representatives: {counts['representatives']}.")
    print(f"Representation details: {args.output / 'summary.html'}")


if __name__ == "__main__":
    main()
