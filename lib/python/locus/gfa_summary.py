#!/usr/bin/env python3
"""Explain how reported TE calls contribute to a GFA, without conflating calls and nodes."""
import argparse
import csv
import html
import json
from collections import Counter, defaultdict
from pathlib import Path

from build_locus_model import InputRun, parse_sample, read_calls, write_tsv

CALL_COLUMNS = ("sample_id", "event_id", "tremolo_id", "chrom", "reported_te", "locus_id",
                "status", "representation", "represented_families", "node_ids", "reported_family_represented")


def summarize(inputs, candidates, loci, exclusions, walks, nodes, links):
    by_call = defaultdict(list)
    representatives = {}
    for row in candidates:
        by_call[(row["sample_id"], row["event_id"])].append(row)
        if row["status"] == "representative" and str(row["node_id"]) != ".":
            representatives[(row["locus_id"], row["family"])] = row
    excluded = {(r["sample_id"], r["event_id"]): r["reason"] for r in exclusions}
    call_rows = []
    seen_calls, seen_samples = set(), set()
    for spec in inputs or []:
        if spec["sample_id"] in seen_samples:
            raise ValueError("Duplicate source sample: " + spec["sample_id"])
        seen_samples.add(spec["sample_id"])
        for call in read_calls(InputRun(spec["sample_id"], Path(spec["te_infos"]))):
            family, _, event = call["te_call"].rpartition("|")
            key = (spec["sample_id"], event)
            if key in seen_calls:
                raise ValueError("Duplicate source event: " + str(key))
            seen_calls.add(key)
            rows = by_call[key]
            contributing = [r for r in rows if r["status"] in ("representative", "same_family_alternative")
                            and (r["locus_id"], r["family"]) in representatives]
            represented = {r["family"]: representatives[(r["locus_id"], r["family"])] for r in contributing}
            own = any(r["status"] == "representative" for r in contributing)
            statuses = {r["status"] for r in rows}
            if represented:
                status = "represented"
            elif key in excluded:
                status = excluded[key]
            elif "excluded_by_te_filter" in statuses or any(r.get("locus_selected") == "no" for r in rows):
                status = "excluded_by_te_filter"
            elif "ambiguous_or_composite_annotation" in statuses:
                status = "ambiguous_or_composite_annotation"
            elif "below_threshold" in statuses:
                status = "below_threshold"
            else:
                status = "no_usable_sequence_evidence"
            call_rows.append(dict(sample_id=spec["sample_id"], event_id=event, tremolo_id=call["tremolo_id"],
                chrom=call["chrom"], reported_te=family, locus_id=rows[0]["locus_id"] if rows else ".",
                status=status, representation="own_sequence" if own else "same_locus_family" if represented else ".",
                represented_families=";".join(sorted(represented)) or ".",
                node_ids=";".join(str(represented[f]["node_id"]) for f in sorted(represented)) or ".",
                reported_family_represented="yes" if family in represented else "no"))
    counts = Counter(r["status"] for r in call_rows)
    call_summary = dict(available=inputs is not None, total=len(call_rows) if inputs is not None else None,
        represented=counts["represented"] if inputs is not None else None,
        excluded=len(call_rows) - counts["represented"] if inputs is not None else None,
        by_status=dict(sorted(counts.items())),
        by_representation=dict(sorted(Counter(r["representation"] for r in call_rows if r["status"] == "represented").items())),
        alternative_family_only=(sum(r["status"] == "represented" and r["reported_family_represented"] == "no"
                                     for r in call_rows) if inputs is not None else None))
    summary = dict(schema_version="1.0.0", calls=call_summary,
        loci=dict(total=len(loci), selected=sum(r["selected"] == "yes" for r in loci),
                  represented=len({r["locus_id"] for r in representatives.values()})),
        graph=dict(nodes=nodes, links=links, walks=len(walks), chromosomes=sorted({r["chrom"] for r in walks})),
        alternatives=dict(representatives=len(representatives),
                          by_family=dict(sorted(Counter(r["family"] for r in representatives.values()).items())),
                          by_strand=dict(sorted(Counter(r["strand"] for r in representatives.values()).items())),
                          candidate_sequences_by_status=dict(sorted(Counter(r["status"] for r in candidates).items()))))
    return summary, sorted(call_rows, key=lambda r: (r["sample_id"], r["event_id"]))


def write_summary(output, summary, call_rows):
    output = Path(output)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    write_tsv(output / "call-representation.tsv", CALL_COLUMNS, call_rows)
    def table(values):
        return "<table><tbody>" + "".join("<tr><th>" + html.escape(str(k)) + "</th><td>" +
            html.escape("Unavailable" if v is None else str(v)) + "</td></tr>" for k, v in values.items()) + "</tbody></table>"
    calls = summary["calls"]
    note = ("Call-level counts are unavailable: source TE_INFOS.bed files were not supplied. "
            "Candidate-sequence counts must not be interpreted as input call counts.") if not calls["available"] else (
            "A call is represented when one of its eligible sequences contributes to a retained family at its locus. "
            "Several calls can share one representative node. Representation may use an alternative family rather than the BED primary label.")
    page = """<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>TrEMOLO GFA representation summary</title><style>
body{font:16px system-ui;max-width:1000px;margin:3rem auto;padding:0 1rem;color:#152a35}
table{border-collapse:collapse;margin:1rem 0;min-width:24rem}th,td{padding:.5rem;border:1px solid #dce5e7;text-align:left}
th{background:#edf6f5}a{color:#087f78}.note{background:#fff4d6;padding:1rem}
</style><h1>GFA representation summary</h1>"""
    page += '<p class="note">' + html.escape(note) + "</p><h2>Input BED calls</h2>"
    page += table({k: calls[k] for k in ("available", "total", "represented", "excluded", "alternative_family_only")})
    page += "<h3>One outcome per call</h3>" + table(calls["by_status"])
    page += "<h3>How represented calls contribute</h3>" + table(calls["by_representation"])
    page += "<h2>Loci</h2>" + table(summary["loci"]) + "<h2>Graph</h2>" + table(summary["graph"])
    page += "<h2>Retained family representatives</h2>" + table(summary["alternatives"]["by_family"])
    page += "<h2>Candidate sequence outcomes (not call counts)</h2>" + table(summary["alternatives"]["candidate_sequences_by_status"])
    page += '<p><a href="call-representation.tsv">Call details TSV</a> · <a href="summary.json">Summary JSON</a> · <a href="sequence-alternatives.tsv">Sequence evidence TSV</a></p></html>\n'
    (output / "summary.html").write_text(page)


def refresh(directory, inputs=None):
    def rows(name):
        with (directory / name).open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    counts = Counter()
    with (directory / "te-graph.gfa").open() as handle:
        for line in handle:
            counts[line.split("\t", 1)[0]] += 1
    summary, calls = summarize(inputs, rows("sequence-alternatives.tsv"), rows("loci.tsv"),
        rows("excluded-events.tsv"), rows("walks.tsv"), counts["S"], counts["L"])
    write_summary(directory, summary, calls)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-dir", type=Path, required=True)
    parser.add_argument("--sample", action="append", type=parse_sample,
                        help="SAMPLE=/path/TE_INFOS.bed, repeat for every original run")
    args = parser.parse_args()
    inputs = [dict(sample_id=r.sample_id, te_infos=str(r.te_infos)) for r in args.sample] if args.sample else None
    print(json.dumps(refresh(args.graph_dir, inputs)["calls"]))


if __name__ == "__main__":
    main()
