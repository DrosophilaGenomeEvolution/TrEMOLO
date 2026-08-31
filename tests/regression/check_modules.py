#!/usr/bin/env python3
"""Validate the refactored module outputs on the bundled regression fixtures."""

import argparse
import csv
import json
from pathlib import Path


def require(condition, message):
    if not condition:
        raise SystemExit(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("population_report", type=Path)
    parser.add_argument("structure_report", type=Path)
    args = parser.parse_args()

    population = json.loads((args.population_report / "report-data.json").read_text())
    summary = population["summary"]
    require(summary["samples"] == 6, "population sample count differs")
    require(summary["loci"] == 214, "population locus count differs")
    require(summary["trajectories"] == 215, "population trajectory count differs")
    require(
        summary["observation_status"] == {"missing": 739, "observed": 551},
        "population missingness semantics differ",
    )
    require(
        all(
            row["frequency"] == "."
            for row in population["observations"]
            if row["observation_status"] == "missing"
        ),
        "a missing population observation was assigned a frequency",
    )
    population_html = (args.population_report / "report.html").read_text()
    require('id="trm-pop-scatter"' in population_html, "population chart missing")

    structure = json.loads((args.structure_report / "report-data.json").read_text())
    summary = structure["summary"]
    expected = {
        "queries": 800,
        "events": 220,
        "hsps": 276,
        "matches": 207,
        "retained_matches": 64,
        "components": 64,
        "final_events": 25,
    }
    for key, value in expected.items():
        require(summary[key] == value, f"structure {key} count differs")
    with (args.structure_report / "insertion-events.tsv").open(newline="") as handle:
        events = {row["event_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    event = events.get("TrEMOLO.INS.78")
    require(event is not None, "expected multi-TE event is absent")
    require(event["final_call"] == "yes", "multi-TE event is not linked to final call")
    require(event["reported_te"] == "R1A1-element", "reported TE differs")
    require(event["multi_component_queries"] == "2", "component support differs")
    require(event["multi_family_queries"] == "1", "multi-family support differs")
    require(event["classification"] == "multi_te_candidate", "event classification differs")
    structure_html = (args.structure_report / "report.html").read_text()
    require('id="trm-structure-svg"' in structure_html, "structure diagram missing")
    require('id="trm-event-body"' in structure_html, "event interpretation table missing")

    print(
        "Module regression: PASSED "
        "(214 loci, 215 trajectories; 220 events, 64 components)"
    )


if __name__ == "__main__":
    main()
