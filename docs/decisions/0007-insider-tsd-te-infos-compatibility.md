# 0007 — Preserve INSIDER TSD and `TE_INFOS.bed` compatibility

Status: accepted

## Context

The legacy `TSD_INSIDER` rule prepared all insertion and deletion flanks,
extracted genome sequence, and called TSDs inside one shell block. The report
rule then rebuilt `TE_INFOS.bed` with hundreds of repeated `grep` processes,
hidden inputs, a compiled `chain_to_id` helper, and side effects unrelated to
the scientific table. A report-rendering or network failure could therefore
prevent an otherwise complete call table from being represented in the DAG.

The reference combined run contains:

- 377 INSIDER TE candidates and 754 ten-base flank intervals;
- 377 formatted flank pairs and 210 INSIDER TSD calls;
- 57 OUTSIDER plus 377 INSIDER records in the final 15-column table;
- 219 calls with a TSD: 9 OUTSIDER and 210 INSIDER.

## Decision

INSIDER TSD preparation and calling are separate Snakemake rules. Genome
sequence is fetched through the declared normalized FASTA index, candidate
order is retained, and the longest exact common motif of at least four bases
is selected with the historical shift formulas. The four legacy text outputs
remain byte-identical.

`TE_INFOS.bed` is now built by one indexed Python aggregation rule whose every
scientific input is declared. It retains the historical source order
(OUTSIDER, INSIDER insertions, INSIDER deletions), values, and formatting. The
old C++ identifier calculation is reproduced with explicit float32 arithmetic,
so existing `ID_TrEMOLO` values do not change. Report rendering is no longer a
prerequisite for producing the public call table.

The rule consumes only the pipelines enabled in the configuration. INSIDER-
only and OUTSIDER-only runs therefore produce the corresponding subset without
requiring placeholder scientific results from the disabled branch.

Empty candidate sets produce empty TSD intermediates and a canonical
header-only `TE_INFOS.bed`; stale data is overwritten. Multiple calls at one
locus are deliberately retained as separate rows and are never collapsed by
the aggregation.

## Consequences

The work-test oracle is reproduced exactly: all four INSIDER TSD hashes and the
`TE_INFOS.bed` SHA-256
`eb545371e5c77ba705012ba952fb54593a477f97cfab465457e24a5d81a0e26a`.
Source-table lookups are indexed instead of repeatedly scanning every table
for every event.

This is a compatibility representation, not the population data model. Its
lossy position-derived identifiers, single-sample frequency fields, and TSD
motif heuristic remain unchanged. Cohort analyses and multi-allelic loci must
use the normalized locus model described in `docs/locus-data-model.md`, where
calls, alleles, and evidence have stable independent identifiers.
