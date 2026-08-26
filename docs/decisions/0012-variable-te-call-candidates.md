# 0012 — Store only ambiguous variable-call TE candidates

## Status

Accepted.

## Context

The legacy report generated `MULTIPLE_TE_BY_ID.txt` by grouping two OUTSIDER
`COUNT_TE_IN_RS.txt` files. It was useful because it exposed several family
labels for one insertion ID, but it also included events later excluded from
the public call table, omitted the final group because of an AWK flush bug, and
did not identify which label was ultimately reported.

Writing every unambiguous call into a new candidate table would duplicate
`TE_INFOS.bed` without adding information. It would also make the exceptional
cases harder to find.

## Decision

The refactored workflow writes `TE_CALL_CANDIDATES.tsv`. An event is included
only when all of the following are true:

- it is present in the final `TE_INFOS.bed` call table;
- candidate evidence is still available from a declared upstream output;
- at least two distinct TE labels remain.

Every candidate for a retained event is written, including the family already
reported in `TE_INFOS.bed`. Rows share a stable `candidate_group_id` and contain
the event coordinates, source, reported and candidate labels, assignment,
rank, candidate count, evidence count and fraction, evidence channels and an
explicit ambiguity type.

`evidence_fraction` is the candidate count divided by the total raw candidate
count for that event. It is a descriptive share of preserved evidence, not a
posterior probability or a biological allele frequency.

The current ambiguity type is `unresolved_family_candidates`. Raw candidate
counts contain no coordinates on the inserted sequence, so the table does not
claim whether alternatives are competing annotations of one segment or
separate components of a composite insertion.

## Current scope

The migrated OUTSIDER direct-alignment and Sniffles evidence populates the
table. The historical INSIDER classifier selects one best family before its
public CSV is written, so no trustworthy alternative set survives at this
stage. INSIDER support must be added when that classifier is redesigned to
preserve per-family alignments and their sequence coordinates.

On the current regression fixture the normalized result contains 10 reported
ambiguous calls and 33 candidate rows. The historical file had 36 rows from a
different, partly unreported set; byte compatibility is therefore neither
expected nor desirable.

## Consequences

The file is compact and directly answers “which reported calls remain
ambiguous?”. Header-only output is a valid result when no ambiguity exists or
when the enabled pipelines do not preserve alternative evidence. Quarto shows
the same candidates in a dedicated table and labels the current primary choice
without presenting alternatives as confirmed biological components.
