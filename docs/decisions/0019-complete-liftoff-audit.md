# 0019 — Keep every integrated INS in the Liftoff audit

Status: accepted.

## Problem

The projection audit was derived only from the mapped GFF. An integrated
insertion with neither flank mapped disappeared from `LIFT_OFF_AUDIT.tsv`,
although it remained present in both reconstructed population genomes.

## Decision

The summarization rule now declares the canonical integrated BED as an input
and passes it through the required `--positions` argument. Its event list is
the source of truth for audit completeness. Each integrated INS receives exactly
one audit row, in canonical BED order, with its original identifier and family.
Starting from the BED also handles an insertion occupying an entire contig,
for which no source GFF flank can be generated.

An event with no mapped flank is `rejected / no_mapped_flanks`, with empty
projected coordinates. A partially mapped or ambiguous pair retains the
existing `missing_or_ambiguous_flank` reason. Existing chromosome, orientation,
ordering and gap checks still determine acceptance of mapped pairs. Projection
rejections are counted in the job log; they are not integration rejections.

Duplicate expected identifiers or mapped identifiers absent from the expected
BED cause an error, avoiding inconsistent audit and public output sets.

## Consequences and validation

No insertion is removed from a reconstructed genome because of mapping failure.
Only accepted projections appear in the public reference BED. On existing
outputs with all events represented in the mapped GFF, the public position set
does not change; its order can change to follow the canonical BED.

Regression validation checks one audit row per integrated INS, family retention,
`no_mapped_flanks` for events absent from the mapped GFF and correspondence of
accepted audit events with the public reference BED. Unit tests cover mixed
complete/partial/absent mappings, a completely empty mapped GFF, an insertion
without source flanks, no expected events, unexpected mappings and duplicate
expected identifiers.
