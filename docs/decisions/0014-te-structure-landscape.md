# 0014 — Make TE structures the primary DB explorer

## Status

Accepted.

## Context

The first standalone DB TE report presented insertion sequences as table rows.
A row had to be selected before its structure appeared in a separate panel.
This supported targeted lookup but made visual discovery difficult: users could
not scan structures, notice an unusual component arrangement and then inspect
its evidence in place.

## Decision

The primary report view is a vertically navigable structure landscape. Every
insertion sequence is represented by a card and no scientific subset is
selected by default. The initial order brings multiple components, retained
alternatives and TE matches before unresolved and no-match sequences, but all
remain present.

Each card exposes, without interaction:

- event and genomic position;
- insertion length and reported TE, when available;
- classification, final-call state and component/alternative counts;
- coordinate-scaled component segments on the insertion sequence;
- TE name, strand, consensus coverage and identity;
- up to three rejected candidates when no component survives filtering.

Chromosome, TE family, classification, final-call state, minimum component
count and free-text filters narrow the same landscape. Alternative sort orders
support genomic navigation, component-first review and length inspection.

Detailed normalized matches and raw BLAST HSPs are disclosed within each card.
They include query and consensus coordinates, coverage, identity, score,
assignment, threshold status and rejection reason. These tables are created
only when opened so the complete default landscape does not instantiate every
evidence row in the browser.

The event aggregation remains a secondary cross-sequence view. Activating an
event focuses the corresponding cards in the landscape rather than rendering a
detached selected-item plot.

## Consequences

The normalized TSV and JSON schemas, thresholds, component selection and event
classifications do not change. The report remains self-contained and requires
no JavaScript framework or network dependency. Rendering more cards increases
the initial DOM size, while lazy evidence tables limit the dominant cost. A
future cohort-scale implementation may require viewport virtualization, but
the current standalone module remains small enough for direct browsing.
