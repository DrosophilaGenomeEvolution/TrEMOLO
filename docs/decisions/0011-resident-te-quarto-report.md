# 0011 — Project the resident TE catalog separately in Quarto

## Status

Accepted.

## Context

`TE_GENOME` now distinguishes BLAST fragments, consensus matches, physical
copies, alternative family assignments and provisional geometric relations.
Keeping these results entirely outside the report made threshold assessment
and multi-TE ambiguity difficult, while merging them into `TE_INFOS.bed` would
incorrectly mix assembly annotation with variant-associated INSIDER/OUTSIDER
calls.

## Decision

When `CHOICE.PIPELINE.TE_GENOME` and `CHOICE.PIPELINE.REPORT` are enabled, the
four audited resident tables are explicit inputs of `prepare_quarto_report`.
The versioned JSON model contains:

- every physical resident copy and its configured classification tier;
- the complete candidate list, with per-candidate coverage and identity, for
  every ambiguous copy;
- aggregate fragment rejection reasons and provisional relation types;
- the configured discovery and classification thresholds.

The Quarto document renders these data in a separate **Resident TE catalog**
section. Interactive sliders simulate `full_length`, `partial` and
`degraded_relic` classification in the browser. They never rewrite pipeline
outputs and do not simulate discovery thresholds, because rejected HSPs cannot
become physical copies without rerunning the annotator. The scatter plot may be
deterministically thinned above 2,500 points, but all counters and tables use
the complete catalog.

The explorer defaults to copies carrying more than one distinct TE label and
shows the primary assignment alongside every retained alternative. A separate
filter also includes multiple matches to the same label. Primary rank remains
a deterministic display choice rather than a biological assertion of
uniqueness.

## Consequences

`TE_INFOS.bed`, its charts, frequencies and regression oracle remain unchanged.
With `TE_GENOME` disabled, the report remains valid and explicitly states that
the resident catalog is unavailable. With multiple targets, target labels keep
GENOME and REFERENCE annotations separate.

The embedded report grows with the number of physical copies. The JSON avoids
embedding every match for unambiguous copies and includes detailed match rows
only where they are necessary to explain ambiguity. Very large pangenome
catalogs will eventually need a paged or external data layer; that belongs to
the future module and cohort redesign rather than this standalone report.
