# 0003 — Preserve OUTSIDER TSD compatibility before strict breakpoint validation

Status: accepted

## Context

The historical `GET_SEQ_TE` rule contains two order-dependent row drops:

- `awk 'NR>1'` discards the first matching direct insertion even though the
  `grep` output has no header (`TrEMOLO.INS.72` in `work_test`);
- another `awk 'NR>1'` discards the first combined record (a terminal 412
  deletion in `work_test`).

It also extracts `SIZE_FLANK + 1` bases on the left side of a read but exactly
`SIZE_FLANK` on the right. With a configured size of 10, this produces 11/10
base read flanks. The TSD caller searches for any exact common motif of at
least four bases inside those windows; only four of the nine reference TSDs
are direct left-suffix/right-prefix junctions.

Changing these behaviours while migrating the DAG would mix workflow and
scientific changes. It would also make differences from the completed
`work_test` oracle difficult to attribute.

## Decision

The first modular implementation preserves the two selection quirks, the
11/10 flank convention, and the existing TSD algorithm. They are isolated in
the preparation script and covered by byte-for-byte regression checks.

Every selected record is additionally written to
`OUTSIDER/FK/TSD_FLANK_ELIGIBILITY.tsv`. The report records whether read and
genome coordinates are usable and gives an explicit reason for rejection.

`TrEMOLO.INS.72` is not restored in compatibility mode. A reconstructed test
shows that it would add a nineteenth usable flank but no tenth TSD.

## Consequences

The migrated reference remains 18 flank records and 9 TSD calls, exactly
matching `work_test`. Workflow cleanup is therefore verifiable independently
of future changes.

A later strict mode should use explicit merged-call membership, symmetric
flanks, and junction-anchored motif criteria. Its output must be versioned or
labelled separately and evaluated as a scientific algorithm change rather
than accepted through compatibility hashes.
