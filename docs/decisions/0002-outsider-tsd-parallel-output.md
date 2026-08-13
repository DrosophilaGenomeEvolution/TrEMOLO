# 0002 — Deterministic OUTSIDER CIGAR flank-candidate output

Status: accepted

## Context

`find_all_type_ins.py` processes chromosomes in parallel. Historically, every
worker wrote its raw flank candidates through the same inherited Python file
handle.
The resulting `INS_FOR_TSD.txt` was neither ordered nor safe for records longer
than an atomic operating-system write.

The `work_test` oracle contains 894 CIGAR insertion records in `SV_INS.bed`, but
only 874 corresponding records in `INS_FOR_TSD.txt`. Sixteen candidates from
2L and four candidates from X were lost. The remaining records also appear in
a scheduler-dependent order.

## Decision

Each chromosome worker writes to a private temporary TSD file. The parent
process merges those files in BAM reference order after all workers finish.

Migration regression requires:

- every legacy raw flank-candidate record to remain present;
- every `SV_INS.bed` record to have exactly one matching candidate record;
- all other OUTSIDER TE detection outputs to remain byte-identical.

## Consequences

The refactored test output contains 894 deterministic raw candidates,
recovering 20 records relative to `work_test`. TE detection and the 57 merged
OUTSIDER calls are unchanged.

The historical workflow never reads `INS_FOR_TSD.txt`: `GET_SEQ_TE` rebuilds
its flanks from the clustered insertion FASTA and the combined BLAST tables.
Consequently, recovering these records cannot alter `TSD_TE.tsv` in the
compatibility pipeline. Connecting them to TSD detection would be a separate
scientific change with its own validation.
