# 0002 — Deterministic OUTSIDER TSD candidate output

Status: accepted

## Context

`find_all_type_ins.py` processes chromosomes in parallel. Historically, every
worker wrote its TSD candidates through the same inherited Python file handle.
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

- every legacy TSD record to remain present;
- every `SV_INS.bed` record to have exactly one matching TSD record;
- all other OUTSIDER TE detection outputs to remain byte-identical.

## Consequences

The refactored test output contains 894 deterministic TSD candidates, recovering
20 records relative to `work_test`. TE detection and the 57 merged OUTSIDER
calls are unchanged. The later TSD/breakpoint migration must explicitly test
whether any recovered candidate produces an additional validated TSD rather
than treating a changed TSD result as an unexplained regression.
