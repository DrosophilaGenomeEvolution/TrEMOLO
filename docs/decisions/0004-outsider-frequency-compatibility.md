# 0004 — Preserve OUTSIDER frequency semantics while replacing orchestration

Status: accepted

## Context

The historical `FREQUENCEv2` rule combines 69 OUTSIDER candidates in the test
dataset: 19 direct `TrEMOLO` insertions, 15 Sniffles events, and 35 hard-clipped
events. It then launches several `getFrequency.py` processes in the background
and concatenates their chunk files.

That orchestration can silently reuse stale chunks or accept failed workers.
Several scientific behaviours are also coupled to its output format:

- evidence records are sorted by `event_id:read_id:type`; `E` therefore sorts
  before insertion/deletion evidence for the same read, and the state machine
  ignores the later evidence after marking that read as visited. In `work_test`,
  this masks 48 candidate/read groups: 28 `I`, 18 `D`, and 2 `H` groups;
- the event type is inferred with `event_id.split(".")[1]`. IDs such as
  `HARD.4.R` (and equivalent `SOFT` IDs) consequently yield a numeric token,
  while the aggregation emits only `INS` and `DEL`. The 35 hard-clipped
  candidates in `work_test` are therefore absent from the frequency table;
- aggregation groups records by event ID alone, rather than by the complete
  candidate identity, TE family, and locus. This cannot safely represent
  ambiguous families or several TE alleles at one population locus;
- a candidate with no retained MAPQ 10 read produces no evidence row and is
  consequently absent from the result. The compatibility table cannot
  distinguish zero coverage (`no_call`) from an event that was never tested.

Correcting these behaviours during the workflow migration would mix scientific
changes with orchestration changes and make comparison with the completed
legacy dataset ambiguous.

## Decision

The first modular frequency implementation preserves the historical text
outputs and scientific semantics. Regression checks require byte-identical
frequency inputs and tables, as well as these test-dataset invariants:

- 69 candidates: 19 `TrEMOLO`, 15 `sniffles`, and 35 `HARD`;
- 34 frequency rows: 29 `INS` and 5 `DEL`, with the 35 `HARD` candidates
  omitted as in the legacy implementation;
- 29 precise-frequency rows, all `INS`;
- 25 frequency events represented in `MERGE_TE_ALL.bed` and 9 retained
  frequency events outside that merged-call set.

Only orchestration is corrected in compatibility mode. The refactored DAG
declares all inputs and outputs, processes the complete candidate table without
persistent chunk artifacts, propagates worker failures, streams the externally
sorted aggregation, and never removes shared `.fai` files globally.

## Consequences

The migrated workflow can be evaluated independently of a future frequency
algorithm: a difference in one of the compatibility files is a migration
regression, not an intentional biological change.

The compatibility table must not be interpreted as a complete population
locus model. A later strict mode should give evidence precedence over `E`,
parse source and SV type explicitly, include clipped-read candidates, and group
by a stable locus/candidate key that can retain multiple TE families or alleles
at the same locus. That mode needs separately named outputs and scientific
validation rather than replacement of the compatibility oracle.
