# 0006 — Treat empty INSIDER results as valid workflow states

Status: accepted

## Context

The migrated INSIDER workflow historically assumed that every assembly
produced structural variants, that both insertion-like and deletion-like
categories were present, and that every sequence search produced at least one
BLAST row. Three normal biological outcomes were consequently reported as
pipeline failures:

- empty Assemblytics `between` and `within` tables failed before TE detection;
- either an empty insertion FASTA or an empty deletion FASTA failed sequence
  extraction, even when the other category was valid;
- a zero-row BLAST table raised `pandas.errors.EmptyDataError` in
  `global_sv.py` and left incomplete outputs.

Removing only the failure checks is unsafe. `filter_gap_SVs.py` interprets the
first BED record specially and would turn a completely empty input into a
synthetic tab-only record. Hidden per-category FASTA files could also retain
stale sequences after a rerun.

## Decision

Absence is represented by a complete set of valid files, rather than by a
missing file or a successful command that left stale content:

| Stage | Empty-result representation |
| --- | --- |
| Assemblytics variants | zero-byte BED |
| insertion/deletion sequences | zero-byte FASTA per absent category |
| BLAST alignments | zero-byte `.bln` per absent category |
| detailed TE calls | canonical header-only TSV |
| combined TE calls | canonical header-only TSV |
| formatted TE positions | zero-byte BED |
| INSIDER compatibility frequency | canonical header-only TSV |

The workflow applies these rules independently to insertion and deletion
branches:

- gap annotation runs only when the concatenated Assemblytics BED contains a
  record; a true zero-SV result remains zero bytes;
- per-category FASTA artifacts are always truncated before extraction;
- a non-empty category BED must yield sequence data, otherwise the rule fails;
- the historical zero-length Assemblytics interval is excluded with an explicit
  log entry; malformed or reversed coordinates remain fatal errors;
- BLAST is skipped only for a zero-byte query. Any BLAST command invoked for a
  non-empty query must succeed;
- `global_sv.py` emits both canonical headers for a zero-byte BLAST table but
  continues to reject malformed non-empty tables;
- normalized genome and reference FASTA indexes are shared, declared DAG
  outputs rather than files created implicitly by concurrent `bedtools` calls.

## Consequences

An assembly with no SV, a sample with only one SV direction, and a valid search
with no TE hit can all reach `insider_te_detection` successfully. Reruns cannot
reuse a sequence or classification from a previous non-empty result.

Downstream code must distinguish a header-only table from a table containing
calls by counting data rows. File-size tests such as `test -s INSERTION.csv`
are not valid presence checks because a canonical empty table is intentionally
non-zero in size.

The `insider_frequency` target continues to require an OUTSIDER BAM: it is a
single-sample compatibility calculation that combines assembly calls with read
evidence, not an INSIDER-only summary. Avoiding read mapping when there are no
candidates would require a separately designed conditional target under the
Snakemake 5.10 compatibility constraint.
