# 0017 — Recover alignment candidates and expose workflow durations

Status: accepted, 2026-10-07.

## Extraction and performance

CIGAR extraction uses reference blocks (default 10,000,000 bp; configurable via
`PARAMS.OUTSIDER_VARIANT.TE_DETECTION.CHUNK_SIZE`, with 0 meaning chromosomes).
Workers own alignments by their reference start, including operations extending
past a block boundary. Results merge in reference-header/block order, independent
of scheduling. Soft grouping happens after block merging, so nearby evidence can
cross block boundaries and IDs remain stable. Chromosome IDs retain header order.

A shared local multiprocessing event replaces per-read Manager RPCs. When no time
limit is configured, workers never poll the event. Waiting is blocking and worker
exceptions propagate. Exceeding TIME_LIMIT fails the job; partial candidates are
never published as completed outputs. Scratch files are isolated and cleaned.

Hard evidence is streamed instead of buffered with full sequences per chromosome.
HARD fragment preparation reads the extracted FASTQ once, retains only needed
fragments, and writes the existing FASTA/index audit artifacts in the same pass.
INS formatting frees completed clusters and uses a set for read-name membership.
The BLAST classifier uses set membership and scalar accesses instead of list
searches and repeated Pandas row construction, preserving its selection semantics.
FASTQ indexing is skipped when clipped detection is disabled.
The second all-soft-reads BLAST has no downstream consumer and is now optional
(`SOFT_CLIP_DIAGNOSTICS`, default false). Its declared output is empty when disabled.

## Correctness changes

The historical extractor reset soft buffers after each clip, losing the evidence
before the flushing threshold. Soft evidence is now retained and grouped within
the configured anchor window. The sequence used for each support record is the
actual current clip, and longest representative sequences are tracked by side.

INS cluster output now includes its final cluster and read-position provenance
for the first read of every cluster. HARD output includes each cluster's first
candidate, singleton clusters, and the final cluster. Empty INS output is valid.
Disabling clips truncates old clipped artifacts rather than touching stale files.

These changes deliberately alter OUTSIDER calls. On the bundled Sniffles 1 fixture,
raw CIGAR INS, TSD, HARD evidence and extracted HARD variants remain byte-identical
to the historical fixture. Corrected downstream calls total 529: 377 unchanged
INSIDER calls and 152 OUTSIDER calls (previously 57). OUTSIDER types are 23 INS,
2 DEL, 119 HARD and 8 SOFT. More retained evidence is not an accuracy benchmark.
The old exact-output regression remains frozen as the migration oracle; use
`check_alignment_candidates.py` to validate BAM provenance and complete clusters.
The frequency formula and historical HARD fragment slicing semantics are unchanged.

## Execution-time report

Report schema 1.3.0 adds `timings`: measured wall seconds per benchmark row, branch,
measurement number and source filename. Benchmark files are explicit report inputs
selected from enabled rules, excluding disabled branches and report generation.
Their dependencies also ensure that enabled integration/Liftoff jobs finish before
report preparation. Input preparation now has its own benchmark.

The standalone Quarto report shows a duration bar chart with branch filtering,
longest-first/name ordering and TSV export. Units adapt between seconds, minutes
and hours, while the TSV preserves exact seconds. Job-time sums are explicitly
not elapsed run time because jobs may overlap. Resumed runs can contain previous
measurements. Missing timings have an empty-state view; malformed measurements fail
rather than being silently presented as zero. No external chart dependency is used.

## Validation

- 120 workflow unit tests pass in `TrEMOLO-3.simg`, including cross-block grouping,
  read ownership, singleton/final clusters, worker errors, timeout, stale outputs,
  direct HARD fragment extraction and duration parsing.
- A fresh 70-job run completes with both pipelines, resident annotation,
  integration, Liftoff and Quarto. Its report includes 66 measured jobs.
- Additional OUTSIDER-only (clipping disabled) and TE_GENOME-only runs complete;
  INSIDER-only DAG construction passes. Optional empty inputs are grouped to avoid
  a Snakemake 5.10 duplicate-output error when both reads and reference are absent.
- Chrome verifies all 66 bars, branch filtering, name ordering and TSV export.
- BAM provenance check: 894 INS records, 432 SOFT records and 867 HARD records;
  complete INS/HARD cluster emission. Report projection regression passes.
- Optimized classifier output matches the old classifier byte for byte on the
  same historical INS/HARD BLAST inputs.
- Three-run synthetic test, 50,000 reads, four workers, clipping disabled:
  median extraction 3.698 s before, 0.356 s after (about 10.4x). INS and TSD bytes
  are identical. This isolates extraction and is not a whole-pipeline speedup;
  the small real fixture does not establish large-dataset performance.

Local artifacts are in `../optimization_validation_20261007_i_sbq428/`, including
`final_config.yml`, `final_run.log`, `final_output/`, `unit-tests.log`, and
`performance.json`. They are not repository-tracked scientific fixtures.
