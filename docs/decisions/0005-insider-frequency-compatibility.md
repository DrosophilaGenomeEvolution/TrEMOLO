# 0005 — Preserve INSIDER empty-site frequency semantics during migration

Status: accepted

## Context

The historical `FREQ_INSIDERv2` rule estimates read support for 201 INSIDER
assembly calls in `work_test`. It combines interval construction, two
`samtools view -L` passes, CIGAR parsing, flank depth, and aggregation in one
shell block. Most of its intermediate files are hidden from Snakemake, worker
resources are not declared, stale files can be reused, and the rule removes
every `.fai` below the work directory.

The resulting table also has scientific behaviours that must not be changed
silently during an orchestration migration:

- the two BAM filters require an alignment to overlap any outer probe and then
  any inner probe. They do not demonstrate coverage of both sides of one
  candidate;
- supplementary alignments are excluded by default, while secondary
  alignments remain. Supplying any custom numeric `-F` option disables the
  implicit `-F 2048` rather than adding to it;
- an empty-site read is a distinct query name with a CIGAR deletion whose
  length is within 30 bases of the assembly-call length. The right edge of the
  deletion must be strictly less than 30 bases from either breakpoint;
- the compatibility CIGAR cursor advances only for `M`, `D`, and `=`. It does
  not advance for `N` or `X`;
- missing flank depth is treated as zero. The two values are averaged with
  integer truncation, so one covered flank can yield a positive depth while a
  depth of one on only one side becomes zero;
- TE-support depth is the mean flank depth minus the empty-site count, without
  clamping. Zero total depth is reported as frequency `0.0000`, not as a
  `no_call`.

On the fixture these behaviours yield 201 rows: 190 at depth zero, ten at
depth one, and one at depth two. Ten rows have frequency 100%, and one covered
event is supported only by an empty-site deletion. The reference output is
byte-identical when recomputed from the refactored mapping BAM.

## Decision

The first modular implementation preserves the historical table and its
single-sample arithmetic. It splits the calculation into explicit rules for:

1. candidate and probe BED preparation;
2. outer and inner BAM subsets with declared indexes and threads;
3. empty-site deletion counting;
4. flank-depth measurement;
5. deterministic table aggregation.

The persistent `INSERTION_TE.bed`, `DEL_NB.bed`, `DEPTH_FK.bed`,
`DEPTH_FK.txt`, and `DEPTH_TE_INSIDER.csv` tables are regression oracles. The
malformed probe rows historically produced from the CSV header are removed;
they never selected an alignment and therefore do not alter a scientific
output. Probes that would start before coordinate zero are omitted instead of
being emitted as invalid BED records; their flank remains missing and therefore
has the historical depth value zero. Header-only candidate input now creates
fresh empty intermediates and a final header-only table. The artificial
dependency on `MERGE_TE_ALL.bed` is removed because that file was never read.

## Consequences

The migrated target is deterministic, propagates failures, exposes scheduler
resources, and does not delete unrelated indexes. A change to a compatibility
table is therefore treated as a migration regression rather than an implicit
algorithm update.

`DEPTH_TE_INSIDER.csv` is not a population allele-frequency model. It has no
`sample_id`, stable `locus_id`, `allele_id`, genotype likelihood, haplotype, or
quality state. Several TE-family calls at the same locus remain separate rows,
but their spanning and empty-site reads are not assigned to competing alleles;
the same evidence can consequently support more than one row.

The frequency sub-DAG accepts a previously generated header-only
`INSERTION.csv`, but the current upstream INSIDER TE-detection rules still
require non-empty insertion and deletion FASTA/BLAST results. A completely
TE-free sample therefore needs a separate upstream empty-result hardening step
before it can reach this target naturally.

A future strict population mode must use the
[locus data model](../locus-data-model.md), preserve all TE alleles per locus
and sample, require independently validated breakpoint-side evidence,
distinguish zero coverage from a measured absence, consume complete CIGAR
reference semantics, and emit separately named outputs. That scientific change
requires new truth sets and must not replace this compatibility oracle in
place.
