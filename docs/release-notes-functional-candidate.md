# Functional release candidate

This candidate brings together the migrated workflow and the first explicit
post-processing command for reference-backed TE insertion graphs. No release
number or tag is assigned by these notes.

## Available behavior

- Main and ambiguous-call report tables have 50 rows per page, numbered navigation
  and direct access to distant pages. Candidate fractions and sequence provenance
  are explained in the report.
- Candidate sequence provenance distinguishes read-extracted insertions from
  VCF-reported sequences and preserves raw alignment evidence.
- GFA representatives use collinear, non-overlapping BLAST bitscores; reverse
  copies use oriented walks. Synthetic levels expose local family alternatives.
- GFA export supports TE-family selection with all admissible alternatives at
  matching loci and a chromosome regex. By default it exports only chromosomes
  present in TE_INFOS; `--all-reference-chromosomes` disables that restriction.
- Each new GFA build writes an HTML/JSON representation balance and per-call
  outcomes, distinguishing original calls from graph representatives.
- Unmapped-read selection uses a bounded memory index and a batched disk fallback.
  Optional samtools collate selection preserves secondary/supplementary exclusions;
  the faster two-pass method remains the default based on the local synthetic test.
- Real-run regression capture records input hashes, exact configuration, tool
  versions and result signatures when these source records are supplied.

## Interpretation and limits

The GFA exporter currently handles OUTSIDER insertion sequence channels. INSIDER
deletions, resident TE alleles, nested/composite allele construction and observed
haplotypes are not encoded. The BED limits accepted events; sequence families and
representatives are evaluated from their raw BLAST evidence, so level 1 need not
match the primary BED family. Walk levels are synthetic, unphased and do not
encode biological frequency. Candidate fractions are shares of preserved evidence,
not allele frequencies or proportions of unique reads.

## Validation records

The user has exercised the pipeline on real datasets. The available dm6/G73 GFA
outputs have a frozen output snapshot. Their original full execution metadata
are still required to turn that snapshot into a reproducible full-run reference.
The existing historical regression baseline remains separate and is not replaced.
Tests cover graph reconstruction/orientation, family and chromosome filters,
per-call summaries, row-order-independent comparisons and extraction equivalence.
Full Snakemake execution must be checked in the intended container environment;
the local development interpreter does not have Snakemake installed.

See [GFA outputs](te-gfa.md) and [real-run regression procedure](real-run-regression.md).
