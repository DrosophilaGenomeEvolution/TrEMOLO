# TrEMOLO 3.0.0-beta.1 — proposed release

Status: beta release preparation, 2026-11-04. These notes do not publish a tag,
GitHub release or downloadable image. `VERSION` contains `3.0.0-beta.1`.

The major version reflects the migrated workflow/report interfaces and
intentional changes to scientific outputs since 2.x. The beta suffix signals
remaining validation and compatibility limits. This is not a claim that
results match the historical migration oracle.

## Migration from 2.x

- Use `workflow/Snakefile`; start in a new output directory.
- Use the current container recipe/image, with Quarto and separate Sniffles 1/2
  environments. The historical 2.5.4b image is not the beta environment.
- Choose Sniffles explicitly; callers and defaults can produce different results.
- Corrected CIGAR/clip candidate extraction, TSD selection, insertion-only genome
  integration and Liftoff auditing intentionally change downstream calls.
- Detailed `PUBLIC/` exports restore original TE names; internal files retain
  prepared identifiers. Frequency tables now have documented headers.
- Module commands build offline Quarto reports. Module 2 requires its adjacent
  `report-assets/`; the old Node server command is retired.
- Freeze original analyses before migration; shared scripts mean this checkout's
  `run.snk` is not a frozen 2.x implementation.

## Available behavior

- Startup identifies version, full source Git commit and modified-checkout state
  for the workflow, module launcher and GFA exporter. Snakemake execution records
  the latest source identity in `log/run-version.json`; archives without Git
  explicitly report an unavailable commit.
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


## Before publishing this beta

- Run the complete tests in the intended scientific image and record its checksum.
- Execute the current checkout on the bundled fixture with a fresh output directory;
  run alignment, integration and report regression checks, then confirm an idle DAG.
- Exercise Sniffles 2 and its no-call case when advertising that backend.
- Preserve the real-data configuration, input hashes, tool versions and commands;
  the available dm6/G73 graph snapshot alone is not a full-run reference.
- Keep complete Savanache visualization/index/SWAP compatibility outside the
  validated scope until its documented indexing limitation is resolved.
- Record the tested source commit and publish the matching image and checksums.
  Create an annotated `v3.0.0-beta.1` tag only after reviewing the final release
  commit. No tag or release is created by updating these notes.

## Audit evidence available during preparation

The available dm6/G73 GFA snapshot comparison and report pagination check pass.
The existing fixture Quarto report passes its semantic checker (529 calls,
58 families, 15 proximity groups); this is not evidence of a fresh execution of
this release commit. Local full-suite execution is limited by missing host
scientific dependencies and the local container runtime installation. The
2026-09-11 clean-run record is historical and predates later corrections.


Preparation checks for the documentation/version change: 46 focused tests ran
successfully (one Snakemake-dependent test skipped), including source identity,
archive fallback, module/GFA CLI version agreement, graph reconstruction,
representation summaries and both report builders. The new population manifest
builds six samples into 214 provisional loci and 215 trajectories. The new
Snakemake startup hook still requires execution in the scientific environment.
