# TrEMOLO 3.0.0-beta.1 — proposed release

Status: beta release preparation; validation updated 2026-11-10. These notes do not publish a tag,
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
- Freeze original analyses before migration. The 2.x workflow has been removed;
  use Git history for previous implementations. SVIM and caller-disabled OUTSIDER
  modes are retired; select Sniffles 1 or 2 explicitly. See the
  [compatibility audit](v3-compatibility-audit.md).

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
Fresh full Snakemake execution now passes in the local scientific image for
Sniffles 1, Sniffles 2 and a Sniffles 2 no-call run. All three include integration,
Liftoff, resident annotation, unmapped reads and Quarto, and have an idle DAG after
completion. The container suite passes all 245 tests without skips. See the
[dated compatibility audit](v3-compatibility-audit.md) and its
[machine-readable validation record](validation/2026-11-10-migration.json).
This validates a modified working checkout, not a clean release commit or a
fresh rebuild of the image.

See [GFA outputs](te-gfa.md) and [real-run regression procedure](real-run-regression.md).


## Before publishing this beta

- [x] Run the complete tests in the intended scientific image and record its checksum.
- [x] Execute the current checkout on the bundled fixture with a fresh output directory;
  run alignment, integration and report regression checks, then confirm an idle DAG.
- [x] Exercise Sniffles 2 and its no-call case when advertising that backend.
- [ ] Preserve the real-data configuration, input hashes, tool versions and commands;
  the available dm6/G73 graph snapshot alone is not a full-run reference.
- [x] Keep complete Savanache visualization/index/SWAP compatibility outside the
  validated scope until its documented indexing limitation is resolved.
- [ ] Record the final reviewed source commit and publish the matching image and checksums.
  Create an annotated `v3.0.0-beta.1` tag only after reviewing the final release
  commit. No tag or release is created by updating these notes.

## Audit evidence available during preparation

Fresh 2026-11-10 validation supersedes the earlier agent-session container
restriction. It produces 529 calls with Sniffles 1 and 526 with Sniffles 2,
including the no-call case where other evidence channels remain active.
Alignment provenance, integration, report semantics and JavaScript interactions
pass; all three DAGs are idle after completion. The Sniffles 2 regression checker
now expects the native caller IDs already preserved by the adapter.

The dm6/G73 GFA snapshot comparison remains separate from these fresh fixture
runs. The 2026-09-11 clean-run record is historical and predates later corrections.
A fresh build from the current recipe and a clean release commit remain untested;
some recipe dependencies are unpinned. Distribute the tested image or revalidate
any newly built image before publishing its checksum.

Earlier preparation checks (before the fresh container runs):
Preparation checks for the documentation/version change: 46 focused tests ran
successfully (one Snakemake-dependent test skipped), including source identity,
archive fallback, module/GFA CLI version agreement, graph reconstruction,
representation summaries and both report builders. The new population manifest
builds six samples into 214 provisional loci and 215 trajectories. The new
Snakemake startup hook still requires execution in the scientific environment.
