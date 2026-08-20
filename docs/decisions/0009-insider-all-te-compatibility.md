# 0009 — Migrate whole-assembly INSIDER TE annotation in compatibility mode

## Status

Superseded by ADR 0010. The byte-compatible implementation remains available
in Git history, but it is no longer part of the active workflow.

## Context

`CHOICE.INSIDER_VARIANT.DETECT_ALL_TE` historically activated
`TE_ALL_IN_ASSEMBLY`. Unlike the main INSIDER branch, this step does not detect
structural differences against the reference. It annotates TE-like regions
across the complete `DATA.GENOME` assembly.

The historical rule performed two BLAST passes:

1. BLAST the complete assembly against the TE consensus database;
2. retain chromosomes matching the historical `CHROM_KEEP` expression;
3. merge first-pass intervals separated by at most 100 bases;
4. extract each merged interval and BLAST it against the same TE database;
5. select calls with `global_sv.py` and publish `POSITION_ALL_TE.bed`.

The rule rebuilt an existing BLAST database, hid every intermediate, treated
BLAST coordinates directly as BED coordinates, and left incomplete outputs on
some empty or failed runs. The option remained documented after the initial
refactor but was silently ignored by `workflow/Snakefile`.

## Decision

The compatibility workflow now reads `DETECT_ALL_TE` as a strict boolean. When
it is true, the default `all` target includes `POSITION_ALL_TE.bed`; the named
target `insider_all_te` is also available. Enabling it while the INSIDER
pipeline is disabled is an explicit configuration error.

The migration is split into declared rules for:

- whole-assembly BLAST;
- 100 bp region merge and FASTA extraction;
- candidate-region BLAST;
- legacy `global_sv.py` classification;
- validated BED formatting and publication of the compatibility symlink.

The shared normalized genome, genome index, TE FASTA and BLAST indexes are
reused. Temporary BLAST and region files are tracked by Snakemake and removed
only after their consumers finish. A run with no hit produces two header-only
TSV files, an empty BED and a valid public symlink.

The whole-assembly annotation remains separate from `TE_INFOS.bed` and the
Quarto report. Those artifacts represent variant-associated INSIDER/OUTSIDER
calls; silently mixing resident TE annotation into them would change their
scientific contract and inflate call counts.

## Compatibility oracle

On the historical test fixture, the migrated outputs are byte-identical:

| Output | Lines | SHA-256 |
|---|---:|---|
| `ALL_TE.csv` | 328 | `18ffb72bb4553bce78ffab9b2a75c1ce28d0fd5836be44686a253d32fb2313d9` |
| `ALL_TE_COMBINE_TE.csv` | 332 | `e52c8c0cf9e1f58150824c821f2a4046ed9cb790ce9dbfd0d5b9a20fdb8801a8` |
| `POSITION_ALL_TE.bed` | 327 | `8d737e18c2235c97164b7a65619956f37cf3ffadacd4560618a84c9d724df143` |

The 327 final records cover 64 TE family labels on two chromosomes. The first
and second BLAST passes take approximately ten seconds on the test assembly;
runtime on a large assembly can be orders of magnitude higher.

## Preserved limitations

This decision migrates behavior; it does not endorse the method as a modern
repeat-annotation model.

- Genome-wide BLAST against consensus sequences is expensive and its defaults
  are not explicitly calibrated for fragmented, diverged or nested repeats.
- First-pass BLAST coordinates are consumed directly as BED intervals even
  though BLAST and BED use different coordinate conventions.
- Merging regions within 100 bp can collapse nearby independent, nested or
  composite TE copies before classification.
- `global_sv.py` selects a best family per merged query. It does not preserve a
  structured set of competing or nested annotations.
- The combined alignment must meet the size threshold, but individual rows in
  `ALL_TE.csv` are filtered only by identity. In the oracle, reported
  `size_per` values therefore extend down to 23.15 despite the configured
  80 percent threshold.
- The public four-column BED discards the strand retained in `ALL_TE.csv` and
  carries neither confidence nor intactness/divergence annotations.
- `CHROM_KEEP` remains under the OUTSIDER parameter namespace solely for
  compatibility.

These limitations should be addressed in a separate scientific redesign with
new output schemas and new regression expectations.
