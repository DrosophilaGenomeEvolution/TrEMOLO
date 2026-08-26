# Refactored workflow

This directory contains the incremental replacement for the legacy `run.snk`.
It currently prepares immutable inputs and contains migrated OUTSIDER variant
calling, TE detection, frequency, and TSD rules, plus INSIDER structural-
variant calling, TE detection, empty-site frequency, and TSD rules. It also
builds the final legacy-compatible `TE_INFOS.bed` call table and a standalone
interactive Quarto report.

When Snakemake is available, run the preparation DAG with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile test/tmp_config.yml \
  --cores 1
```

The preparation writes normalized, stable filenames under `WORK_DIRECTORY/INPUT`
and records original paths and SHA-256 checksums in `input_manifest.json`.

The first migrated scientific segment can be inspected or executed explicitly:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile test/tmp_config.yml \
  --cores 8 outsider_variant_calling
```

It covers Minimap2 indexing and mapping, BAM preparation, mapping statistics,
and Sniffles 1 variant calling. The default `all` target now completes every
enabled scientific branch through `TE_INFOS.bed`; it additionally renders the
report when `CHOICE.PIPELINE.REPORT` is true.

The migrated INSIDER structural-variant segment is available with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 insider_variant_calling
```

Shared FASTA and BLAST indexes for the normalized TE database are built once:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 1 te_database_indexes
```

INSIDER TE detection reuses these indexes and can be run with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 insider_te_detection
```

Resident whole-assembly TE annotation is an independent optional target:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 te_genome
```

It is included in the default target when `CHOICE.PIPELINE.TE_GENOME` is true;
the old `insider_all_te` target and `DETECT_ALL_TE` choice remain temporary
aliases. `PARAMS.TE_GENOME.TARGETS` selects `GENOME` and optionally `REFERENCE`.
One TE-versus-assembly BLAST is converted into audited fragments, collinear
matches, physical copies, alternative family assignments, provisional
relations, BED and GFF3. The permissive discovery thresholds are independent
from full-length classification. See
`docs/decisions/0010-te-genome-resident-annotation.md`.

This resident annotation is deliberately not merged into the
variant-associated `TE_INFOS.bed`. When both `TE_GENOME` and `REPORT` are
enabled, its audited tables feed a separate resident-catalog section in the
Quarto report.

Insertion and deletion branches are independent. A sample with no qualifying
SV, no sequence in one category, or no BLAST hit now produces canonical empty
or header-only outputs instead of failing or reusing stale artifacts. The
empty-result contract is recorded in
`docs/decisions/0006-insider-empty-result-contract.md`.

When both INSIDER calls and OUTSIDER read mappings are available, the
legacy-compatible INSIDER empty-site frequency table can be generated with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 insider_frequency
```

This target replaces `FREQ_INSIDERv2` with declared interval, BAM, deletion,
depth, and aggregation steps. It preserves `DEPTH_TE_INSIDER.csv` while
propagating command failures and retaining shared FASTA indexes. Its
single-sample and multi-allele limitations are documented in
`docs/decisions/0005-insider-frequency-compatibility.md`.

INSIDER genome flanks and legacy-compatible TSD calls are available with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 insider_tsd
```

This target replaces `TSD_INSIDER`, declares the normalized genome index, and
preserves byte for byte the 377 candidate pairs and 210 historical TSD calls.

OUTSIDER TE detection is split into four explicit evidence sources: CIGAR
insertions, Sniffles calls, soft-clipped reads, and hard-clipped reads. To stop
after the two primary sources, run:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 outsider_primary_te_detection
```

To include clipped-read evidence and reproduce the historical 100 bp merge:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 outsider_te_detection
```

All four branches reuse the single normalized TE FASTA and its shared indexes.
Unlike the legacy rules, they declare their intermediate outputs, do not append
concurrently to `SV_SIZE.tsv`, and do not delete unrelated FASTA indexes. Raw
CIGAR flank candidates are also written per chromosome before deterministic
merging, avoiding the data loss caused by concurrent writes to one file.

Legacy-compatible OUTSIDER read frequencies are available with:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 outsider_frequency
```

This target replaces the background chunk processes from `FREQUENCEv2` with a
declared, deterministic DAG and propagates worker failures. It preserves the
historical `FREQUENCY_TE_INS.tsv` and `FREQUENCY_TE_INS_PRECISE.tsv` bytes,
including the current candidate-selection semantics; the known compatibility
limits are recorded in
`docs/decisions/0004-outsider-frequency-compatibility.md`.

OUTSIDER read/genome flanks and legacy-compatible TSD calls are available as a
separate target:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 outsider_tsd
```

This target replaces the monolithic `GET_SEQ_TE` and `TSD_OUTSIDER` rules. It
declares every sequence, table, and merged-call dependency; preserves the
historical `ALL_FK_REPORT_FT*.bed` and `TSD_TE.tsv` bytes; and writes
`OUTSIDER/FK/TSD_FLANK_ELIGIBILITY.tsv` with one acceptance or rejection reason
per candidate. It does not delete shared `.fai` files.

The public call table for the enabled INSIDER/OUTSIDER pipelines can then be
generated without invoking the legacy R Markdown report or its network
requests:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 te_infos
```

`TE_INFOS.bed` keeps one row per detected event—including distinct TE calls at
the same locus—and reproduces the historical 15-column table exactly. The
compatibility choices and population-model boundary are documented in
`docs/decisions/0007-insider-tsd-te-infos-compatibility.md`.

Ambiguous family evidence for reported variable calls is normalized separately:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 1 te_call_candidates
```

`TE_CALL_CANDIDATES.tsv` deliberately omits unambiguous calls. For each retained
event it contains every candidate, identifies the family currently reported in
`TE_INFOS.bed`, ranks alternatives, and records evidence counts, fractions and
channels. The migrated OUTSIDER evidence is populated now; INSIDER alternatives
require a later classifier redesign because the compatibility classifier keeps
only one family per event.

The new report is generated directly with Quarto:

```bash
snakemake --snakefile workflow/Snakefile \
  --configfile tests/workflow/refactor_config.yml \
  --cores 8 report
```

It produces `REPORT/report-data.json`, the generated `REPORT/report.qmd`, and a
single self-contained `REPORT/report.html`. The report has interactive source,
type, chromosome, family, TSD and frequency filters; genome and family views;
the complete call table; and an explicitly provisional nearby-call view. It
also exposes ambiguous variable-call candidates, the resident-copy tiers, a
browser-only threshold calibration view, and all alternative TE assignments
when `TE_GENOME` is enabled. It does not use R, `curl`, a web server, or external
JavaScript/CDN assets.

The refactored container definition pins Quarto 1.9.36. The previously built
Snakemake 5.10 image does not contain it: that image can generate the QMD/data,
but the source must be rendered by a host Quarto installation or by a rebuilt
container. Configure a non-default executable with `TOOLS.QUARTO`; customize
the report with `REPORT.TITLE`, `REPORT.AUTHOR`, and `REPORT.LOCUS_WINDOW`.

Migration equivalence against a completed legacy work directory is checked by:

```bash
python3 tests/regression/check_variant_calling.py \
  ../work_test ../work_refactor_test
```

The INSIDER caller architecture decision is recorded in
`docs/decisions/0001-insider-sv-backends.md`.
The OUTSIDER frequency compatibility decision is recorded in
`docs/decisions/0004-outsider-frequency-compatibility.md`.
The INSIDER frequency compatibility decision is recorded in
`docs/decisions/0005-insider-frequency-compatibility.md`.
The INSIDER empty-result contract is recorded in
`docs/decisions/0006-insider-empty-result-contract.md`.
The INSIDER TSD and final call-table compatibility decision is recorded in
`docs/decisions/0007-insider-tsd-te-infos-compatibility.md`.
The direct Quarto report migration is recorded in
`docs/decisions/0008-quarto-report.md`.
The resident-catalog report integration is recorded in
`docs/decisions/0011-resident-te-quarto-report.md`.
The variable-call ambiguity table is recorded in
`docs/decisions/0012-variable-te-call-candidates.md`.
The whole-assembly INSIDER TE annotation decision is recorded in
`docs/decisions/0009-insider-all-te-compatibility.md`.
