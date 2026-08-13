# Refactored workflow

This directory contains the incremental replacement for the legacy `run.snk`.
It currently prepares immutable inputs and contains migrated OUTSIDER variant
calling plus INSIDER structural-variant calling and TE detection rules.

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
and Sniffles 1 variant calling. The default `all` target remains input
preparation only during the incremental migration.

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

Migration equivalence against a completed legacy work directory is checked by:

```bash
python3 tests/regression/check_variant_calling.py \
  ../work_test ../work_refactor_test
```

The INSIDER caller architecture decision is recorded in
`docs/decisions/0001-insider-sv-backends.md`.
