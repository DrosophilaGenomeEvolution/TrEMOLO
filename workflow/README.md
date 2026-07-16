# Refactored workflow

This directory contains the incremental replacement for the legacy `run.snk`.
At this stage it validates the pipeline choices and prepares immutable inputs;
scientific rules still run through the legacy workflow.

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
