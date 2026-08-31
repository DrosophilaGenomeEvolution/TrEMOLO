# Population frequency trajectories

This module compares TE allele frequencies across samples, generations or
timepoints. Version 2 builds the normalized TrEMOLO locus model first and never
turns a missing call into frequency zero.

## Recommended sample manifest

Use a tab-separated file with one row per sample or replicate:

```text
sample_id	timepoint	replicate	te_infos
G11_R1	11	1	/path/G11_R1/TE_INFOS.bed
G11_R2	11	2	/path/G11_R2/TE_INFOS.bed
G17_R1	17	1	/path/G17_R1/TE_INFOS.bed
```

`timepoint` accepts a number or `G<number>`. Relative `te_infos` paths are
resolved from the manifest directory. Sample IDs must be unique.

The historical list format remains accepted:

```text
/path/work_G11:G11
/path/work_G17:G17
```

It cannot describe replicates and should only be used for compatibility.

## Run

```bash
./module build 1 \
  -i samples.tsv \
  -g reference.fasta \
  -o POPULATION_REPORT
```

Options:

- `--locus-window BP`: maximum complete anchor span of a provisional locus;
  default `20`. Calls cannot chain beyond this span.
- `--trend-epsilon VALUE`: minimum frequency change per timepoint used to call
  increasing/decreasing trajectories; default `0.001`.
- `-c/--chrom` is accepted for command compatibility, but filtering now occurs
  interactively in the report.

The reference FASTA is mandatory. Its SHA-256 checksum becomes the
`reference_id` shared by every input run.

## Outputs

```text
POPULATION_REPORT/
  LOCUS_MODEL/
    manifest.json
    loci.tsv
    alleles.tsv
    components.tsv
    observations.tsv
    evidence.tsv
  population-observations.tsv
  population-trajectories.tsv
  report-data.json
  report.qmd
  report.html
  index.html
```

`population-observations.tsv` uses explicit states:

- `observed`: a compatible call and a frequency are available;
- `observed_unquantified`: a call exists but has no valid frequency;
- `missing`: no compatible call was emitted for that sample/allele.

`missing` is not biological absence. Confirmed empty-site evidence will require
a later population genotyping stage and will receive a distinct status.

Trajectory classes are calculated from replicate means at each observed
timepoint. With fewer than two observed timepoints the result is
`insufficient_data`; changes in both directions are labelled `variable`.

The Quarto report is self-contained, uses no R, CDN, pandas or Node server, and
retains interactive chromosome, TE, trend and observation-count filters.
The refactored Singularity definition includes Quarto; with an older image,
run the wrapper from a host installation of Quarto instead.
