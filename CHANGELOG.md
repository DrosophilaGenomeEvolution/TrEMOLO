# Changelog

## 3.0.0-beta.1 — in preparation

Proposed beta; no publication date or release tag assigned.

- Automatically decompress gzip FASTQ inputs during preparation when OUTSIDER
  clipped-read analysis requires random-access indexing; preserve source files
  and record the preparation method in the input manifest.

- Remove the 2.x generated workflow (`run.snk`, `rules.snk`, `template.snk`)
  and its private `env.sh`; retain shared scientific helpers and historical oracles.
- Retire SVIM and caller-disabled OUTSIDER modes; document the compatibility audit.
- Rename development-labelled configuration and documentation for v3.
- Migrate the supported entry point to `workflow/Snakefile` and Quarto reports.
- Correct OUTSIDER evidence extraction, TSD candidate selection, insertion
  integration and Liftoff auditing; scientific outputs can differ from 2.x.
- Add resident TE and unmapped-read catalogues, explicit candidate provenance,
  configurable windows and Sniffles 2 support.
- Replace module report interfaces and add reference-backed insertion GFA export.
- Add contextual offline report help, detector/event counts, chromosome-window
  heatmaps with filter-impact tooltips, and independent structure threshold previews.
- Add an interactive frequency-versus-depth report view with matched support
  counts, linear/log scales and adjustable guides; preserve original frequencies.
- Identify the version, source commit and modified-checkout state at startup.
- Reorganize user documentation, minimal configurations and interpretation guides.

See [full beta notes, limitations and validation status](docs/release-notes-functional-candidate.md).
Historical 2.x release notes remain available in the repository history.
