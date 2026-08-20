# 0008 — Replace the refactored report directly with Quarto

## Status

Accepted during the Snakemake 5.10 compatibility migration.

## Context

The historical `REPORT` rule mixed construction of `TE_INFOS.bed`, report-data
aggregation, JavaScript generation, R Markdown mutation with `sed`, package
installation, Bookdown rendering and cleanup of scientific intermediates. It
also queried a remote author service with `curl`. Consequently the report was
not reproducible offline, undeclared source changes did not reliably invalidate
the rule, and report generation could modify or delete pipeline outputs.

The legacy report remains available on the old workflow and can therefore be
used as a visual and semantic reference while the refactored workflow adopts a
clean report architecture directly.

## Decision

The refactored workflow does not migrate the R Markdown implementation. It
generates:

- `REPORT/report-data.json`, a versioned, deterministic presentation model;
- `REPORT/report.qmd`, a complete generated Quarto source with data, CSS and
  JavaScript embedded;
- `REPORT/report.html`, a standalone HTML document rendered by Quarto.

The data builder consumes the already migrated, byte-compatible
`TE_INFOS.bed`, the normalized genome index and input manifest. OUTSIDER mapping
statistics and the SV VCF are declared only when OUTSIDER is enabled. The
report does not recalculate or alter calls, TSDs, sizes or frequencies.

INSIDER empty-site frequency requires the OUTSIDER read mapping by design. In
an INSIDER-only configuration it is therefore omitted rather than silently
triggering the disabled OUTSIDER branch; the corresponding legacy-compatible
frequency field remains unavailable.

The interactive interface uses dependency-free browser JavaScript embedded in
the Quarto document. Vue is not added at this stage: the filters, plots, table
and TSV download do not require an application framework, and avoiding another
runtime keeps the report offline and portable. A component framework can be
introduced later if cohort editing or server-backed interaction warrants it,
without changing `report-data.json`.

Quarto is pinned to version 1.9.36 in the container definition. The existing
pre-built image does not contain Quarto; it can still prepare the report source,
which may be rendered with a host Quarto installation. A rebuilt image runs the
whole target internally.

## Compatibility contract

Report equivalence is semantic, not byte-for-byte HTML equivalence:

- every `TE_INFOS.bed` row appears once in `report-data.json` and the call
  explorer;
- call, source, chromosome, family, type, strand and TSD summaries are derived
  directly from that table;
- the displayed frequency preserves the historical preference for
  `FREQ_WITH_CLIPPED`, falling back to `FREQ`;
- the frequency landscape preserves the historical position-versus-frequency
  view, one chromosome at a time, with points grouped by the legacy `TYPE`
  field and all values driven by the current report filters;
- the report embeds the input `TE_INFOS.bed` SHA-256 checksum;
- the HTML loads no external script or stylesheet and performs no network
  request;
- report generation never removes INSIDER, OUTSIDER or other scientific files.

The regression checker recomputes the report summaries independently from
`TE_INFOS.bed`. The HTML itself is not compared with the Bookdown output because
the interface is deliberately replaced.

The historical timeline parsed scheduler log text and had no declared timing
inputs. It is not reproduced in this first Quarto report. A later execution
profile must consume a declared benchmark manifest rather than scrape log-file
ordering.

## Nearby calls and population loci

The report exposes complete-span groups of breakpoint anchors within
`REPORT.LOCUS_WINDOW` (100 bp by default). These rows are explicitly named
**proximity candidates**. They demonstrate that several calls, TE families or
evidence sources can coexist in one genomic neighborhood, and no row is
collapsed in `TE_INFOS.bed`.

Proximity is not sufficient to infer a shared biological locus, alternative
alleles, a nested TE or a composite insertion. Cohort reports must use the
normalized locus model and cross-sample evidence described in
`docs/locus-data-model.md`; the presentation layer must not silently promote
these candidates to resolved loci.

## Consequences

When `CHOICE.PIPELINE.REPORT` is true, the default `all` target now includes the
complete scientific chain and `REPORT/report.html`. With reporting disabled,
`all` still completes the enabled scientific pipeline through `TE_INFOS.bed`.
The explicit `report` target is always available.
