# 0022 — Expose genomic analysis windows in configuration

Status: accepted.

## Problem and decision

Migrated BEDTools commands and several Python entry points still used fixed
clustering distances, BAM windows and frequency margins. `REPORT.LOCUS_WINDOW`
controls only report proximity groups and cannot configure those calculations.

`PARAMS.WINDOWS` now exposes seventeen named distances. Every omitted option
retains its historical value. The resolver rejects unknown names, negative
values, strings, floats and booleans before constructing the DAG. The full
parameter table and boundary conventions are in `docs/window-parameters.md`;
tracked example configs declare the default values explicitly.

The OUTSIDER frequency region builder and evidence counter share one window.
INSIDER interval preparation and depth lookup share one depth margin. Clustering
within sources and matching between sources have independent BEDTools distances.
Sniffles read extraction exposes matching distance, fetch radius and sequence
flanks. Its BAM fetch interval is clipped to chromosome bounds so larger
configured windows cannot request negative reference coordinates. Zero is allowed
and retains each command's strict or inclusive comparison semantics.

Analysis windows remain separate from report grouping, VCF sequence-length
filters, TSD motif/flank settings, Liftoff and existing resident-TE chain settings.
Changing them can intentionally change retained events and frequencies. Different
configurations should use separate output directories, particularly with the
bundled Snakemake 5, which does not automatically rerun every parameter change.

## Validation

152 workflow tests, 8 module tests and 6 locus tests pass. Synthetic indexed BAMs
verify that wider Sniffles fetch/match windows recover an otherwise absent read,
that sequence flank lengths control extracted fragments and that the frequency
window changes breakpoint evidence from E to I. Other tests cover default
compatibility, independent overrides, zero, invalid values, typo rejection and
the INSIDER-deletion exclusion distance.

Two full 71-job runs (default and seventeen custom values) and a 43-job Sniffles2
run with clipping disabled complete. A printed-command dry run verifies every
custom value reaches the expected command; a negative value fails before execution.
With defaults, final call/candidate tables, frequencies, TSDs, integration audit,
Liftoff audit and both rebuilt genomes are byte-identical to the preceding run.
Default runs retain 529 calls and 23 integrated insertions; custom windows retain
529 calls with a changed source/type distribution and 22 integrated insertions.
Report projections and CIGAR/genome integration checks pass in all applicable
runs. A subsequent dry run reports that no work remains.
