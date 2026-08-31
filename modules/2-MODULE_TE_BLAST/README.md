# Insertion structure explorer

This module normalizes the BLAST alignments between OUTSIDER insertion
sequences and the TE database. It is designed to expose competing TE matches,
possible composite insertions and possible repeated or nested structures
without silently converting them into confirmed biological components.

## Run

```bash
./module build 2 WORK_DIRECTORY [OUTPUT_DIRECTORY]
```

The default output is `WORK_DIRECTORY/MODULE_TE_BLAST`. The module consumes:

- `OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln`;
- `OUTSIDER/TrEMOLO_SV_TE/INS/SV_SIZE.tsv`;
- `1-UTILS/TE_SIZE.tsv`;
- `TE_INFOS.bed`, when available, to mark final reported events;
- `TE_CALL_CANDIDATES.tsv`, when available, to identify known alternatives.

Scientific thresholds can be passed after the output directory:

```bash
./module build 2 work_test STRUCTURE_REPORT \
  --min-pident 90 \
  --min-aligned-bp 80 \
  --min-consensus-coverage 20 \
  --max-component-overlap 0.2
```

All raw HSPs are preserved regardless of thresholds. Thresholds only determine
the normalized match status and the provisional component proposal.

## Outputs

```text
STRUCTURE_REPORT/
  insertion-hsps.tsv
  insertion-matches.tsv
  insertion-components.tsv
  insertion-structures.tsv
  insertion-events.tsv
  manifest.json
  report-data.json
  report.qmd
  report.html
  index.html
```

- `insertion-hsps.tsv` contains every BLAST HSP in 0-based half-open query and
  subject coordinates, with explicit filter reasons.
- `insertion-matches.tsv` chains HSPs by query, TE and orientation and records
  union coverage, weighted identity, score, final-call status and assignment.
- `insertion-components.tsv` contains a conservative, non-overlapping set of
  coordinate-supported matches for each query sequence.
- `insertion-structures.tsv` describes each supporting insertion sequence.
- `insertion-events.tsv` aggregates all supporting sequences belonging to the
  same TrEMOLO event. This is the preferred table for finding events supported
  by several query sequences.

Event classifications distinguish:

- `multi_te_candidate`: multiple TE families occur on one supporting sequence;
- `multi_te_supported`: the pattern occurs on at least two supporting
  sequences;
- `repeated_or_rearranged_te_candidate/supported`: several components have the
  same TE label;
- `single_candidate`, `unresolved`, and `no_match`.

These labels remain hypotheses. A composite or nested biological allele should
only be promoted after read/assembly evidence confirms linkage and breakpoints.

The report is a self-contained Quarto HTML file. The historical Node server,
global repository symlink, 500 MB monolithic JSON and CDN dependencies are no
longer used.
The refactored Singularity definition includes Quarto; with an older image,
run the wrapper from a host installation of Quarto instead.
