# 0013 — Refactor optional modules around normalized scientific tables

## Status

Accepted.

## Context

The historical population module clustered `TE_INFOS.bed` rows with a shell
pipeline, treated missing calls as zero in the graph, wrote derived files into
input work directories and relied on CDN-hosted JavaScript. Its cluster state
also mixed the event-ID and bedtools-cluster columns.

The historical BLAST module declared scientific thresholds without applying
them. It serialized all alignments into one JSON file, lost remainders when
splitting that file, required a mutable repository symlink and rendered HTML
from a synchronous Node server.

## Decision

Both modules now write versioned, tabular scientific outputs before rendering
an interface. Their Quarto reports embed the normalized JSON projection and use
local JavaScript and CSS only.

Population analysis builds the cohort locus model with complete-span candidate
clustering. Every sample/allele combination receives `observed`,
`observed_unquantified` or `missing`; missing is never inferred as frequency
zero or confirmed absence. Reference identity is the input FASTA SHA-256.

Insertion structure analysis preserves every HSP, marks filter decisions,
chains matches per query/TE/orientation and proposes components only when their
aligned query segments satisfy the configured overlap limit. Query structures
are aggregated by TrEMOLO event so repeated sequence support is explicit.

## Consequences

The old shell-generated graphs, server-side HTML, Node server and example JSON
are removed. The historical entry-point scripts remain as wrappers around the
new builders. Reports no longer require network access or persistent services.

Locus assignments and insertion components remain provisional. Future
breakpoint validation and population empty-site genotyping can refine these
tables without changing the interface contract or rewriting legacy
`TE_INFOS.bed`.
