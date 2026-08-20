# 0010 — Replace ALL_TE with a resident-TE annotation model

## Status

Accepted as the first scientific redesign after the compatibility migration.

## Context

The historical `DETECT_ALL_TE` rule annotated the complete assembled genome by
BLASTing the genome against the TE library twice. Between the two searches it
merged every first-pass interval separated by at most 100 bases. It then chose
one best family for each merged region and required 80 percent consensus
coverage for the combined call.

That contract is inappropriate for resident elements. Old, fragmented and
interrupted copies remain biologically useful even when they cover much less
than 80 percent of a modern consensus. The fixed 100 bp merge can also collapse
independent copies, while a single best-family column hides homologous family
matches and composite or nested structures.

The old report attempted to preserve some of that information in
`MULTIPLE_TE_BY_ID.txt`. That file counted insertion/family pairs but discarded
alignment coordinates, so it could not distinguish competing names for one
sequence segment from several TE components. Its AWK implementation also
failed to flush the final insertion group.

## Decision

Resident annotation is now an independent `TE_GENOME` pipeline rather than an
INSIDER sub-option. `CHOICE.PIPELINE.TE_GENOME` enables it, while the historical
`CHOICE.INSIDER_VARIANT.DETECT_ALL_TE` remains a temporary alias. The module can
operate without INSIDER or OUTSIDER. `PARAMS.TE_GENOME.TARGETS` defaults to
`["GENOME"]`; `REFERENCE` can be added when a prepared reference is available.

The search direction is reversed and only one BLAST is performed:

```text
TE consensuses (query) -> complete target assembly (BLAST database)
```

The explicit 14-column output includes consensus and target lengths. No
candidate-region extraction or second BLAST is performed.

Discovery and classification are separate:

- `MIN_PIDENT`, `MIN_ALIGNED_BP` and `MAX_EVALUE` define significant fragments;
- `FULL_LENGTH_COVERAGE`, `HIGH_CONFIDENCE_PIDENT` and `PARTIAL_COVERAGE` assign
  `full_length`, `partial` or `degraded_relic` labels;
- consensus coverage is the union of aligned consensus intervals, not the
  genomic span between the first and last HSP;
- the full-length threshold never removes partial or degraded annotations.

Accepted HSPs are chained only within the same target sequence, TE family and
strand, and only when genomic and consensus coordinates are collinear. The gap
tolerance combines a base-pair minimum with a fraction of consensus length.
This replaces the family-independent 100 bp region merge.

Family matches with high reciprocal genomic overlap are treated as alternative
annotations of one physical component. Complete-linkage grouping is required:
every member must overlap every other member sufficiently, preventing
single-linkage chains from merging neighbouring copies. One primary match is
retained for compatibility and display, but every alternative remains in the
normalized tables and GFF3.

Contained or partially overlapping components that are not alternative family
matches remain separate copies. Their relations are recorded as
`nested_candidate`, `overlap_candidate` or `same_family_overlap`, always with
`provisional` confidence. Geometry alone does not prove a nested insertion.

## Outputs

Each target is written below `TE_GENOME/<TARGET>/`:

| file | contract |
|---|---|
| `ALL_TE_FRAGMENTS.tsv` | every BLAST HSP, including rejected rows and rejection reasons |
| `ALL_TE_MATCHES.tsv` | collinear family-specific chains and primary/alternative assignment |
| `ALL_TE_COPIES.tsv` | physical components with all candidate families and completeness tier |
| `ALL_TE_RELATIONS.tsv` | provisional containment and overlap relations between components |
| `POSITION_ALL_TE.bed` | one standard BED6 row per physical component |
| `ALL_TE.gff3` | parent components and every family match for genome-browser display |

The public `POSITION_ALL_TE.bed` symlink continues to point to the GENOME BED
when GENOME is selected. Stable opaque identifiers are hashes of normalized
coordinates and evidence rather than row numbers, so adding an unrelated later
record does not renumber existing annotations.

## Initial fixture result

With the default resident thresholds on the current Drosophila fixture:

- 11,589 BLAST HSPs are retained in the audit table;
- 10,746 pass discovery filters;
- 8,336 family-specific matches form 6,224 physical components;
- 364 components are `full_length`, close to the 327 strict legacy calls;
- 320 of 327 legacy intervals overlap a new full-length component, and 338 of
  364 new full-length components overlap a legacy interval;
- all 327 legacy intervals overlap at least one component in the complete new
  catalogue, so the permissive model adds evidence without dropping the old
  strict regions;
- 1,818 are `partial` and 4,042 are `degraded_relic`;
- 1,335 components retain more than one family match;
- 7,301 geometric relations are explicitly provisional.

These counts are a regression fixture, not evidence that the default thresholds
are universally optimal.

## Consequences and limits

The new outputs are intentionally not byte-compatible with legacy ALL_TE. The
variant-associated INSIDER/OUTSIDER calls, `TE_INFOS.bed`, their frequencies,
TSD values and the Quarto report remain unchanged.

The current Python implementation holds parsed HSPs in memory and is therefore
expected to need further streaming or indexed grouping work for very large
repeat-rich genomes. Primary-family ranking is a deterministic display
convenience, not a claim that alternatives are biologically false. Nested and
composite structures must later be validated with sequence continuity or read
evidence before receiving a resolved status.

The same fragment/match/component distinction is intended for subsequent
INSIDER and OUTSIDER redesigns, but this decision changes only TE_GENOME.
