# TrEMOLO locus data model

Status: draft v0.1. This model complements the legacy `TE_INFOS.bed`; it does
not replace or modify that file.

## Why a separate model

`TE_INFOS.bed` stores one detected TE event per row. A population analysis
needs to distinguish concepts that cannot safely be represented by adding more
columns to that row:

- a genomic locus shared between samples;
- alternative insertion alleles at that locus;
- one or several TE components within an allele;
- the observation and frequency of each allele in each sample;
- the original calls supporting the interpretation.

In particular, a TE family is **not** part of locus identity. Two different TE
families can be alternative alleles of the same locus. Conversely, two nearby
insertions are not necessarily the same locus.

## Output layout

The normalized files live beside `TE_INFOS.bed`:

```text
LOCUS_MODEL/
  manifest.json
  loci.tsv
  alleles.tsv
  components.tsv
  observations.tsv
  evidence.tsv
```

All TSV files are UTF-8, tab-separated, include one header row, use 0-based
half-open genomic coordinates, and use `.` for an unknown value. IDs are opaque
and must not be parsed to recover biological attributes.

### `manifest.json`

Records `schema_version`, TrEMOLO version and commit, reference identity,
coordinate system, sample metadata, creation command, and checksums of input
files. Combining runs is forbidden when their reference identity differs unless
an explicit liftover has been performed.

### `loci.tsv`

| column | meaning |
|---|---|
| `locus_id` | stable identifier independent of TE family and sample |
| `chrom` | reference sequence |
| `start` | left boundary of the locus uncertainty interval |
| `end` | right boundary of the locus uncertainty interval |
| `anchor_left` | consensus left breakpoint, if known |
| `anchor_right` | consensus right breakpoint, if known |
| `status` | `resolved`, `ambiguous`, or `provisional` |
| `confidence` | normalized confidence in `[0,1]`, or `.` |
| `method` | locus construction method and version |

One row represents one genomic site, never one TE family.

### `alleles.tsv`

| column | meaning |
|---|---|
| `allele_id` | stable allele identifier |
| `locus_id` | parent locus |
| `allele_type` | `reference`, `insertion`, `deletion`, or `complex` |
| `structure` | `single`, `nested`, `tandem`, `composite`, `absence`, or `unknown` |
| `length` | allele length in bp, or `.` |
| `sequence_id` | identifier in an optional allele FASTA, or `.` |
| `confidence` | normalized confidence in `[0,1]`, or `.` |

A locus can have zero or one reference allele and any number of non-reference
alleles. Two alleles at the same locus may contain different TE families.

### `components.tsv`

| column | meaning |
|---|---|
| `component_id` | component identifier |
| `allele_id` | parent allele |
| `rank` | component order on the allele, starting at 1 |
| `parent_component_id` | containing component for a nested TE, otherwise `.` |
| `te_name` | matched TE name |
| `te_family` | normalized family, or `.` |
| `te_class` | normalized class/superfamily, or `.` |
| `start_on_allele` | 0-based component start on the allele, or `.` |
| `end_on_allele` | half-open component end on the allele, or `.` |
| `strand` | `+`, `-`, or `.` |
| `identity` | alignment identity in percent, or `.` |
| `coverage` | TE consensus coverage in percent, or `.` |

Several rows with one `allele_id` explicitly represent several TEs in the same
allele. `parent_component_id` represents nesting rather than mere proximity.

### `observations.tsv`

| column | meaning |
|---|---|
| `sample_id` | sample or population identifier |
| `locus_id` | observed locus |
| `allele_id` | observed allele |
| `frequency` | estimated allele frequency in `[0,1]`, or `.` |
| `support_reads` | supporting read count, or `.` |
| `total_reads` | informative read count, or `.` |
| `genotype` | optional genotype such as `0/1`, otherwise `.` |
| `quality` | observation quality, or `.` |
| `filter` | `PASS` or a semicolon-separated reason list |

For a pooled population, `frequency` is primary and `genotype` is `.`. For an
individual, both may be provided. Frequencies at one locus are allele-specific;
absence must not be inferred merely because a caller emitted no row.

### `evidence.tsv`

| column | meaning |
|---|---|
| `call_id` | original TrEMOLO call identifier |
| `sample_id` | source sample |
| `locus_id` | assigned locus |
| `allele_id` | assigned allele, or `.` when ambiguous |
| `component_id` | supported component, or `.` |
| `source` | `INSIDER` or `OUTSIDER` |
| `source_file` | path relative to the run directory |
| `assignment` | `exact`, `compatible`, `ambiguous`, or `rejected` |
| `reason` | machine-readable assignment explanation |

This table preserves traceability and permits locus construction to be rerun
without losing the original calls.

## Locus construction rules

The first implementation must use candidate generation followed by compatibility
testing. A fixed distance alone is insufficient.

1. Generate candidate pairs on the same chromosome using breakpoint uncertainty.
2. Compare left and right breakpoints when both exist.
3. Compare TSD, insertion length, orientation and insertion sequence when known.
4. Treat TE family as allele/component evidence, never as locus identity.
5. Keep incompatible nearby calls in separate loci.
6. Permit compatible calls with different TE families in the same locus but in
   distinct alleles.
7. Assign a composite or nested allele only with sequence/read evidence joining
   its components. Co-localization alone is not sufficient.
8. Record uncertain assignments as `ambiguous`; never resolve them silently.

The legacy module's `bedtools cluster -d 200` can be used only to generate
candidates. It cannot make the final locus assignment.

The initial cohort implementation uses complete-span clustering on normalized
breakpoint anchors. All anchors in a cluster must fit inside `--locus-window`;
single-linkage chaining is forbidden. Such loci remain `provisional` or
`ambiguous`. This is candidate construction, not final biological breakpoint
resolution.

## Required invariants

- Every allele references exactly one existing locus.
- Every component references exactly one existing allele.
- Component ranks are unique within an allele.
- A `parent_component_id` belongs to the same allele and cannot form a cycle.
- Every observation references an existing `(locus_id, allele_id)` pair.
- `0 <= frequency <= 1` and `support_reads <= total_reads` when both are known.
- Each original accepted call has at least one evidence row.
- IDs remain stable when unrelated samples are added to a cohort.

## Illustrative multi-allelic locus

```text
locus L0001
  allele A0: reference/absence
  allele A1: roo
  allele A2: copia
  allele A3: roo containing a nested copia
```

`A1` and `A2` each have one component. `A3` has two components, with the copia
component pointing to the roo component through `parent_component_id`. Each
sample or pool receives separate observation rows for the alleles it supports.

## Compatibility policy

During migration, `TE_INFOS.bed` remains the public legacy output. The locus
tables are derived outputs under an explicit schema version. Reports should
prefer the locus model when present and fall back to `TE_INFOS.bed` otherwise.
