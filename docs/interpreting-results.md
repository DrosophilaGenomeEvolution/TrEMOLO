# Interpreting results

## Calls, absence and frequency

`TE_INFOS.bed` contains reported variant-associated events. Several calls can
occur at one locus; row counts are not counts of independent biological alleles.
An absent call does not establish TE absence. Low coverage, thresholds and
unsupported evidence can all prevent a call.

Module 1 distinguishes `observed`, `observed_unquantified` and `missing`.
Missing values are not converted to zero. Its trajectories use replicate means
at observed timepoints. Every input must use the same reference coordinate
system; supplying one reference FASTA assigns a checksum but does not prove
that all input BED files actually use that reference.

Read-derived frequency estimates depend on the detection channel and support
rules. See [frequency columns and INS/DEL semantics](frequency-tables.md).
Candidate evidence fractions are shares of retained evidence, not biological
allele frequencies or necessarily proportions of unique reads.

## Families and insertion structures

The BED's primary TE family is the reported annotation, not a guarantee that
the insertion contains only that family. Candidate tables retain alternative
annotations and the sequence/read evidence supporting them. Several families
on one query differ from families supported on separate insertion sequences.

Module 2 describes provisional components and composite/nested/repeated
structures. Its `supported` labels summarize sequence evidence; they do not
replace independent confirmation of linkage, breakpoints or read independence.
Rejected HSPs remain inspectable. A family filter can therefore return a card
whose provisional components do not include that family.

## Graphs and reconstructed genomes

GFA levels rank selected family representatives by alignment evidence. They do
not encode frequency or observed phased haplotypes. A level-1 family can differ
from the BED primary label. Graph thresholds are independent of pipeline YAML
thresholds. Review `summary.html`, `call-representation.tsv` and exclusions;
a graph may retain reference walks while omitting unsupported events.

The exporter handles supported OUTSIDER insertion channels, not INSIDER deletion
alleles, resident copies, nested allele construction or observed haplotypes.
See [GFA rules](te-gfa.md).

Integrated population genomes intentionally concatenate same-site alternatives.
Inspect `INTEGRATION_TE.tsv` for skipped or rejected events. With partial
integration enabled, a successful job does not mean every eligible insertion
was integrated. Rejected Liftoff projections do not remove insertions from the
reconstructed genomes.

## Coordinates

- Variant calls are located on the input mapping/comparison assembly unless
  the output explicitly identifies reference-projected coordinates.
- BED intervals use zero-based half-open coordinates. Insertion anchors are
  represented as intervals; they are not the span of inserted sequence.
- Liftoff's intermediate GFF uses one-based inclusive coordinates.
- Shifted BED outputs refer to their corresponding reconstructed genome.
- Module 2 HSP coordinates refer to insertion queries and TE database sequences.
- Unmapped-read evidence is positioned within reads, not genome chromosomes.

Do not combine these coordinate systems without an explicit conversion.
Nearby calls and nonzero locus windows define provisional groups; proximity
alone does not establish a shared biological locus.
