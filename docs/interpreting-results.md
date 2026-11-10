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

Use the report's **Frequency versus local depth** plot to distinguish a high
frequency based on few reads from one supported by a larger denominator. The
adjustable guides do not establish statistical confidence. OUTSIDER uses an
adjusted read denominator; INSIDER uses mean flanking depth. These are local
support measures, not genome-wide enrichment. Missing and zero-depth estimates
are explicitly omitted from that plot.

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


## Resident consensus coverage and reporting tiers

Consensus coverage is `100 × union of aligned consensus intervals / database TE
sequence length`; overlapping intervals count once. For a 5,000 bp consensus,
1,000 covered bp means 20% coverage, regardless of read depth. Identity measures
matching bases in the alignments and is a separate quantity.

The interactive preview first assigns full-length when both its coverage and
identity thresholds are met. Otherwise, coverage at or above the partial
threshold assigns partial; remaining copies are degraded/relic. Thus a
high-coverage but lower-identity copy can be partial too. These are reporting
tiers, not experimentally confirmed activity or structural completeness.
The discovery filters constrain the evidence already available; the preview
sliders do not rerun discovery. Cards show all reported copies, independently
of resident-table filters.

The chromosome overview uses a shared base-pair scale, with actual reference
lengths beside each bar and exact bp values on hover. Its scale stays fixed
when filtering calls. Unknown lengths are labelled explicitly; dashed bars
show observed spans rather than claiming a chromosome length.


## Contextual report help and exploration

`REPORT/report-help.html` explains exact labels including `accepted`,
`nested_candidate`, `same_family_overlap` and `overlap_candidate`. Contextual links
in the report open the matching section. Share this companion beside the main HTML.

The detector summary counts final `TE_INFOS.bed` rows by detector and event type;
Sniffles VCF totals remain separate because they describe raw SVs before TE
annotation. HARD/SOFT channel calls are counted separately from INS/DEL columns.

The additional Magma heatmap counts call anchors in adjustable genomic windows,
using all reported calls for a fixed colour maximum. Blue outlines identify
windows losing calls to the active call filters. Tooltips show selected and
baseline counts by INSIDER/OUTSIDER source. Empty windows do not establish absence
or callability. Terminal insertions at chromosome length are assigned to the
last window. A display safeguard increases window width if necessary; the actual
width is always reported. The individual-call chromosome plot is unchanged.

Structure preview sliders apply to both endpoints of already reported relations,
using primary-match identity and consensus coverage. At 0%/0% there is no extra
threshold. They do not rerun discovery, select alternative primary labels or
regroup copies. The fixed global relation summary remains visible, and this
preview is independent of the call and resident-table filters.
