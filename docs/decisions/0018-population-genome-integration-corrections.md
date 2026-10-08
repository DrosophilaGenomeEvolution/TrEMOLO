# 0018 — Correct population genome integration

Status: accepted, supersedes the sequence and coordinate compatibility choices
in 0015.

## Population representation

The user clarified that the input reads represent a population of individuals
with divergent TE content. The reconstructed genome deliberately combines all
retained insertions, including alternatives at the same breakpoint. Their stable
input order remains the concatenation order. This output is a population
representation, not an inferred individual haplotype. No allele selection or
same-site deduplication is introduced.

## Insertion sequences and DEL events

Classified candidate FASTAs contain read flanks for BLAST and TSD analysis. New
companion FASTAs, with identical identifiers, contain only CIGAR insertion
sequences or the caller's consensus sequence. The direct extraction uses query
positions and the configured flank size to handle truncated flanks accurately.
The Sniffles read extractor writes the exact CIGAR I slice, including correct
consumption of X and N operations. Integration reads these companion FASTAs.
Classifications continue to use the original candidate FASTAs.

Only INS events enter either reconstructed genome. DEL events are recorded as
`excluded / non_insertion_event`; their deleted sequence is already present in
the starting assembly and must not be inserted again. This branch does not
apply population deletions to the starting assembly.

## Liftoff coordinates

An inserted BED interval `[start,end)` yields GFF flanks
`[max(0,start-flank)+1,start]` and
`[end+1,min(chromosome_length,end+flank)]`. Empty flanks at contig boundaries are
omitted instead of generating invalid GFF. Missing-flank projection rejection
remains available when the other flank maps.

Mapped pairs must share chromosome, TE family, strand and compatible order.
Forward pairs use the inner BED boundaries `left.end` and `right.start-1`;
reverse pairs use `right.end` and `left.start-1`. The existing maximum gap applies
to the signed difference of these boundaries. Small ordered overlaps are retained
because target-site duplication can cause overlap. Opposite strands, reversed
order, invalid intervals and unknown strands are rejected. Audit columns include
the flank and projected strands. Public BED schemas are unchanged.

## Completeness

An eligible INS rejected during sequence resolution now fails the step before
publishing genome outputs. Audit and counts are written before failure; Snakemake
may remove declared audit outputs of a failed job, while rejection details remain
in its log. An explicit `CHOICE.OUTSIDER_VARIANT.INTEGRATION_ALLOW_PARTIAL: true`
adds `--allow-partial`, retaining the previous partial-output ability. Report
schema 1.4.0 adds an integration summary, counts and rejection reasons. Expected
insertions exclude intentionally excluded DEL events.

## Scope and compatibility

Historical NEO FASTA and BED bytes intentionally change. The previous comparator
in `check_variant_calling.py` remains a frozen migration oracle; use
`check_outsider_integration.py` for the corrected reconstruction. The current
candidate provenance and report checkers remain applicable.

TSD-based breakpoint choice, arbitration between candidate sequences, genome
memory usage remain separate work. Auditing events with both Liftoff flanks
entirely unmapped is subsequently addressed by decision 0019.
