# INSIDER backend options — 2026-09-11

Recommendation: retain Assemblytics as the compatibility backend, then evaluate
SyRI as the next backend on chromosome-level assemblies. This follows ADR 0001;
none of the alternatives below is implemented in TrEMOLO by this change.
The age of a caller alone does not establish its error rate on TE insertions.

| Candidate | Fit for INSIDER | Qualification |
| --- | --- | --- |
| SyRI | Chromosome-level assembly comparisons with synteny, inversions and rearrangement context | Requires paired chromosome orientation and appropriate assembly quality; large rearrangements can increase resource use |
| SVIM-asm | A more direct INS/DEL caller from haploid or diploid assembly alignments | Useful comparison backend, including contig alignments; does not replace the same global synteny interpretation |
| PAV | Broader analysis of phased assemblies and complex structural variants | A larger integration effort; prioritize when phased haplotypes and complex variants become primary requirements |

SyRI's repository documents chromosome-level inputs, same-strand homologous
chromosomes, and recent Python/pandas modernization. The fit for TrEMOLO above
is an engineering recommendation based on those requirements, not a measured
accuracy ranking. [SyRI documentation](https://github.com/schneebergerlab/syri)

SVIM-asm accepts haploid or diploid genome-genome alignments and calls insertions,
deletions and other SV classes. [SVIM-asm documentation](https://github.com/eldariont/svim-asm)

PAV provides assembly-based variant discovery for phased genomes.
[PAV documentation](https://github.com/EichlerLab/pav)

## Proposed integration sequence

1. Define a shared SV table containing caller/version, original event ID,
   reference and query coordinates, orientation, event type, sequence and
   confidence. Feed TE detection from this table rather than caller-specific
   Assemblytics artifacts.
2. Add experimental `SV_CALLER: syri`, keeping `assemblytics` as default.
   Check chromosome identity/orientation before running. Do not automatically
   rearrange assembly blocks to match the reference: that can erase the
   biological rearrangements being evaluated. Any orientation transformation
   must retain a coordinate mapping to the original assembly.
3. Compare TE recall, sequence recovery, breakpoint agreement, family assignment,
   tandem-repeat confusion and ambiguous mappings on curated and synthetic loci.
   Include long TEs, nested copies, inversions, absent events and fragmented input.
4. Add SVIM-asm if contig-level assemblies or an independent local-SV comparison
   are needed. Annotate caller agreement before considering a consensus mode.

A union of callers can inflate false positives; an intersection can lose real
TEs. Neither should become the default without this evaluation. INSIDER
empty-site frequency and TSD calculations also need validation against the new
coordinate contract before the default backend changes.
