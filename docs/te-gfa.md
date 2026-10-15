# TE insertion graph

`lib/python/locus/build_te_gfa.py` builds a reference-backed GFA 1.1 from
reported insertion events and their raw insertion/read evidence. It does not
construct consensus sequences. This first version supports OUTSIDER direct
and Sniffles insertion channels; other events are recorded as excluded.
It is an explicit post-processing command, not yet an automatic Snakemake target.

```bash
python3 lib/python/locus/build_te_gfa.py \
  --reference /path/reference.fasta \
  --te-database /path/original_TE_database.fasta \
  --reference-name reference \
  --run pool=/path/analysis \
  --output /path/analysis/GFA
```

Repeat `--run SAMPLE=DIR` to combine runs with the same reference. The command
checks each run's `INPUT/genome.fasta` against the supplied reference sequences.
Use the original TE database, retaining its complete FASTA headers, as used
when preparing the run. `INPUT/te_headers.tsv` restores the original TE labels
from the internal BLAST identifiers. Old analyses missing `SEQUENCE_PROVENANCE.tsv`
must rerun insertion extraction before this command can use their channels.

To select TE names, add either `--te-names 'roo,copia,I-element'` or
`--te-names-file /path/TE_names.txt` (one exact original TE name per line).
Without either option, all loci remain eligible. A locus is retained if its
reported primary label or any qualifying alternative matches a selected name.
All supported family representatives at that locus remain in the graph,
including those absent from the selection. Existing ambiguity and alignment
filters still apply. Other loci leave the reference sequence unchanged; their
sequence evidence is preserved with `excluded_by_te_filter` where applicable.
The tables also record `locus_selected`/`selected`, and the manifest records
the name selection and selected locus count. An empty list is rejected.

To keep only the main dm6 chromosomes, add:

```bash
--chrom-regex '^([23][RL]|[X4Y])$'
```

This is a Python regular expression searched on the FASTA sequence name alone,
without its description or following tabular columns. Use `^` and `$` for exact
names. Unlike `CHROM_KEEP`, commas are not converted into alternatives: use `|`.
Without the option, every reference chromosome/contig is retained, including those
without TE calls. With it, excluded chromosomes contribute no segments, links,
walks, loci or sequence candidates. Their BED events are recorded in
`excluded-events.tsv` with reason `excluded_by_chrom_filter`. Chromosome and TE
name filters can be combined. The manifest records the regex and retained sequence
names. Invalid patterns or patterns matching no reference sequence fail before
output files are written. Input-reference consistency is still checked against
the complete reference before applying this export filter.

For explicit input files, use `--inputs inputs.json` instead of `--run`:

```json
[
  {
    "sample_id": "pool",
    "te_infos": "/path/analysis/TE_INFOS.bed",
    "reference": "/path/analysis/INPUT/genome.fasta",
    "te_headers": "/path/analysis/INPUT/te_headers.tsv",
    "channels": [
      {
        "name": "outsider_direct_alignment",
        "alignment": "/path/analysis/OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln",
        "insertion_fasta": "/path/analysis/OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.insertions.fasta",
        "provenance": "/path/analysis/OUTSIDER/TrEMOLO_SV_TE/INS/SEQUENCE_PROVENANCE.tsv"
      }
    ]
  }
]
```

Paths in JSON are absolute or relative to the working directory.
`reference` and `te_headers` are optional for manually prepared inputs; without
`reference`, the caller is responsible for matching coordinates to the supplied
reference. Without `te_headers`, BLAST names must already match the original DB.

## Selection and ranking

Only insertion events present in `TE_INFOS.bed` are considered. Raw alignments
are loaded for all their sequences, including events with a single candidate
family. The ambiguous-call evidence TSV alone would omit such events.

For each query and family, dynamic programming chooses the highest-bitscore
collinear HSP chain with one orientation and no overlap on either query or TE.
Scores are summed only within that chain, not across reads. Identity is weighted
by aligned HSP length. Coverage uses the aligned consensus intervals divided
by the database sequence length. HSP identity must pass `--min-identity`;
the chain must pass `--min-coverage` and `--min-aligned-bp` (defaults: 80%, 80%,
100 bp). These are explicit graph thresholds; they are not parsed from run YAML.

A sequence with qualifying matches to multiple families is marked
`ambiguous_or_composite_annotation` and excluded from graph alternatives.
The check also covers identical insertion sequences in multiple reads or
channels at the same locus. This initial conservative policy leaves composite
structures in the table rather than inventing several independent alternatives.

Within a locus, sequences are grouped by original TE label. One representative
is selected by decreasing chain bitscore, coverage, identity, and number of
distinct reads supporting that exact insertion sequence within that family.
Sample/channel/query/family identifiers break remaining ties. Representatives
are then ordered by the same rule to determine levels. Bitscores across families
can favor longer TEs. Levels describe this ranking, not biological frequency.

The full extracted insertion is used, including any non-TE bases within it;
read flanks are absent from the insertion FASTA. A reverse representative is
reverse-complemented for its segment and traversed with `<ID`, reconstructing
the observed genomic orientation. A VCF-reported sequence remains explicitly
`vcf_reported` with no assigned read; the exporter does not create a consensus
but does not change how Sniffles produced that input sequence either.

## Loci and walks

The exporter reuses complete-span breakpoint clustering from the cohort model.
The default `--locus-window 0` groups exact anchors only; a nonzero window makes
candidate loci and uses their minimum anchor for graph placement. Such loci
are marked `provisional_anchor`; proximity is not proof of a shared biological
locus. TE family and orientation do not determine locus grouping. Conflicting
orientations within one family are simplified to the selected representative
and remain visible in the sequence table.

The reference walk reconstructs every supplied chromosome exactly. At level 1,
each represented locus uses its highest-scoring family. At level 2 it uses its
second family, and so on. After the last alternative at a locus, higher levels
retain that last choice. Levels are generated per chromosome; chromosomes with
no insertion alternatives have only a reference walk. No walk concatenates two
alternatives at the same locus. The reference is a separate walk, not a claim
that reference absence was observed in this population.

These walks are synthetic and unphased. Comparing family alternatives between
the same anchors exposes replacements (SWAPs); comparison to the reference
exposes insertions. There is no special SWAP record in GFA.

## Outputs

- `te-graph.gfa`: header `H VN:Z:1.1 RS:Z:<reference>`, globally unique positive
  numeric segments, chromosome blocks with `S`, `W`, then oriented `L` records
  with `0M`. Walk end coordinates are reconstructed sequence lengths.
- `sequence-alternatives.tsv`: qualifying and below-threshold chain summaries,
  scores, provenance, representative status, node and level assignments.
  HSPs below the identity threshold do not enter these chain summaries.
- `loci.tsv`: selected anchors, uncertainty and alternative counts.
- `excluded-events.tsv`: unsupported events and insertions lacking usable evidence.
- `walks.tsv`: reference/synthetic interpretation and walk lengths.
- `manifest.json`: input checksums, thresholds, ranking and graph counts.

The record structure follows the local Savanache full-GFA requirements and the
[GFA specification](https://gfa-spec.github.io/GFA-spec/GFA1.html).
The local Savanache `check_format_gfa` accepts the two-chromosome test graph
without errors or warnings, and its text/binary indexes are generated.
The small synthetic [example GFA](te-gfa-example.gfa) has three alternatives
at the first locus, one at the second, and a reverse-oriented representative.
The local index wrapper leaves the first chromosome's SQLite position table
empty during this test, so complete application compatibility is not yet
established. Unit tests validate sequences and links, not visualization or the
Savanache SWAP caller.
This version does not encode INSIDER deletions, resident TE alleles, nested
structures, or observed haplotypes.
