# 0015 — Make OUTSIDER integration valid and audit Liftoff projections

Status: accepted during the Snakemake 5.10 migration.

## Context

The historical `TE_TOWARD_GENOME` rule generated two reconstructions:

- `NEO_GENOME.fasta` inserts the sequence observed in the read/SV call;
- `PSEUDO_GENOME_TE_DB_ID.fasta` was intended to insert the oriented canonical
  TE sequence.

The test oracle exposes two defects in the canonical reconstruction. Every
inserted sequence was prefixed with its event identifier, including digits and
punctuation, so the result was not a DNA FASTA. In addition, a malformed shell
pipe in the reverse-complement branch reduced seven negative-strand insertions
to the identifier alone (15 or 16 non-nucleotide characters) instead of their
TE sequence. Preserving these bytes would retain an invalid biological object
and give Liftoff a corrupted source assembly.

The old `LIFT_OFF` rule also hid partial mappings and discordant chromosomes in
shell loops. A left and right flank mapped to different chromosomes could be
promoted to the regular output, and one branch mixed a coordinate from one
chromosome with the name of the other.

## Decision

`INTEGRATE_TE_TO_GENOME` is now an explicit optional DAG branch. It integrates
all non-`SOFT`/non-`HARD` merged calls in coordinate order and produces:

- a canonical reconstruction in the legacy
  `PSEUDO_GENOME_TE_DB_ID.fasta` path, containing nucleotide sequence only;
- an observed reconstruction in `NEO_GENOME.fasta`;
- shifted BED coordinates for both reconstructions;
- `INTEGRATION_TE.tsv`, with source, strand, lengths, status, and rejection
  reason for each merged call.

On the bundled oracle, 25 calls are integrated. The observed reconstruction
and its BED remain byte-identical to the legacy result. The canonical result is
an intentional scientific correction: all 25 canonical sequences are present,
negative-strand sequences are reverse-complemented, and event identifiers are
kept only in BED/audit metadata.

When both INSIDER and OUTSIDER are enabled, Liftoff projects two 100 kb flanks
per integrated insertion from the corrected canonical reconstruction to the
reference. `LIFT_OFF_AUDIT.tsv` distinguishes concordant projections, excessive
gaps, missing/ambiguous flanks, and discordant chromosomes. Discordant pairs
remain in `BAD_POS_TE_LIFT.bed` and are never promoted to
`POS_TE_OUTSIDER_ON_REF.bed`.

The defaults retain the historical flank and gap thresholds:

```yaml
PARAMS:
  OUTSIDER_VARIANT:
    LIFT_OFF:
      FLANK_SIZE: 100000
      MAX_GAP: 20000
```

## Consequences

The legacy filename ending in `_DB_ID` is retained for downstream path
compatibility even though identifiers are no longer embedded in DNA. Exact
comparison of that canonical FASTA to `work_test` is neither expected nor
desirable. The observed genome is still an exact migration oracle.

Liftoff output can vary with assembly repetitiveness and tool version; the
public call set is therefore validated through explicit projection states and
coordinate constraints rather than an assumed byte-identical GFF3. The image
currently provides Liftoff 1.6.3.
