# Reference-guided assembly correction

TrEMOLO provides two standalone Bash utilities for making an assembly follow
the order, chromosome assignment and orientation supported by a reference:

- `correction_genome_transloc_sim_chrom.sh` reorders or reverses blocks inside
  their original chromosome;
- `correction_genome_transloc_diff_chrom.sh` can additionally transfer blocks
  between chromosomes.

Both scripts use the same filtering, planning, reconstruction and audit code.
They preserve every base of the query exactly once. Query contigs that are not
selected, do not match `--chrom_regex`, or have no accepted anchor are copied
unchanged. Existing FASTA indexes are never removed or overwritten.

## Requirements

- Bash;
- minimap2, unless a PAF file is supplied;
- samtools;
- bedtools;
- standard POSIX utilities (`awk`, `sort`, `cut`, `cmp`).

## Basic use

Intrachromosomal correction:

```bash
lib/bash/correction_genome_transloc_sim_chrom.sh \
  -c chromosomes.txt \
  --size_min 20000 \
  --min-mapq 20 \
  --min-identity 80 \
  reference.fasta assembly.fasta corrected.fasta
```

Interchromosomal correction:

```bash
lib/bash/correction_genome_transloc_diff_chrom.sh \
  -c chromosomes.txt \
  --size_min 20000 \
  --min-mapq 20 \
  --min-identity 80 \
  reference.fasta assembly.fasta corrected.fasta
```

An existing alignment can be passed as a fourth positional argument. Its PAF
orientation must be `reference` as target and `assembly` as query:

```bash
lib/bash/correction_genome_transloc_sim_chrom.sh \
  [options] reference.fasta assembly.fasta corrected.fasta mapping.paf
```

The chromosome correspondence file contains one one-to-one pair per line. A
colon or a tab can be used, and either orientation is accepted when names are
unambiguous:

```text
query_chr2L:reference_chr2L
query_chr2R:reference_chr2R
```

Without `--chrom`, a deterministic one-to-one correspondence is inferred from
the cumulative length of accepted PAF alignments. Supplying the file is
recommended when chromosome names differ or repetitive mappings make the
correspondence ambiguous.

## Anchor filters

`--size_min` (also named `--min-anchor-size`) is applied before correction to
the query span, reference span and PAF alignment-block length. It therefore
controls the minimum evidence used to define a corrected block; its default is
20 kb in both scripts.

Additional filters are:

- `--min-mapq`, default `20`;
- `--min-identity`, default `80` percent, computed from PAF matches/block size;
- `--max-anchor-overlap`, default `1000` bp;
- `--chrom_regex`, default `.`, applied to query chromosome names.

Secondary PAF alignments (`tp:A:S`) are rejected. Remaining candidates are
selected deterministically, prioritizing anchors supporting the declared
chromosome pair, then alignment length, MAPQ and identity. Query and reference
overlaps above the configured limit are excluded.

## Planning and audit files

Use `--report-only` to inspect the planned changes without writing a corrected
FASTA. By default, audit files use `<output>.correction` as their prefix; this
can be changed with `--audit-prefix`.

| File suffix | Content |
| --- | --- |
| `.anchors.tsv` | Every PAF row and its accepted/rejected reason |
| `.selected_anchors.tsv` | Deterministic non-overlapping anchor selection |
| `.chromosomes.tsv` | Query/reference chromosome pairs and alignment score |
| `.corrections.tsv` | Source interval, destination, orientation, order and action |
| `.minimap2.log` | minimap2 stderr when the alignment is generated |

Possible correction actions are `KEEP`, `MOVE`, `REVERSE`,
`MOVE_AND_REVERSE`, `TRANSFER`, and `TRANSFER_AND_REVERSE`.

Before completing, each script checks that source intervals are gap-free and
non-overlapping, every input base is represented once, total assembly length is
unchanged, and chromosome names remain present in their original order. The
intrachromosomal script additionally requires every chromosome length to stay
unchanged.
