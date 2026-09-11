# 0016 — Select Sniffles generations explicitly

Status: accepted, 2026-09-11.

## Selection and installation

`CHOICE.OUTSIDER_VARIANT.CALL_SV` accepts `sniffles1`, `sniffles2`, and the
backwards-compatible `sniffles` alias for generation 1. Invalid choices fail
at DAG construction instead of silently running another caller. This applies
to `workflow/Snakefile`; the archived `run.snk` entry point is unchanged.

The image retains 1.0.12b at `sniffles` (also `sniffles1`) and installs 2.8.1
at `sniffles2`. Generation 2 runs in `/opt/sniffles2` with Python 3.12 because
its Python >=3.10 and NumPy >=2.2 requirements conflict with the legacy stack.
The environment is not prepended to PATH; a dedicated executable symlink selects
it without replacing Python, NumPy, or the original Sniffles binary.

`TOOLS.SNIFFLES1` and `TOOLS.SNIFFLES2` override executable paths. The runner
checks the actual major version before calling and records the version and
argument list in `OUTSIDER/VARIANT_CALLING/sniffles-run.json`.
`PARAMS.OUTSIDER_VARIANT.SNIFFLES.MIN_SUPPORT` is a positive integer, default 1.
Other defaults remain caller-specific: selecting generation 2 is a scientific
change, not a promise of byte-identical output.

Generation 1 retains the historical command, including its version-specific
`--report_seq`/`--report-seq`. Generation 2 uses `--input`, `--vcf`, `--reference`,
`--threads`, `--minsupport`, `--output-rnames`, and `--allow-overwrite`.
Genome and BAM indexes are declared inputs. YAML files are explicit caller
dependencies so editing the caller choice invalidates the VCF even on
Snakemake 5.10. Use separate work directories for comparisons; when changing
only command-line configuration overrides on older Snakemake, force the caller
with `--forcerun call_sniffles_outsider`.

## VCF contract

The sequence and read-name parsers accept `--caller sniffles2`, or detect
`##source=Sniffles2` in automatic mode. VCF format version alone cannot identify
the caller: Sniffles 2 and SVIM both emit VCF 4.2. The legacy parsing branches
remain unchanged for historical outputs.

- Read INFO by keys: `SVTYPE`, `SVLEN`, `SUPPORT`, `END`, precision flags and
  `RNAMES`, regardless of ordering or a trailing semicolon.
- For sequence-resolved INS/DEL, remove the reference anchor using SVLEN and
  the opposing allele. Unanchored insertion ALT without a reference is also
  accepted when its length already equals SVLEN.
- VCF POS numerically equals the zero-based insertion breakpoint after the
  preceding base; it is passed unchanged to the current downstream contract.
- Keep original IDs within `sniffles.<TYPE>.<original_id>` so all existing
  evidence, frequency, TSD, and integration joins keep the same schema.
- Preserve N bases in generation 2 sequences. Apply the existing strict
  `length > -m` threshold and normalized INS/DEL exclusions.
- Skip symbolic alleles without sequence and report their number in the log;
  never manufacture DNA from `<INS>`/`<DEL>`. Reject malformed sequence lengths,
  duplicate IDs and multiallelic records. BND/INV/DUP do not feed TE extraction.
- Empty VCF/sequence/BLAST branches emit valid empty or header-only artifacts.

## References

- [Sniffles 2.8.1 release](https://github.com/fritzsedlazeck/Sniffles/releases/tag/v2.8.1)
- [Versioned VCF writer](https://github.com/fritzsedlazeck/Sniffles/blob/v2.8.1/src/sniffles/vcf.py)
- [Versioned CLI options](https://github.com/fritzsedlazeck/Sniffles/blob/v2.8.1/src/sniffles/config.py)
- [Python dependencies](https://github.com/fritzsedlazeck/Sniffles/blob/v2.8.1/setup.cfg)

## Validation

Both generations completed fresh 70-job runs with all scientific branches and
Quarto enabled, using the existing image and isolated temporary Sniffles 2
installation (Python 3.13 for the container execution; also smoke-tested with
Python 3.12 on the host). The image definition itself has not been rebuilt.
Sniffles 1 passed the historical exact-output regression (434 final calls).
Sniffles 2.8.1 produced 49 SVs, 5 retained consensus insertion sequences and
431 final TE calls. Consensus lengths, anchors, IDs, coordinates and support
were compared to the native VCF; the complete report regression also passed.
No Sniffles 2 candidate in this small fixture passes the TE classifier thresholds;
the 431 calls come from INSIDER and other OUTSIDER evidence. A separate synthetic
integration test verifies a normalized native Sniffles 2 ID and sequence through
genome reconstruction. These counts do not measure comparative biological accuracy.

An additional fresh 70-job run with MIN_SUPPORT=1000000 validated the zero-SV
branch through integration, Liftoff and Quarto. Unit tests cover both command
syntaxes, version mismatch, caller failure, anchored/unanchored alleles,
deletions, symbolic/empty inputs, malformed lengths, IDs and read-name parsing.
