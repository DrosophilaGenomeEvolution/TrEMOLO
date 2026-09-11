# Historical `work_test` regression

This test treats `../work_test` from the `master` branch as the first complete
TrEMOLO regression oracle without
copying its 712 MB of generated files into the repository.

Run the semantic comparison from the repository root:

```bash
python3 tests/regression/check_work_test.py
```

To check another output directory:

```bash
python3 tests/regression/check_work_test.py /path/to/output
```

The default comparison checks stable scientific properties: total and
INSIDER/OUTSIDER call counts, chromosomes, TE families, event types, TSD count,
output presence, and line counts. Use `--strict` to additionally require an
exact `TE_INFOS.bed` SHA-256 checksum.

During workflow migration, compare the refactored scientific intermediates to
the same oracle with:

```bash
python3 tests/regression/check_variant_calling.py \
  ../work_test ../work_refactor_test
```

This comparison requires exact hashes for stable variant-calling and TE
detection outputs. For `INS_FOR_TSD.txt`, every legacy record must remain
present and every CIGAR insertion must have one matching raw flank candidate.
The legacy multiprocessing code lost 20 such records in the current oracle.
This file was not consumed by the historical `TSD_OUTSIDER` rule; its recovery
therefore does not change `TSD_TE.tsv` during compatibility migration.

The migrated OUTSIDER flank/TSD chain additionally requires byte-identical
`SIZE_SEQ.tsv`, `ALL_FK_REPORT_FT*.bed`, and `TSD_TE.tsv` files. Its explicit
`TSD_FLANK_ELIGIBILITY.tsv` report must contain the expected 18 accepted and
8 rejected candidates with explicit reasons.

The migrated INSIDER flank/TSD chain pins all four historical text artifacts:
377 merged candidates, 754 flank intervals, 377 formatted flank pairs, and 210
TSD calls. The final `TE_INFOS.bed` is also required byte for byte, including
its 434 records, 15-column schema, source order, TSD shifts, frequencies, and
legacy TrEMOLO identifiers.

The optional `INTEGRATE_TE_TO_GENOME` branch intentionally departs from the
invalid canonical legacy FASTA: it no longer inserts event IDs as DNA and now
reverse-complements negative-strand canonical TEs. Its observed `NEO_GENOME`
reconstruction remains byte-identical to `work_test`. Liftoff projections are
reported with per-event states in `LIFT_OFF_AUDIT.tsv`; discordant chromosome
pairs cannot enter the public reference BED. See decision 0015 for the exact
contract.

The Quarto report is checked as a lossless presentation of that final table:

```bash
python3 tests/regression/check_quarto_report.py ../work_refactor_test
```

This verifies every summary against `TE_INFOS.bed`, the embedded input checksum,
the interactive call and nearby-candidate views, and the absence of external
script or stylesheet dependencies. The HTML is intentionally not compared byte
for byte with the old R Markdown report because the interface has been replaced.

The optional module refactor is validated after generating both standalone
reports from the bundled population fixture and `work_refactor_test`:

```bash
python3 tests/regression/check_modules.py \
  /tmp/tremolo-population-refactor \
  /tmp/tremolo-structure-cli
```

The migrated OUTSIDER frequency chain is checked byte for byte from
`TE_SIZE.tsv` and its combined 69-candidate table through `COUNT_READS.txt`,
`FREQUENCY_TE_INS.tsv`, and `FREQUENCY_TE_INS_PRECISE.tsv`. Semantic checks
also pin the 34 legacy frequency calls (29 INS and 5 DEL), including the 25
merged calls and 9 compatibility-only extras. Known scientific limitations of
that compatibility result are documented in
`docs/decisions/0004-outsider-frequency-compatibility.md`.

## Current oracle status

Both pipelines completed in this reference run: 377 calls are labelled INSIDER
and 57 are labelled OUTSIDER. The historical R Markdown report attempts an
external `curl`; the refactored Quarto report does not perform network access.

## Clean end-to-end validation

Run from the parent directory of the checkout. Copy
`TrEMOLO/tests/workflow/refactor_config.yml` to a separate YAML file and change
`DATA.WORK_DIRECTORY` to a new, empty output path. Keep the bundled input paths.
The fixture enables INSIDER, OUTSIDER, resident annotation, genome integration,
Liftoff, and the Quarto report.

```bash
singularity exec --bind /opt/quarto:/opt/quarto:ro TrEMOLO-test.simg snakemake \
  --snakefile TrEMOLO/workflow/Snakefile \
  --configfile /path/to/validation.yml --cores 8 all
python3 TrEMOLO/tests/regression/check_variant_calling.py work_test /path/to/new-output
python3 TrEMOLO/tests/regression/check_quarto_report.py /path/to/new-output
```

The bind is required only for an older image without Quarto. Adjust
`TOOLS.QUARTO` and `TOOLS.LIFTOFF` to the paths inside your image. The local
Snakemake 5.10 image has Liftoff at
`/opt/conda/envs/liftoff_env/bin/liftoff`, outside its default PATH.
Do not reuse a completed output directory as evidence of a clean run.
