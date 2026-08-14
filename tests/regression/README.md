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

The migrated OUTSIDER frequency chain is checked byte for byte from
`TE_SIZE.tsv` and its combined 69-candidate table through `COUNT_READS.txt`,
`FREQUENCY_TE_INS.tsv`, and `FREQUENCY_TE_INS_PRECISE.tsv`. Semantic checks
also pin the 34 legacy frequency calls (29 INS and 5 DEL), including the 25
merged calls and 9 compatibility-only extras. Known scientific limitations of
that compatibility result are documented in
`docs/decisions/0004-outsider-frequency-compatibility.md`.

## Current oracle status

Both pipelines completed in this reference run: 377 calls are labelled INSIDER
and 57 are labelled OUTSIDER. The generated report currently attempts an
external `curl`; its network failure did not prevent completion and is recorded
as a known reporting limitation.
