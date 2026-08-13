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
present and every CIGAR insertion must have one matching TSD candidate. The
legacy multiprocessing code lost 20 such records in the current oracle.

## Current oracle status

Both pipelines completed in this reference run: 377 calls are labelled INSIDER
and 57 are labelled OUTSIDER. The generated report currently attempts an
external `curl`; its network failure did not prevent completion and is recorded
as a known reporting limitation.
