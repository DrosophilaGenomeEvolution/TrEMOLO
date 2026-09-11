# Migration validation — 2026-09-11

The refactored `workflow/Snakefile` completed all 70 jobs from a fresh output
directory on the bundled fixture. INSIDER, OUTSIDER, TE_GENOME, OUTSIDER genome
integration, Liftoff and the Quarto report were enabled.

## Environment and execution

- Local image: `../TrEMOLO-test.simg`, Snakemake 5.10.0, Python 3.8.10.
- Liftoff 1.6.3: `/opt/conda/envs/liftoff_env/bin/liftoff`.
- Host Quarto 1.9.36 mounted read-only at `/opt/quarto`.
- Eight cores; clean execution completed in approximately 66 seconds.
- Configuration: `../migration_validation_20260911/config.yml`.
- Outputs: `../migration_validation_20260911/clean_output`.
- Log: `../migration_validation_20260911/clean_run.log`.

These paths identify local validation artifacts, not tracked repository files.
The first attempt exposed Liftoff missing from the image PATH. The fixture now
sets its executable explicitly, and the container definition exports its bin
directory. Validation restarted in a separate empty directory after this fix.
The image itself was not rebuilt; the successful run uses the existing image
and the read-only Quarto bind documented in the README.

## Results

- 101 workflow unit tests, 7 module tests and 5 locus-model tests passed.
- `check_variant_calling.py` passed against `../work_test`.
- The final 434-row `TE_INFOS.bed` is byte-identical to the historical oracle.
- All 25 canonical/observed OUTSIDER integrations passed validation; the
  observed reconstructed genome remains byte-identical to the oracle.
- Liftoff produced 10 accepted projections and 13 audited rejected mappings.
- `check_quarto_report.py` passed: 434 calls, 56 families, 7 proximity groups.
- A subsequent dry run without an explicit target reported `Nothing to be done`.

This validates compatibility on the bundled fixture, not additional datasets
or a freshly built container. SyRI and consensus backends remain unimplemented.
See `tests/regression/README.md` for clean-run instructions.
