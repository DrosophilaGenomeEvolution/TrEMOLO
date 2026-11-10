# TrEMOLO v3 compatibility audit — 2026-11-10

## Conclusion and scope

Source inspection covers the executable branches selected by the 2.x launcher,
comparison with the explicit v3 DAG, existing compatibility decisions and local
tests. The main scientific branches have replacements. This is **not** a claim
of complete option parity or a fresh end-to-end scientific validation.

The v3 entry point is `workflow/Snakefile`. The generated-workflow launcher
`run.snk`, monolithic `rules.snk`, generation template `template.snk` and their
private shell environment `env.sh` are removed. Their previous implementations
remain in Git history. Shared scientific helpers, historical regression oracles,
input fixtures and existing analysis directories are retained.

## Branch coverage

| 2.x behavior / rules | v3 implementation | Compatibility |
| --- | --- | --- |
| Input preparation and `RENAME_PSEUDO_TE` | `common.smk`, `te_database.smk`, `te_names.smk` | Immutable normalized inputs, manifest/checksums, original names in public exports |
| `mapping`, `samtools`, `sniffles` | `outsider.smk` | Sniffles 1 retained; Sniffles 2 added; caller-specific results |
| `TrEMOLO_SV_TE`, `DETECTION_TE`, `SOFT_TE`, `HARD_TE`, `MERGE_TE` | `outsider_te.smk` | Four evidence sources; corrected cluster/candidate retention intentionally changes calls |
| `SV_INSIDER`, `TE_INSIDER` | `insider.smk` | Assemblytics backend; insertion/deletion and empty-result handling |
| `FREQUENCEv2`, `FREQ_INSIDERv2` | `outsider_frequency.smk`, `insider_frequency.smk` | Compatibility algorithms retained; documented headers and deterministic scheduling |
| `GET_SEQ_TE`, `TSD_OUTSIDER`, `TSD_INSIDER` | `outsider_tsd.smk`, `insider_tsd.smk` | Flanks and TSDs covered; selection fixes intentionally change historical OUTSIDER rows |
| `TE_ALL_IN_ASSEMBLY` | `te_genome.smk` | Replaced by resident-copy annotation; not byte-equivalent; separate from variant calls |
| `TE_TOWARD_GENOME`, `LIFT_OFF` | `outsider_integration.smk` | Insertion-only integration, oriented sequence and audited projection corrections |
| Final table and `REPORT` | `te_infos.smk`, `report.smk` | 15-column `TE_INFOS.bed` retained; offline Quarto replaces R/Bookdown |
| Optional modules | `module` and module builders | Offline reports replace previous interfaces/server |

INSIDER frequency estimation still requires the OUTSIDER mapping evidence.
The resident catalogue is not a replacement for variant-associated calls.
Historical byte comparisons remain useful only for stages whose scientific
semantics have not changed; see decisions 0004–0010 and 0017–0021.

## Explicitly retired or unavailable behavior

- `CALL_SV: svim` and `CALL_SV: no_sniffles` existed in the generated 2.x
  workflow. v3 rejects them; supported choices are `sniffles1`, `sniffles2` and
  the `sniffles` alias. These modes have **not** been migrated. Analyses requiring
  them must use a preserved previous checkout/environment.
- `FIND_SV_ON_REF` and `FIND_TE_ON_REF` were commented out in the 2.x launcher's
  generated instruction chain. Their removal does not remove an enabled default
  branch. `TE_INSIDERv2` was likewise not selected by that launcher.
- Dynamic Snakefile generation, progress-spinner processes, R/Bookdown reporting
  and its runtime network requests are retired implementation interfaces.
- SyRI and consensus remain planned backends, not migrated functionality.
- Full Savanache application/index/SWAP compatibility remains outside validated scope.

## Fresh validation in the scientific image — 2026-11-10

The agent successfully executed Singularity outside its restricted sandbox.
Earlier agent-session launch failures were sandbox restrictions, not evidence
of a broken host installation. The following checks supersede the earlier
missing-host-dependency and blocked-container records.

Image: `TrEMOLO-3.simg`, SHA-256
`30d76b50770ddd407d6f5bba74986d60c3bb4e17d0901d845e54c0964fce684a`.
Source: `39f91cac490edf07251a18727f6e10ee4709b832` **with uncommitted changes**.
This is validation of the working checkout, not a clean release commit.
The [machine-readable record](validation/2026-11-10-migration.json) records
runtime source fingerprints, exact configurations, tool versions, input hashes
and final table checksums. Local full outputs and logs are under
`/tmp/tremolo-audit-20261110/`; preserve them separately if needed.

- Full container suite: **245 tests passed, no skips**.
- Four standalone JavaScript suites pass: exploration, frequency/depth,
  resident/genome UX and pagination.
- Three fresh `all` executions completed, **74/74 steps each**. They enable
  INSIDER, OUTSIDER, resident annotation, unmapped-read analysis, integration,
  Liftoff and Quarto. The fixture config is `tests/workflow/v3_config.yml`,
  with distinct new output paths and the caller/support changes below.

| Fresh run | Sniffles SVs | Final TE calls | Families | Proximity groups |
| --- | --- | --- | --- | --- |
| Sniffles 1, 8 cores | 286 | 529 | 58 | 15 |
| Sniffles 2, 4 cores | 49 | 526 | 58 | 13 |
| Sniffles 2, `MIN_SUPPORT: 1000000`, 4 cores | 0 | 526 | 58 | 13 |

A Sniffles no-call result does not disable the independent CIGAR/clipped evidence
channels; the final table can therefore remain populated.

- Alignment provenance passes: 894 INS, 432 SOFT and 867 HARD records, with all
  INS/HARD clusters retained (Sniffles 1 run).
- Integration checks pass for all three runs: 23 integrated INS / 2 excluded
  events for Sniffles 1; 19 integrated INS / 0 excluded for both Sniffles 2 runs.
  Exact CIGAR sequences, flank removal, reconstructed genomes and Liftoff audit
  coordinates are checked.
- Sniffles 2 VCF normalization passes for both runs: 5 retained consensus
  sequences in the ordinary run; 0 in the no-call run. The regression checker
  was corrected to expect native `Sniffles2.INS.*` IDs, as already required by
  the adapter's unit tests; sequence/coordinate/support comparisons remain intact.
- Semantic Quarto regression and actual JavaScript initialization/filter checks
  pass on all three freshly generated reports. These are automated checks,
  not a manual visual review in a browser.
- All three post-run dry runs report **Nothing to be done**.
- The README's complete minimal YAML, with fixture paths substituted, passes
  a Snakemake dry run. Local documentation links resolve; `git diff --check` passes.
- Recipe help/labels no longer advertise v2.4 or an uninstalled Snakemake 7.32.4.

## Version and release gate

`3.0.0-beta.1` is appropriate: workflow/module/report interfaces and scientific
outputs change since 2.x. A major version communicates those changes; a beta
suffix keeps the remaining validation explicit. This audit does not publish a
release or create a tag. The bundled end-to-end, normalization, integration, report and idle-DAG gates now
pass in the local scientific image. Before publishing the beta, commit the reviewed
checkout and identify the exact distributed image with its checksum. A fresh
image build from the recipe was not performed in this audit: some dependencies
remain unpinned, so rebuilding is not guaranteed to reproduce this tested image.
Before a stable `3.0.0`, preserve and validate representative real-data execution
metadata as described in [release notes](release-notes-functional-candidate.md).
Full Savanache compatibility and SyRI/consensus remain outside this validation.
Paths containing whitespace are not reliably supported: some existing shell
fragments interpolate paths without quoting (including optional input arguments
in `common.smk`). The getting-started documentation now specifies paths without
spaces. Fixing this throughout all scientific rules remains separate work.
