# Freezing a real analysis for regression

`tests/regression/real_run.py` records observed outputs and compares subsequent
runs to them. It does not execute the pipeline, regenerate expected results or
silently approve scientific changes. Keep baselines outside the analysis output
directory so a rerun cannot replace them. Existing baselines cannot be overwritten
by the `snapshot` command; use a new filename after reviewing an intended change.

## Record the execution environment

Run this inside the same container/job environment used for the analysis:

```bash
python3 TrEMOLO/tests/regression/real_run.py versions \
  --tool minimap2=minimap2 --tool samtools=samtools \
  --tool sniffles=sniffles --tool blastn=blastn \
  --tool snakemake=snakemake --tool quarto=quarto \
  --output versions.json
```

The recorder stores Python's version, the capture checkout's Git commit and dirty
state when Git is available, executable commands, version output and return codes.
Use a clean checkout for a release reference; a commit alone does not preserve
uncommitted changes. It returns failure if a tool cannot be queried. Supply paths for
executables not on PATH, and select the tools actually used by the run. A failed
version query cannot be accepted as a complete version record. Versions queried
on another machine after the run do not establish the original environment.

## Freeze a complete run record

```bash
python3 TrEMOLO/tests/regression/real_run.py snapshot \
  --run /path/analysis \
  --gfa main=/path/analysis/GFA \
  --config /path/exact-config.yaml --versions versions.json \
  --input reads=/path/reads.fastq.gz \
  --input genome=/path/genome.fasta \
  --input te_database=/path/TE_database.fasta \
  --command 'singularity exec TrEMOLO.simg snakemake all --snakefile TrEMOLO/workflow/Snakefile --configfile exact-config.yaml --cores 8' \
  --output /path/baselines/real-run-v1.json
```

Include all scientific inputs actually used, adding an assembly/reference or the
container image with additional `--input NAME=FILE` arguments as needed. Large
inputs are hashed without copying them. The exact YAML text and version record
are embedded in the baseline. Record the actual pipeline command and separately
retain the build-GFA command when it is run as a post-processing step; graph
selection/threshold settings are also copied from its manifest.

The baseline records available TE_INFOS/candidate/evidence tables, the GFA and its
sequence/locus/exclusion/walk tables. GFA file checksums cover node sequences,
orientations, paths and links. Table checksums ignore row order but require exact
field values, keeping duplicate rows significant. Consequently a changed sequence,
position, family, orientation or support count is detected even if aggregate
totals remain the same. Numeric formatting changes also require review.

`reproducibility.complete` means that the run, configuration, version record,
explicit input inventory and command were all supplied. It is an inventory
check, not independent proof that the supplied metadata describe the original
execution or that every scientific input/tool has been listed.

## Compare a rerun

```bash
python3 TrEMOLO/tests/regression/real_run.py check \
  --baseline /path/baselines/real-run-v1.json \
  --run /path/new-analysis --gfa main=/path/new-analysis/GFA
```

Use the same named graph/input roles; their paths may change. Supply `--input`,
`--config` and `--versions` to check files at new locations. Otherwise the original
locations in the baseline are used. Configuration contents are compared exactly.
Use `--strict` to require byte-identical tables as well. The checker lists changed
fields/files and returns a nonzero exit status when results differ or a required
file cannot be read. The recorded pipeline Git revision/dirty state is descriptive
and is not compared: a regression run intentionally tests a new code revision.
It does not assume that a difference is necessarily a bug;
review the scientific reason before creating a replacement baseline.

## Available dm6/G73 output snapshot

`tests/regression/baselines/dm6-g73-available-gfa.json` freezes the locally
available ZAM, gtwin and ZAM+gtwin GFA outputs under `9-WORK`. It deliberately
records `available_outputs_only`: the complete source run, exact YAML, original
tool versions, raw inputs and execution command are not available locally.
Input hashes already present in the graph manifests are preserved as recorded
provenance, without claiming that the inaccessible `/scratch/...` inputs were
verified here. This snapshot can detect output changes but cannot reproduce the
whole biological analysis on its own.

These frozen graphs predate the BED-only chromosome default. A comparable rebuild
must use `--all-reference-chromosomes` with their original selections/thresholds.
A rebuild with the new default can intentionally differ; review it and create a
new baseline rather than replacing the historical output snapshot silently.

```bash
python3 tests/regression/real_run.py check \
  --baseline tests/regression/baselines/dm6-g73-available-gfa.json
```
