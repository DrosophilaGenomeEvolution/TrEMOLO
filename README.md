# TrEMOLO

<img src="images/TrEMOLO9.png" alt="TrEMOLO logo">

**Transposable Elements MOvement detection using LOng reads**

TrEMOLO detects TE-associated variation using long reads and genome assemblies.
INSIDER compares assemblies; OUTSIDER maps reads to an assembly to identify
variation within the sequenced sample or population.

**Current development version: `3.0.0-beta.1`.** This is a proposed beta release;
no release tag or new downloadable image is published by this documentation.
The supported entry point is `workflow/Snakefile`. Read the
[beta release notes](docs/release-notes-functional-candidate.md) before migrating
from 2.x: corrected evidence extraction and integration intentionally change
some scientific results. Use a fresh output directory.

## Quick start

Use an image built from the current [Singularity recipe](Singularity).
The local development image is named `TrEMOLO-3.simg`; that filename
alone does not establish its version. The old 2.5.4b download is a legacy image
and is not the recommended environment for this beta.

The commands below run **from the directory containing the TrEMOLO checkout
and the image**. They use the bundled small fixture. Start with a fresh
`work_refactor_test` directory; do not remove an existing analysis to rerun it.

```bash
git clone https://github.com/DrosophilaGenomeEvolution/TrEMOLO.git
# If the checkout already exists, use it instead of cloning again.
# Build only if you do not already have a current image:
sudo singularity build TrEMOLO-3.simg TrEMOLO/Singularity

singularity exec TrEMOLO-3.simg python3 TrEMOLO/lib/python/tremolo_version.py
singularity exec TrEMOLO-3.simg snakemake \
  --snakefile TrEMOLO/workflow/Snakefile \
  --configfile TrEMOLO/tests/workflow/refactor_config.yml \
  --cores 8 all
```

Open `work_refactor_test/REPORT/report.html` in a browser. It needs no server or
network connection. Apptainer can replace `singularity` where installed.
For data outside the container's default mounts, bind the relevant directory
and use paths visible inside the container.

For your own data, copy an [example configuration](examples/config), replace
input/output paths and pass it with `--configfile`. YAML paths resolve from the
working directory of the command. Prefer absolute paths for your own analyses.
`--cores` limits the total scheduled CPU allocation; `PARAMS.THREADS` requests
threads for individual jobs.

## Choose an analysis

| Goal | Required scientific inputs | Mode or command | Main result |
| --- | --- | --- | --- |
| Compare assemblies | Assembly, reference, TE database | `INSIDER_VARIANT` | Variant-associated TE calls |
| Detect within-sample variation | Reads, mapping assembly, TE database | `OUTSIDER_VARIANT` | Variant-associated calls and read frequencies |
| Catalogue resident TEs | Assembly, TE database | `TE_GENOME` | Copy/fragment tables, BED, GFF3 |
| Examine unmapped reads | Reads, assembly, TE database | `UNMAPPED_TE` | Read-coordinate TE evidence |
| Compare timepoints/replicates | Call tables on the same reference, reference FASTA | Module 1 | Frequency trajectories |
| Inspect insertion structures | OUTSIDER direct-channel BLAST and size tables | Module 2 | Insertion evidence explorer |
| Export an insertion graph | Calls, raw insertion evidence, reference, original TE database | GFA exporter | GFA, representation balance and audit tables |

Optional features have independent settings. See the
[configuration reference](docs/configuration.md), [workflow targets](workflow/README.md)
and [documentation index](docs/README.md).

## Version and reproducibility

`VERSION` is the single source of the application version. Workflow startup,
module startup and the GFA command display the version, full Git commit and a
modified-checkout marker when applicable. Git is queried in the source checkout,
not the container's bundled checkout. An archive without Git reports the commit
as unavailable; it never invents a commit identifier.

```bash
python3 TrEMOLO/lib/python/tremolo_version.py
python3 TrEMOLO/lib/python/tremolo_version.py --json
TrEMOLO/module --version
python3 TrEMOLO/lib/python/locus/build_te_gfa.py --version
```

An executing Snakemake workflow writes `WORK_DIRECTORY/log/run-version.json`.
A dry run displays the banner without creating this file. The file describes the
latest execution, including resumed jobs; it does not establish the origin of
all pre-existing outputs. Version identity does not replace input checksums,
configuration, tool versions or the container image checksum. See
[recording a reproducible run](docs/real-run-regression.md).

## Installation and requirements

The current image includes the scientific stack, Sniffles 1 and 2, Liftoff and
Quarto. Check the actual image before an analysis:

```bash
singularity exec TrEMOLO-3.simg snakemake --version
singularity exec TrEMOLO-3.simg sniffles --help
singularity exec TrEMOLO-3.simg sniffles2 --version
singularity exec TrEMOLO-3.simg /opt/quarto/bin/quarto --version
singularity exec TrEMOLO-3.simg /opt/conda/envs/liftoff_env/bin/liftoff --version
```

The migrated workflow and modules use Quarto and browser JavaScript for reports.
R/Bookdown requirements belong to the historical `run.snk` report. A persistent
Node server is not required. For native installations, the recipe records the
scientific dependencies; the tested legacy workflow environment uses Snakemake
5.10 and Python 3.8. Do not assume compatibility with every newer Snakemake release.

Editing the recipe does not update an existing image. With an older image lacking
Quarto, see the [workflow environment notes](workflow/README.md); building the
current image is the recommended route for a new installation.

## Results and interpretation

Start with `TE_INFOS.bed`, `REPORT/report.html` and the original-name exports in
`PUBLIC/`. The [output guide](docs/outputs.md) describes their schemas and optional
branch outputs. The [interpretation guide](docs/interpreting-results.md) explains
missing calls, frequency, alternative families and coordinate systems.

Integration can skip unresolved insertions with warnings by default; inspect
`INTEGRATION_TE.tsv` and the report. Set `INTEGRATION_ALLOW_PARTIAL: false` for
strict failure. The reconstructed population genomes concatenate same-site
alternatives and should not be interpreted as observed haplotypes.

## Modules

Run these commands from the parent of the checkout, using paths visible inside
the container. Module 1 expects a TSV manifest; its input BED files must share
the same reference coordinates.

```bash
singularity exec TrEMOLO-3.simg TrEMOLO/module build 1 \
  -i samples.tsv -g reference.fasta -o POPULATION_REPORT

singularity exec TrEMOLO-3.simg TrEMOLO/module build 2 \
  work_refactor_test STRUCTURE_REPORT
```

- [Module 1: population frequency trajectories](modules/1-FREQUENCY-MULTI-GENERATIONS/README.md).
  Its HTML is self-contained. Missing calls remain missing.
- [Module 2: insertion structure explorer](modules/2-MODULE_TE_BLAST/README.md).
  Share the whole output directory: `report-assets/` must remain beside the HTML.
  Structure classifications are provisional hypotheses.

![Module 2 insertion evidence explorer](modules/2-MODULE_TE_BLAST/img/insertion-structure-real-example.png)

## Build a GFA

This is an explicit post-processing command, independent of the `all` target:

```bash
singularity exec TrEMOLO-3.simg python3 TrEMOLO/lib/python/locus/build_te_gfa.py \
  --reference work_refactor_test/INPUT/genome.fasta \
  --te-database TrEMOLO/test/canonical_TE.fa \
  --reference-name fixture \
  --run fixture=work_refactor_test \
  --output work_refactor_test/GFA
```

Open `GFA/summary.html` to inspect represented and excluded calls. Supported
OUTSIDER insertions become alternatives on synthetic, unphased walks. Defaults
are 80% identity, 80% TE coverage and 100 aligned bp; they are independent of
pipeline thresholds. Only chromosomes present in the input BED are exported by
default. See [GFA selection, limitations and output files](docs/te-gfa.md).

## Migration and validation

The historical `run.snk` remains for compatibility; it is not the supported beta
entry point. Shared parsing scripts have changed, so running it from this
checkout does not reproduce a frozen 2.x analysis.

See [regression procedures](tests/regression/README.md),
[historical migration validation](docs/migration-validation.md) and
[beta release validation status](docs/release-notes-functional-candidate.md).
The historical 434-call oracle predates scientific corrections; counts from a
current run can differ intentionally. SyRI/consensus backends and full Savanache
application compatibility are not included in the beta's validated scope.

## Other utilities

[Assembly correction utilities](docs/assembly-correction.md) reorder and orient
assembly blocks with audited base retention. Historical strategy diagrams are
available under [images/HOW-WORKING-TrEMOLO](images/HOW-WORKING-TrEMOLO).

## Citation and licence

Documentation is licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
See [CeCILL-C](Licence_CeCILL-C_V1-en.txt) and [GPLv3](LICENSE) for software licences.

Please cite Mohamed et al., **TrEMOLO: accurate transposable element allele
frequency estimation using long-read sequencing data combining assembly and
mapping-based approaches**, Genome Biology 24, 63 (2023),
[doi:10.1186/s13059-023-02911-2](https://doi.org/10.1186/s13059-023-02911-2).

See also Mohamed et al., **A Transposon Story: From TE Content to TE Dynamic
Invasion of Drosophila Genomes Using the Single-Molecule Sequencing Technology
from Oxford Nanopore**, Cells 9, 1776 (2020),
[doi:10.3390/cells9081776](https://doi.org/10.3390/cells9081776).
[Paper datasets](https://dataverse.ird.fr/dataverse/tremolo_data).
