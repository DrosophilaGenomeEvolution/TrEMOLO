# Configuration reference

Use the minimal examples in [examples/config](../examples/config) to start a new analysis.
Paths in YAML are resolved from the execution working directory, not the YAML directory.
Use a new output directory when changing pipeline versions or scientific settings.
Use paths without spaces: some shell fragments still interpolate paths without quoting.

For a complete copy-and-edit example and launch command, start with
[Run with your own data](../README.md#run-with-your-own-data). The longer YAML
below illustrates optional settings; it is not the minimum needed to run.

| Mode | Enable under `CHOICE.PIPELINE` | Required `DATA` inputs |
| --- | --- | --- |
| OUTSIDER | `OUTSIDER_VARIANT: true` | `GENOME`, `TE_DB`, `SAMPLE` |
| INSIDER | `INSIDER_VARIANT: true` | `GENOME`, `TE_DB`, `REFERENCE` |
| Combined | Both variant branches | `GENOME`, `TE_DB`, `REFERENCE`, `SAMPLE` |
| Resident catalogue only | `TE_GENOME: true`, both variant branches false | `GENOME`, `TE_DB` |

Set `WORK_DIRECTORY` to a new output path and `REPORT: true` to generate HTML.
For OUTSIDER, `GENOME` is the read-mapping assembly. For INSIDER, it is the
sample assembly compared with `REFERENCE`. INSIDER empty-site frequency needs
OUTSIDER read mappings; assembly-only runs do not estimate it.


In the v3 workflow, `CHROM_KEEP` always searches the chromosome name
alone. For exact main dm6 chromosomes, use
`CHROM_KEEP: '^([23][RL]|[X4Y])$'`; descriptions, whitespace and other BED columns
are not searched. See [chromosome filtering](chromosome-filtering.md) for
parameter scopes, regex syntax and compatibility with older comma-separated lists.

The migrated workflow can also catalogue TE matches in reads with no genomic
alignment using `CHOICE.PIPELINE.UNMAPPED_TE: true` (disabled by default).
Results are written under `UNMAPPED_TE/`, with original TE labels and coordinates
within each read. See [configuration and column descriptions](unmapped-te.md).

Ambiguous variable TE calls retain their sequence/read provenance in
`TE_CALL_CANDIDATES.tsv` and the companion `TE_CALL_CANDIDATE_EVIDENCE.tsv`.
The REPORT distinguishes annotations on the same query from evidence on different
insertion sequences; see [column descriptions and limitations](te-call-candidates.md).

A reference-backed TE insertion GFA can be built as an explicit post-processing
step, with representatives ranked by BLAST bitscore and synthetic walks by level.
See [inputs, selection rules and output files](te-gfa.md).
Each graph build includes an HTML/JSON representation balance and per-call
outcomes. See also [real-run regression capture](real-run-regression.md)
and the [functional release candidate notes](release-notes-functional-candidate.md).

TrEMOLO uses Snakemake. Set parameters in YAML using the examples below or the minimal configurations. Optional settings omitted from a minimal example use workflow defaults. `INTERMEDIATE_FILE` is a legacy setting; the migrated workflow declares its intermediate outputs in the DAG.

```yaml
# all path can be relative or absolute depending of your tree.
#It is advised to only use absolute path if you are not familiar with computer science or the importance of folder trees structure.
DATA:
    GENOME:          "/path/to/genome_file.fasta"      #genome (fasta file) [required]
    TE_DB:           "/path/to/database_TE.fasta"      #Database of TE (a fasta file) [required]
    REFERENCE:       "/path/to/reference_file.fasta"   #reference genome (fasta file) only if INSIDER_VARIANT = True [optional]
    SAMPLE:          "/path/to/reads_file.fastq"       #long reads (a fastq[.gz] file) only if OUTSIDER_VARIANT = True [optional]
    # REFERENCE and SAMPLE are required only by their enabled branches.
    WORK_DIRECTORY:  "/path/to/directory"         #name of output directory [optional, will be created as 'TrEMOLO_OUTPUT']

# For a resident-only catalogue, neither REFERENCE nor SAMPLE is required.

CHOICE:
    PIPELINE:
        OUTSIDER_VARIANT: True  # outsiders, TE not in the assembly - population variation
        INSIDER_VARIANT: True   # insiders, TE in the assembly
        TE_GENOME: False        # resident TE annotation across complete assemblies
        REPORT: True            # for getting a report.html file with graphics
    OUTSIDER_VARIANT:
        CALL_SV: "sniffles1"    # sniffles1, sniffles2; sniffles is the legacy alias
        INTEGRATE_TE_TO_GENOME: True # (True, False) Re-build the assembly with the OUTSIDER integrated in
        INTEGRATION_ALLOW_PARTIAL: True # warn and skip rejected INS; False enables strict failure
        CLIPPED_READS: False # (True, False) Processing of clipped reads (SOFT, HARD)


PARAMS:
    THREADS: 8 #number of threads for some task
    TE_GENOME:
        TARGETS: ["GENOME"] # optionally add REFERENCE
        CHROM_KEEP: "."
        MIN_PIDENT: 65
        MIN_ALIGNED_BP: 80
        MAX_EVALUE: 1e-10
        FULL_LENGTH_COVERAGE: 80 # classification only; partial copies are retained
        PARTIAL_COVERAGE: 20
    OUTSIDER_VARIANT:
        MINIMAP2:
            PRESET_OPTION: 'map-ont' # minimap2 option is map-ont by default (map-pb, map-ont)
            OPTION: '' # more option of minimap2 can be specified here
        SAMTOOLS_VIEW:
            PRESET_OPTION: ''
        SAMTOOLS_SORT:
            PRESET_OPTION: ''
        SAMTOOLS_CALLMD:
            PRESET_OPTION: ''
        TSD:
            SIZE_FLANK: 15  # flanking sequence size for calculation of TSD; put value > 4
        LIFT_OFF:
            FLANK_SIZE: 100000 # flank projected on each side of an integrated TE
            MAX_GAP: 20000     # maximum accepted interval/overlap on the reference
        TE_DETECTION:
            CHROM_KEEP: "." # regular expresion for chromosome filtering; for instance for Drosophila  "2L,2R,3[RL],X" ; Put "." to keep all chromosome
            GET_SEQ_REPORT_OPTION: "-m 30" #sequence recovery file in the vcf
            TIME_LIMIT: 0 # time limit for processing SV detection (rule TrEMOLO_SV_TE); put value >= 0 (hours); 0 means no time limit
        PARS_BLN_OPTION: "--min-size-percent 80 --min-pident 80 -k 'INS|DEL'" # option for TrEMOLO/lib/python/parse_blast_main.py - don't put -c option
    INSIDER_VARIANT:
        PARS_BLN_OPTION: "--min-size-percent 80 --min-pident 80" # parameters for validation of insiders
        MINIMAP2:
            PRESET_OPTION: 'asm5' # minimap2 preset option is asm5 by default (asm5, asm10, asm20 etc)
            OPTION: '--cs'


```

When `INTEGRATE_TE_TO_GENOME` is enabled, TrEMOLO writes both the canonical
reconstruction (`OUTSIDER/TE_TOWARD_GENOME/PSEUDO_GENOME_TE_DB_ID.fasta`) and
the reconstruction based on observed SV sequences (`NEO_GENOME.fasta`). The
corresponding shifted BED files and `INTEGRATION_TE.tsv` make every inserted,
rejected or excluded event explicit. Only INS events are added; DEL events are
excluded because their sequence is already present in the input assembly.
Observed sequences exclude the read flanks used for BLAST/TSD. Alternatives
at the same site are concatenated intentionally for this population
representation. An unresolved INS is skipped with a warning by default; valid
insertions are reconstructed and the integration step succeeds. Set
`CHOICE.OUTSIDER_VARIANT.INTEGRATION_ALLOW_PARTIAL: false` to fail instead of
accepting a partial reconstruction. The standalone integration script uses the
same default; `--strict` enables failure and `--allow-partial` remains supported.
The report shows eligible, integrated,
rejected and excluded counts. If INSIDER is also enabled, the two flanks of each
integrated insertion are projected to the reference with Liftoff; concordant
calls are written to `POS_TE_OUTSIDER_ON_REF.bed`, while uncertain mappings and
their reasons remain available in `BAD_POS_TE_LIFT.bed` and
`OUTSIDER/INSIDER_VR/LIFT_OFF_AUDIT.tsv`. Liftoff flanks use one-based inclusive
GFF coordinates; output BED coordinates are zero-based and half-open, with
orientation-aware validation of the flank pairs.
The Liftoff audit contains one row per integrated INS, including
`rejected / no_mapped_flanks` when neither flank is found. These events remain
in the reconstructed genomes and have no public projected reference position.

The v3 workflow supports two Sniffles generations. Select one in YAML:

```yaml
CHOICE:
    OUTSIDER_VARIANT:
        CALL_SV: sniffles2   # sniffles1 (default), sniffles2; sniffles aliases sniffles1
PARAMS:
    OUTSIDER_VARIANT:
        SNIFFLES:
            MIN_SUPPORT: 1
TOOLS:
    SNIFFLES1: sniffles     # optional executable override
    SNIFFLES2: sniffles2    # optional executable override
```

Merge these settings into the existing sections; do not duplicate YAML keys.
`TOOLS.SNIFFLES1` and `TOOLS.SNIFFLES2` may be omitted: their defaults are
`sniffles` and `sniffles2`, respectively. Only the selected caller needs to be
available on PATH (or at its configured path) when OUTSIDER variant calling runs.
`REPORT.AUTHOR` is optional as well: omit it, leave it blank, or set `""` to
generate a report without an author.
The container definition preserves Sniffles 1.0.12b and adds Sniffles 2.8.1 in
an isolated Python 3.12 environment. Rebuild the image to obtain `sniffles2`.
Use `workflow/Snakefile` for this choice; the 2.x entry point has been removed.
`OUTSIDER/VARIANT_CALLING/sniffles-run.json` records the actual version and command.
Sniffles 2 uses `--input`, `--reference`, `--minsupport`, and `--output-rnames`;
its REF/ALT sequences and supporting reads are normalized before TE analysis.
Symbolic alleles without sequence are counted in the extraction log and skipped.
The algorithms and other defaults differ, so equal support thresholds do not
imply identical calls. Use separate work directories when comparing callers.
See [the Sniffles backend contract](decisions/0016-sniffles-backends.md).

The main parameters are:

- `GENOME` : Assembly of the sample of interest (or mix of samples), fasta file.
- `TE_DB`  : A **Multifasta** file containing the canonical sequence of transposable elements. You can add also copy sequences but results will be more complex to interpretate.
- `REFERENCE` : Fasta file containing the reference genome of the species of interest.
- `WORK_DIRECTORY` : Directory that will contain the output files. If the directory does not exist it will be created;  default value is **TrEMOLO_OUTPUT**.
- `SAMPLE` : File containing the reads used for the sample assembly.


Use [insider.yaml](../examples/config/insider.yaml) for assembly-only INSIDER,
[outsider.yaml](../examples/config/outsider.yaml) for OUTSIDER, or
[combined.yaml](../examples/config/combined.yaml) for both. These minimal v3
examples are the recommended starting points; the root `config*.yaml` files
retain historical settings.
To analyse **INSIDER**, only the `REFERENCE` , the `GENOME`, the `TE_DB` and the `WORK_DIRECTORY` are required.
To analyse **OUTSIDER**, only the `SAMPLE` , the `GENOME`, the `TE_DB` and the `WORK_DIRECTORY` are required.


