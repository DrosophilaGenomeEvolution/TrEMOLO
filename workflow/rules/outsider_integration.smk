"""Optional OUTSIDER genome integration and reference projection rules."""

OUTSIDER_INTEGRATION_DIR = f"{WORKDIR}/OUTSIDER/TE_TOWARD_GENOME"
OUTSIDER_CANONICAL_GENOME = (
    f"{OUTSIDER_INTEGRATION_DIR}/PSEUDO_GENOME_TE_DB_ID.fasta"
)
OUTSIDER_OBSERVED_GENOME = f"{OUTSIDER_INTEGRATION_DIR}/NEO_GENOME.fasta"
OUTSIDER_CANONICAL_BED = f"{OUTSIDER_INTEGRATION_DIR}/TRUE_POSITION_TE_PSEUDO.bed"
OUTSIDER_OBSERVED_BED = f"{OUTSIDER_INTEGRATION_DIR}/TRUE_POSITION_TE_NEO.bed"
OUTSIDER_INTEGRATION_AUDIT = f"{OUTSIDER_INTEGRATION_DIR}/INTEGRATION_TE.tsv"
OUTSIDER_CANONICAL_PUBLIC_BED = (
    f"{WORKDIR}/POSITION_TE_OUTSIDER_IN_PSEUDO_GENOME.bed"
)
OUTSIDER_OBSERVED_PUBLIC_BED = f"{WORKDIR}/POSITION_TE_OUTSIDER_IN_NEO_GENOME.bed"
OUTSIDER_CANONICAL_INDEX = OUTSIDER_CANONICAL_GENOME + ".fai"

OUTSIDER_LIFT_DIR = f"{WORKDIR}/OUTSIDER/INSIDER_VR"
OUTSIDER_LIFT_INPUT_GFF = f"{OUTSIDER_LIFT_DIR}/INOUTSIDER.gff"
OUTSIDER_LIFT_FEATURES = f"{OUTSIDER_LIFT_DIR}/feature.txt"
OUTSIDER_LIFT_GFF = f"{OUTSIDER_LIFT_DIR}/output_INOUT.gff3"
OUTSIDER_LIFT_UNMAPPED = f"{OUTSIDER_LIFT_DIR}/unmapped_features.txt"
OUTSIDER_LIFT_INTERMEDIATES = f"{OUTSIDER_LIFT_DIR}/intermediate_files"
OUTSIDER_LIFT_IDS = f"{OUTSIDER_LIFT_DIR}/ID_FK_TE.txt"
OUTSIDER_LIFT_GOOD = f"{OUTSIDER_LIFT_DIR}/POS_TE_LIFT.bed"
OUTSIDER_LIFT_BAD = f"{OUTSIDER_LIFT_DIR}/BAD_POS_TE_LIFT.bed"
OUTSIDER_LIFT_AUDIT = f"{OUTSIDER_LIFT_DIR}/LIFT_OFF_AUDIT.tsv"
OUTSIDER_ON_REFERENCE = f"{WORKDIR}/POS_TE_OUTSIDER_ON_REF.bed"
POSITION_TE_ON_REFERENCE = f"{WORKDIR}/POSITION_TE_ON_REF.bed"

OUTSIDER_INTEGRATION_FINAL_OUTPUTS = [
    OUTSIDER_CANONICAL_GENOME,
    OUTSIDER_OBSERVED_GENOME,
    OUTSIDER_CANONICAL_BED,
    OUTSIDER_OBSERVED_BED,
    OUTSIDER_CANONICAL_PUBLIC_BED,
    OUTSIDER_OBSERVED_PUBLIC_BED,
    OUTSIDER_INTEGRATION_AUDIT,
]
OUTSIDER_LIFT_FINAL_OUTPUTS = [
    OUTSIDER_ON_REFERENCE,
    POSITION_TE_ON_REFERENCE,
    OUTSIDER_LIFT_BAD,
    OUTSIDER_LIFT_AUDIT,
]

LIFTOFF_PARAMS = OUTSIDER_PARAMS.get("LIFT_OFF", {})
OUTSIDER_LIFT_FLANK = LIFTOFF_PARAMS.get("FLANK_SIZE", 100000)
OUTSIDER_LIFT_MAX_GAP = LIFTOFF_PARAMS.get("MAX_GAP", 20000)
if type(OUTSIDER_LIFT_FLANK) is not int or OUTSIDER_LIFT_FLANK <= 0:
    raise ValueError("PARAMS.OUTSIDER_VARIANT.LIFT_OFF.FLANK_SIZE must be > 0")
if type(OUTSIDER_LIFT_MAX_GAP) is not int or OUTSIDER_LIFT_MAX_GAP < 0:
    raise ValueError("PARAMS.OUTSIDER_VARIANT.LIFT_OFF.MAX_GAP must be >= 0")
LIFTOFF_EXECUTABLE = config.get("TOOLS", {}).get("LIFTOFF", "liftoff")


rule integrate_outsider_te_into_genome:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/integrate_outsider_te.py"),
        genome=PREPARED_GENOME,
        te_database=PREPARED_TE_DATABASE,
        merged=OUTSIDER_MERGED_BED,
        sniffles_calls=OUTSIDER_SNIFFLES_CSV,
        direct_calls=OUTSIDER_INS_CSV,
        sniffles_fasta=OUTSIDER_SNIFFLES_FASTA,
        direct_fasta=OUTSIDER_INS_FASTA,
        tsd=OUTSIDER_TSD,
    output:
        canonical_genome=OUTSIDER_CANONICAL_GENOME,
        observed_genome=OUTSIDER_OBSERVED_GENOME,
        canonical_bed=OUTSIDER_CANONICAL_BED,
        observed_bed=OUTSIDER_OBSERVED_BED,
        canonical_public=OUTSIDER_CANONICAL_PUBLIC_BED,
        observed_public=OUTSIDER_OBSERVED_PUBLIC_BED,
        audit=OUTSIDER_INTEGRATION_AUDIT,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_integrate_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_integrate_te.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_INTEGRATION_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} \
            --genome {input.genome:q} \
            --te-database {input.te_database:q} \
            --merged-bed {input.merged:q} \
            --sniffles-calls {input.sniffles_calls:q} \
            --direct-calls {input.direct_calls:q} \
            --sniffles-fasta {input.sniffles_fasta:q} \
            --direct-fasta {input.direct_fasta:q} \
            --canonical-genome {output.canonical_genome:q} \
            --observed-genome {output.observed_genome:q} \
            --canonical-bed {output.canonical_bed:q} \
            --observed-bed {output.observed_bed:q} \
            --canonical-public-bed {output.canonical_public:q} \
            --observed-public-bed {output.observed_public:q} \
            --audit {output.audit:q} > {log:q} 2>&1
        """


rule index_outsider_integrated_genome:
    input:
        fasta=OUTSIDER_CANONICAL_GENOME,
    output:
        index=OUTSIDER_CANONICAL_INDEX,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_integrated_genome_index.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_integrated_genome_index.tsv",
    shell:
        "samtools faidx {input.fasta:q} > {log:q} 2>&1"


rule prepare_outsider_liftoff_flanks:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/prepare_outsider_liftoff.py"),
        positions=OUTSIDER_CANONICAL_BED,
        genome_index=OUTSIDER_CANONICAL_INDEX,
    output:
        gff=OUTSIDER_LIFT_INPUT_GFF,
        features=OUTSIDER_LIFT_FEATURES,
    params:
        flank=OUTSIDER_LIFT_FLANK,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_prepare_liftoff.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_prepare_liftoff.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_LIFT_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} \
            --positions {input.positions:q} \
            --genome-index {input.genome_index:q} \
            --output-gff {output.gff:q} \
            --feature-file {output.features:q} \
            --flank-size {params.flank} > {log:q} 2>&1
        """


rule lift_outsider_flanks_to_reference:
    input:
        reference=PREPARED_REFERENCE,
        reference_index=PREPARED_REFERENCE_INDEX,
        genome=OUTSIDER_CANONICAL_GENOME,
        genome_index=OUTSIDER_CANONICAL_INDEX,
        annotation=OUTSIDER_LIFT_INPUT_GFF,
        features=OUTSIDER_LIFT_FEATURES,
    output:
        gff=OUTSIDER_LIFT_GFF,
        unmapped=OUTSIDER_LIFT_UNMAPPED,
        intermediates=directory(OUTSIDER_LIFT_INTERMEDIATES),
    params:
        executable=LIFTOFF_EXECUTABLE,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_liftoff.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_liftoff.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_LIFT_DIR} {output.intermediates:q}
        if test -s {input.annotation:q}; then
            {params.executable:q} -p {threads} \
                -g {input.annotation:q} -f {input.features:q} \
                -o {output.gff:q} -u {output.unmapped:q} \
                -dir {output.intermediates:q} \
                {input.reference:q} {input.genome:q} > {log:q} 2>&1
        else
            printf '##gff-version 3\n' > {output.gff:q}
            : > {output.unmapped:q}
            printf 'No integrated insertion eligible for Liftoff.\n' > {log:q}
        fi
        test -f {output.gff:q}
        test -f {output.unmapped:q}
        """


rule summarize_outsider_liftoff:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/summarize_outsider_liftoff.py"),
        lifted=OUTSIDER_LIFT_GFF,
        insider=f"{INSIDER_TE_DIR}/INSERTION_TE_ON_REF.bed",
    output:
        good=OUTSIDER_LIFT_GOOD,
        bad=OUTSIDER_LIFT_BAD,
        mapped_ids=OUTSIDER_LIFT_IDS,
        public=OUTSIDER_ON_REFERENCE,
        combined=POSITION_TE_ON_REFERENCE,
        audit=OUTSIDER_LIFT_AUDIT,
    params:
        max_gap=OUTSIDER_LIFT_MAX_GAP,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_summarize_liftoff.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_summarize_liftoff.tsv",
    shell:
        """
        set -euo pipefail
        python3 {input.script:q} \
            --lifted-gff {input.lifted:q} \
            --insider-bed {input.insider:q} \
            --good-bed {output.good:q} \
            --bad-bed {output.bad:q} \
            --mapped-ids {output.mapped_ids:q} \
            --public-bed {output.public:q} \
            --combined-bed {output.combined:q} \
            --audit {output.audit:q} \
            --max-gap {params.max_gap} > {log:q} 2>&1
        """
