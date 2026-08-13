"""OUTSIDER flank extraction and target-site duplication rules.

The legacy ``GET_SEQ_TE`` rule mixed candidate selection, FASTA indexing,
sequence extraction, formatting, and cleanup in one shell block.  These rules
keep its historical candidate-selection semantics while declaring every input
and the complete flank directories as Snakemake outputs.
"""

OUTSIDER_ET_FASTA_DIR = f"{WORKDIR}/OUTSIDER/ET_FIND_FA"
OUTSIDER_GENOME_INDEX = PREPARED_GENOME + ".fai"
OUTSIDER_FK_DIR = f"{WORKDIR}/OUTSIDER/FK"
OUTSIDER_FK_READS_DIR = f"{OUTSIDER_FK_DIR}/READS"
OUTSIDER_FK_GENOME_DIR = f"{OUTSIDER_FK_DIR}/GENOME"
OUTSIDER_TSD_DIR = f"{WORKDIR}/OUTSIDER/TSD"

OUTSIDER_TSD_CANDIDATES = f"{OUTSIDER_FK_DIR}/TSD_CANDIDATES.tsv"
OUTSIDER_TSD_FLANK_ELIGIBILITY = f"{OUTSIDER_FK_DIR}/TSD_FLANK_ELIGIBILITY.tsv"
OUTSIDER_SEQUENCE_SIZES = f"{WORKDIR}/OUTSIDER/SIZE_SEQ.tsv"
OUTSIDER_FK_READS = f"{OUTSIDER_FK_DIR}/ALL_FK_REPORT_FT_READS.bed"
OUTSIDER_FK_GENOME = f"{OUTSIDER_FK_DIR}/ALL_FK_REPORT_FT_GENOME.bed"
OUTSIDER_FK_ALL = f"{OUTSIDER_FK_DIR}/ALL_FK_REPORT_FT.bed"
OUTSIDER_TSD = f"{OUTSIDER_TSD_DIR}/TSD_TE.tsv"


rule outsider_tsd:
    input:
        OUTSIDER_TSD,


rule outsider_tsd_flanks:
    input:
        OUTSIDER_FK_ALL,
        OUTSIDER_TSD_FLANK_ELIGIBILITY,


rule index_outsider_tsd_genome:
    input:
        fasta=PREPARED_GENOME,
    output:
        fai=OUTSIDER_GENOME_INDEX,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_genome_fasta_index.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_genome_fasta_index.tsv",
    shell:
        "samtools faidx {input.fasta:q} > {log:q} 2>&1"


rule prepare_outsider_tsd_flanks:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/prepare_outsider_tsd.py"),
        genome=PREPARED_GENOME,
        genome_index=OUTSIDER_GENOME_INDEX,
        sniffles_calls=OUTSIDER_SNIFFLES_CSV,
        sniffles_combined=OUTSIDER_SNIFFLES_COMBINE,
        sniffles_fasta=OUTSIDER_SNIFFLES_FASTA,
        merged=OUTSIDER_MERGED_BED,
        direct_calls=OUTSIDER_INS_CSV,
        direct_combined=OUTSIDER_INS_COMBINE,
        direct_fasta=OUTSIDER_INS_FASTA,
        sequence_sizes=OUTSIDER_SV_SIZE,
    output:
        te_fastas=directory(OUTSIDER_ET_FASTA_DIR),
        read_artifacts=directory(OUTSIDER_FK_READS_DIR),
        genome_artifacts=directory(OUTSIDER_FK_GENOME_DIR),
        candidates=OUTSIDER_TSD_CANDIDATES,
        validation=OUTSIDER_TSD_FLANK_ELIGIBILITY,
        sizes=OUTSIDER_SEQUENCE_SIZES,
        reads=OUTSIDER_FK_READS,
        genome=OUTSIDER_FK_GENOME,
        combined=OUTSIDER_FK_ALL,
    params:
        flank=OUTSIDER_SIZE_FLANK,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_prepare_tsd_flanks.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_prepare_tsd_flanks.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} \
            --genome {input.genome:q} \
            --genome-index {input.genome_index:q} \
            --sniffles-calls {input.sniffles_calls:q} \
            --sniffles-combined {input.sniffles_combined:q} \
            --sniffles-fasta {input.sniffles_fasta:q} \
            --merged-bed {input.merged:q} \
            --direct-calls {input.direct_calls:q} \
            --direct-combined {input.direct_combined:q} \
            --direct-fasta {input.direct_fasta:q} \
            --sequence-sizes {input.sequence_sizes:q} \
            --flank-size {params.flank} \
            --compatibility-mode legacy \
            --te-fasta-dir {output.te_fastas:q} \
            --read-artifact-dir {output.read_artifacts:q} \
            --genome-artifact-dir {output.genome_artifacts:q} \
            --candidates {output.candidates:q} \
            --validation {output.validation:q} \
            --all-sequence-sizes {output.sizes:q} \
            --read-flanks {output.reads:q} \
            --genome-flanks {output.genome:q} \
            --combined-flanks {output.combined:q} \
            > {log:q} 2>&1
        """


rule call_outsider_tsd:
    input:
        script=str(PIPELINE_ROOT / "lib/python/TSD/getTSDreads.py"),
        utils=str(PIPELINE_ROOT / "lib/python/TSD/utils.py"),
        flanks=OUTSIDER_FK_ALL,
    output:
        tsd=OUTSIDER_TSD,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_call_tsd.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_call_tsd.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_TSD_DIR}
        python3 {input.script:q} -t {threads} \
            {input.flanks:q} {output.tsd:q} > {log:q} 2>&1
        """
