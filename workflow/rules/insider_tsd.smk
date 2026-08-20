"""Legacy-compatible INSIDER target-site duplication calling."""

INSIDER_TSD_DIR = f"{WORKDIR}/INSIDER/TSD"
INSIDER_TSD_MERGED_BED = f"{INSIDER_TSD_DIR}/tmp_merge_ins_del.bed"
INSIDER_TSD_FLANK_BED = f"{INSIDER_TSD_DIR}/FK_TE.bed"
INSIDER_TSD_FLANKS = f"{INSIDER_TSD_DIR}/FK_TE_FT.bed"
INSIDER_TSD = f"{INSIDER_TSD_DIR}/TSD_TE.tsv"


rule insider_tsd:
    input:
        INSIDER_TSD,


rule prepare_insider_tsd_flanks:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/prepare_insider_tsd.py"),
        genome=PREPARED_GENOME,
        genome_index=PREPARED_GENOME_INDEX,
        insertion=INSIDER_INSERTION_TE_BED,
        deletion=INSIDER_DELETION_TE_BED,
    output:
        merged=INSIDER_TSD_MERGED_BED,
        flank_bed=INSIDER_TSD_FLANK_BED,
        flanks=INSIDER_TSD_FLANKS,
    params:
        flank=OUTSIDER_SIZE_FLANK,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_prepare_tsd_flanks.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_prepare_tsd_flanks.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_TSD_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} \
            --genome {input.genome:q} \
            --genome-index {input.genome_index:q} \
            --insertion-bed {input.insertion:q} \
            --deletion-bed {input.deletion:q} \
            --flank-size {params.flank} \
            --merged-bed {output.merged:q} \
            --flank-bed {output.flank_bed:q} \
            --formatted-flanks {output.flanks:q} \
            > {log:q} 2>&1
        """


rule call_insider_tsd:
    input:
        script=str(PIPELINE_ROOT / "lib/python/TSD/getTSDgenome.py"),
        utils=str(PIPELINE_ROOT / "lib/python/TSD/utils.py"),
        flanks=INSIDER_TSD_FLANKS,
    output:
        tsd=INSIDER_TSD,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_call_tsd.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_call_tsd.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_TSD_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} -t {threads} \
            {input.flanks:q} {output.tsd:q} > {log:q} 2>&1
        """
