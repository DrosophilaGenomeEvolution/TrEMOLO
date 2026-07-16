OUTSIDER_MAPPING_DIR = f"{WORKDIR}/OUTSIDER/MAPPING"
OUTSIDER_VARIANT_DIR = f"{WORKDIR}/OUTSIDER/VARIANT_CALLING"
MINIMAP2_INDEX = f"{OUTSIDER_MAPPING_DIR}/genome.mmi"
MAPPING_SAM = f"{OUTSIDER_MAPPING_DIR}/SAMPLE_mapping_GENOME.sam"
MAPPING_BAM = f"{OUTSIDER_MAPPING_DIR}/SAMPLE_mapping_GENOME_MD.sorted.bam"
MAPPING_STATS = f"{OUTSIDER_MAPPING_DIR}/stats.txt"
SV_VCF = f"{OUTSIDER_VARIANT_DIR}/SV.vcf"


def option_without_threads(value):
    return re.sub(r"(?:-t|--threads|-@)\s+\d+", "", value or "").strip()


MINIMAP2_PRESET = OUTSIDER_PARAMS.get("MINIMAP2", {}).get("PRESET_OPTION") or "map-ont"
MINIMAP2_OPTIONS = option_without_threads(
    OUTSIDER_PARAMS.get("MINIMAP2", {}).get("OPTION", "")
)
SAMTOOLS_VIEW_OPTIONS = option_without_threads(
    OUTSIDER_PARAMS.get("SAMTOOLS_VIEW", {}).get("PRESET_OPTION", "")
)
SAMTOOLS_SORT_OPTIONS = option_without_threads(
    OUTSIDER_PARAMS.get("SAMTOOLS_SORT", {}).get("PRESET_OPTION", "")
)
SAMTOOLS_CALMD_OPTIONS = option_without_threads(
    OUTSIDER_PARAMS.get("SAMTOOLS_CALLMD", {}).get("PRESET_OPTION", "")
)


rule outsider_variant_calling:
    input:
        SV_VCF,


rule minimap2_index_outsider:
    input:
        genome=PREPARED_GENOME,
    output:
        index=MINIMAP2_INDEX,
    params:
        preset=MINIMAP2_PRESET,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/minimap2_index.log",
    benchmark:
        f"{WORKDIR}/benchmarks/minimap2_index.tsv",
    shell:
        """
        mkdir -p {OUTSIDER_MAPPING_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        minimap2 -x {params.preset} -d {output.index:q} {input.genome:q} \
            > {log:q} 2>&1
        """


rule map_reads_outsider:
    input:
        genome=PREPARED_GENOME,
        reads=PREPARED_SAMPLE,
        index=MINIMAP2_INDEX,
    output:
        sam=MAPPING_SAM,
    params:
        preset=MINIMAP2_PRESET,
        options=MINIMAP2_OPTIONS,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/minimap2_mapping.log",
    benchmark:
        f"{WORKDIR}/benchmarks/minimap2_mapping.tsv",
    shell:
        """
        minimap2 -ax {params.preset} -t {threads} {params.options} \
            {input.index:q} {input.reads:q} > {output.sam:q} 2> {log:q}
        test -s {output.sam:q}
        """


rule prepare_bam_outsider:
    input:
        sam=MAPPING_SAM,
        genome=PREPARED_GENOME,
    output:
        bam=MAPPING_BAM,
        stats=MAPPING_STATS,
    params:
        view=SAMTOOLS_VIEW_OPTIONS,
        sort=SAMTOOLS_SORT_OPTIONS,
        calmd=SAMTOOLS_CALMD_OPTIONS,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/samtools.log",
    benchmark:
        f"{WORKDIR}/benchmarks/samtools.tsv",
    shell:
        """
        set -euo pipefail
        view_bam={output.bam:q}.view.bam
        sorted_bam={output.bam:q}.sorted.bam
        samtools view -h --threads {threads} {params.view} -F 4 -b \
            {input.sam:q} -o "$view_bam" 2> {log:q}
        samtools sort --threads {threads} {params.sort} "$view_bam" \
            -o "$sorted_bam" 2>> {log:q}
        samtools calmd --threads {threads} {params.calmd} -b "$sorted_bam" \
            {input.genome:q} > {output.bam:q} 2>> {log:q}
        samtools stats -@ {threads} {output.bam:q} \
            | grep '^SN' | cut -f 2-6 > {output.stats:q}
        """


rule call_sniffles_outsider:
    input:
        bam=MAPPING_BAM,
    output:
        vcf=SV_VCF,
        bai=MAPPING_BAM + ".bai",
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/sniffles.log",
    benchmark:
        f"{WORKDIR}/benchmarks/sniffles.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_VARIANT_DIR}
        samtools index -@ {threads} {input.bam:q} {output.bai:q} 2> {log:q}
        version=$(sniffles -h 2>&1 | awk '/Version/ {{print $2; exit}}')
        case "$version" in
            1.0.10)
                sniffles -t {threads} --report_seq -s 1 -m {input.bam:q} \
                    -v {output.vcf:q} -n -1 2>> {log:q}
                ;;
            1.0.11|1.0.12|1.0.12b)
                sniffles -t {threads} --report-seq -s 1 -m {input.bam:q} \
                    -v {output.vcf:q} -n -1 2>> {log:q}
                ;;
            *)
                echo "Unsupported Sniffles version: $version" >> {log:q}
                exit 2
                ;;
        esac
        test -s {output.vcf:q}
        """
