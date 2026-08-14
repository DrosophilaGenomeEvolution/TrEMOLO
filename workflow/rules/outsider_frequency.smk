"""Legacy-compatible OUTSIDER allele-frequency estimation.

The historical ``FREQUENCEv2`` rule prepared candidates, sliced the BAM,
started untracked background workers, and aggregated their results in one
shell block.  This module exposes each data dependency and keeps the legacy
frequency semantics isolated until a separately validated population model is
introduced.
"""

OUTSIDER_FREQUENCY_DIR = f"{WORKDIR}/OUTSIDER/FREQUENCY"
OUTSIDER_FREQUENCY_CALLS = f"{OUTSIDER_FREQUENCY_DIR}/FILTER_BLAST_INS.csv"
OUTSIDER_FREQUENCY_CALLS_SORTED = (
    f"{OUTSIDER_FREQUENCY_DIR}/FILTER_BLAST_INS.sorted.csv"
)
OUTSIDER_FREQUENCY_SV_SIZES = f"{OUTSIDER_FREQUENCY_DIR}/SV_SIZE.tsv"
OUTSIDER_FREQUENCY_REGIONS = f"{OUTSIDER_FREQUENCY_DIR}/POSITION_START_TE.bed"
OUTSIDER_FREQUENCY_BAM = f"{OUTSIDER_FREQUENCY_DIR}/MAPPING_POSTION_TE.bam"
OUTSIDER_FREQUENCY_BAI = OUTSIDER_FREQUENCY_BAM + ".bai"
OUTSIDER_FREQUENCY_COUNTS_RAW = (
    f"{OUTSIDER_FREQUENCY_DIR}/COUNT_READS.unsorted.txt"
)
OUTSIDER_FREQUENCY_COUNTS = f"{OUTSIDER_FREQUENCY_DIR}/COUNT_READS.txt"
OUTSIDER_FREQUENCY_COUNTS_BY_EVENT = (
    f"{OUTSIDER_FREQUENCY_DIR}/COUNT_READS.by_event.txt"
)
OUTSIDER_FREQUENCY_READ_SUPPORT = (
    f"{OUTSIDER_FREQUENCY_DIR}/COUNT_TE_IN_RS.txt"
)
OUTSIDER_FREQUENCY = f"{OUTSIDER_FREQUENCY_DIR}/FREQUENCY_TE_INS.tsv"
OUTSIDER_FREQUENCY_PRECISE = (
    f"{OUTSIDER_FREQUENCY_DIR}/FREQUENCY_TE_INS_PRECISE.tsv"
)
OUTSIDER_FREQUENCY_HEADER = "\t".join(
    (
        "sseqid", "qseqid", "pident", "size_per", "size_el", "mismatch",
        "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    )
)


def legacy_frequency_size_percent(options):
    """Derive the insertion threshold used by FREQUENCEv2."""

    match = re.search(r"(?:^|\s)--min-size-percent(?:=|\s+)(\d+)", options or "")
    if not match:
        raise ValueError(
            "PARAMS.OUTSIDER_VARIANT.PARS_BLN_OPTION must define "
            "--min-size-percent for OUTSIDER frequency estimation"
        )
    return max(int(match.group(1)) - 10, 5)


rule outsider_frequency:
    input:
        OUTSIDER_FREQUENCY,
        OUTSIDER_FREQUENCY_PRECISE,


rule combine_outsider_frequency_calls:
    input:
        direct=OUTSIDER_INS_CSV,
        sniffles=OUTSIDER_SNIFFLES_CSV,
        soft=OUTSIDER_SOFT_CSV,
        hard=OUTSIDER_HARD_CSV,
    output:
        calls=OUTSIDER_FREQUENCY_CALLS,
    params:
        header=OUTSIDER_FREQUENCY_HEADER,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_combine_frequency_calls.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_combine_frequency_calls.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_FREQUENCY_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        : > {output.calls:q}
        for source in {input.direct:q} {input.sniffles:q} {input.soft:q} {input.hard:q}; do
            if test -s "$source"; then
                head -n 1 "$source" > {output.calls:q}
                break
            fi
        done
        if ! test -s {output.calls:q}; then
            printf '%s\n' {params.header:q} > {output.calls:q}
        fi
        for source in {input.direct:q} {input.sniffles:q} {input.soft:q} {input.hard:q}; do
            awk 'NR>1' "$source" >> {output.calls:q}
        done
        awk 'END {{print NR-1 " frequency candidates"}}' {output.calls:q} > {log:q}
        """


rule sort_outsider_frequency_calls:
    input:
        calls=OUTSIDER_FREQUENCY_CALLS,
    output:
        calls=OUTSIDER_FREQUENCY_CALLS_SORTED,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_sort_frequency_calls.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_sort_frequency_calls.tsv",
    shell:
        """
        set -euo pipefail
        head -n 1 {input.calls:q} > {output.calls:q}
        awk 'NR>1 && OFS="\t" {{split($2,position,":"); print position[1],position[3],position[4]+1,$0}}' \
            {input.calls:q} | bedtools sort | cut -f 4- >> {output.calls:q} \
            2> {log:q}
        """


rule measure_outsider_frequency_variants:
    input:
        sniffles=OUTSIDER_SNIFFLES_FASTA,
        direct=OUTSIDER_INS_FASTA,
        soft=OUTSIDER_SOFT_FASTA,
        hard=OUTSIDER_HARD_FASTA,
    output:
        sizes=OUTSIDER_FREQUENCY_SV_SIZES,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_measure_frequency_variants.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_measure_frequency_variants.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_FREQUENCY_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        : > {output.sizes:q}
        awk '/^>/ {{header=substr($0,2)}} /^[^>]/ && OFS="\t" {{print header,length($0)}}' \
            {input.sniffles:q} >> {output.sizes:q}
        awk '/^>/ {{header=substr($0,2)}} /^[^>]/ && OFS="\t" {{print header,length($0)}}' \
            {input.direct:q} >> {output.sizes:q}
        awk '/^>/ {{header=substr($0,2)}} /^[^>]/ && OFS="\t" {{print header,length($0)}}' \
            {input.soft:q} >> {output.sizes:q}
        awk '/^>/ {{header=substr($0,2)}} /^[^>]/ && OFS="\t" {{print header,length($0)}}' \
            {input.hard:q} >> {output.sizes:q}
        printf 'Measured variant sequences from four OUTSIDER sources.\n' > {log:q}
        """


rule build_outsider_frequency_regions:
    input:
        calls=OUTSIDER_FREQUENCY_CALLS_SORTED,
    output:
        regions=OUTSIDER_FREQUENCY_REGIONS,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/outsider_build_frequency_regions.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_build_frequency_regions.tsv",
    shell:
        """
        set -euo pipefail
        awk -v window=50 'NR>1 && OFS="\t" {{split($2,position,":"); print position[1],(position[3]-window>0 ? position[3]-window : 1),position[4]+window}}' \
            {input.calls:q} | bedtools merge -d 1 > {output.regions:q} \
            2> {log:q}
        """


rule subset_outsider_frequency_bam:
    input:
        bam=MAPPING_BAM,
        bai=MAPPING_BAM + ".bai",
        regions=OUTSIDER_FREQUENCY_REGIONS,
    output:
        bam=OUTSIDER_FREQUENCY_BAM,
        bai=OUTSIDER_FREQUENCY_BAI,
    params:
        view=SAMTOOLS_VIEW_OPTIONS,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_subset_frequency_bam.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_subset_frequency_bam.tsv",
    shell:
        """
        set -euo pipefail
        if test -s {input.regions:q}; then
            samtools view {params.view} --threads {threads} -q 10 -F 4 -h -b \
                {input.bam:q} -L {input.regions:q} > {output.bam:q} 2> {log:q}
        else
            samtools view {params.view} --threads {threads} -H -b \
                {input.bam:q} > {output.bam:q} 2> {log:q}
        fi
        samtools index -@ {threads} {output.bam:q} {output.bai:q} \
            2>> {log:q}
        """


rule count_outsider_frequency_reads:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/count_outsider_frequency.py"),
        bam=OUTSIDER_FREQUENCY_BAM,
        bai=OUTSIDER_FREQUENCY_BAI,
        calls=OUTSIDER_FREQUENCY_CALLS_SORTED,
        te_sizes=TE_SIZE_TABLE,
        sv_sizes=OUTSIDER_FREQUENCY_SV_SIZES,
    output:
        counts=temp(OUTSIDER_FREQUENCY_COUNTS_RAW),
    params:
        size_percent=lambda wildcards: legacy_frequency_size_percent(
            OUTSIDER_BLAST_FILTER
        ),
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_count_frequency_reads.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_count_frequency_reads.tsv",
    shell:
        """
        set -euo pipefail
        python3 {input.script:q} --per-size {params.size_percent} \
            --threads {threads} {input.bam:q} {input.calls:q} \
            {input.te_sizes:q} {input.sv_sizes:q} {output.counts:q} \
            > {log:q} 2>&1
        """


rule sort_outsider_frequency_read_counts:
    input:
        counts=OUTSIDER_FREQUENCY_COUNTS_RAW,
    output:
        counts=OUTSIDER_FREQUENCY_COUNTS,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_sort_frequency_read_counts.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_sort_frequency_read_counts.tsv",
    shell:
        """
        set -euo pipefail
        LC_ALL=C sort -k 3 {input.counts:q} > {output.counts:q} \
            2> {log:q}
        """


rule order_outsider_frequency_counts_by_event:
    input:
        counts=OUTSIDER_FREQUENCY_COUNTS,
    output:
        counts=temp(OUTSIDER_FREQUENCY_COUNTS_BY_EVENT),
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_order_frequency_counts_by_event.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_order_frequency_counts_by_event.tsv",
    shell:
        """
        set -euo pipefail
        awk '{{print $3":"$2":"$1, $0}}' {input.counts:q} \
            | LC_ALL=C sort | cut -d ' ' -f 2- > {output.counts:q} \
            2> {log:q}
        """


rule summarize_outsider_frequency:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/summarize_outsider_frequency.py"),
        counts=OUTSIDER_FREQUENCY_COUNTS_BY_EVENT,
        sniffles_support=OUTSIDER_SNIFFLES_TE_READ_COUNTS,
        direct_support=OUTSIDER_INS_READ_COUNTS,
    output:
        support=OUTSIDER_FREQUENCY_READ_SUPPORT,
        frequency=OUTSIDER_FREQUENCY,
        precise=OUTSIDER_FREQUENCY_PRECISE,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_summarize_frequency.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_summarize_frequency.tsv",
    shell:
        """
        set -euo pipefail
        python3 {input.script:q} {input.counts:q} \
            {input.sniffles_support:q} {input.direct_support:q} \
            {output.frequency:q} {output.precise:q} \
            --combined-support-output {output.support:q} --pre-sorted \
            > {log:q} 2>&1
        """
