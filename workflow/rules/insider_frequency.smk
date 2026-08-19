"""Legacy-compatible INSIDER empty-site frequency estimation.

The historical ``FREQ_INSIDERv2`` rule mixed interval preparation, two BAM
subsets, CIGAR parsing, depth measurement, and table aggregation in one shell
block.  These rules expose the real dependencies while keeping the historical
single-sample frequency semantics for regression compatibility.
"""

INSIDER_FREQUENCY_DIR = f"{WORKDIR}/INSIDER/FREQ_INSIDER"
INSIDER_MAPPING_TO_REF_DIR = f"{WORKDIR}/OUTSIDER/MAPPING_TO_REF"
INSIDER_FREQUENCY_INSERTIONS = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/INSERTION_TE.bed"
)
INSIDER_FREQUENCY_OUTER_FLANKS = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/FLANK_TE.bed"
)
INSIDER_FREQUENCY_INNER_FLANKS = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/FLANK_TE_IN.bed"
)
INSIDER_FREQUENCY_DEPTH_FLANKS = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/DEPTH_FK.bed"
)
INSIDER_FREQUENCY_OUTER_BAM = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/FLANK_TE_OUTER.bam"
)
INSIDER_FREQUENCY_SPANNING_BAM = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/SPANNING_TE.bam"
)
INSIDER_FREQUENCY_EMPTY_COUNTS = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/DEL_NB.bed"
)
INSIDER_FREQUENCY_FLANK_DEPTH = (
    f"{INSIDER_MAPPING_TO_REF_DIR}/DEPTH_FK.txt"
)
INSIDER_FREQUENCY = f"{INSIDER_FREQUENCY_DIR}/DEPTH_TE_INSIDER.csv"


def legacy_insider_exclude_flag(options):
    """Keep the historical supplementary-alignment exclusion policy."""

    if re.search(r"(?:^|\s)-F\s+[0-9]+", options or ""):
        return ""
    return "-F 2048"


INSIDER_FREQUENCY_EXCLUDE_FLAG = legacy_insider_exclude_flag(
    SAMTOOLS_VIEW_OPTIONS
)


rule insider_frequency:
    input:
        INSIDER_FREQUENCY,


rule prepare_insider_frequency_intervals:
    input:
        script=str(
            PIPELINE_ROOT / "lib/python/workflow/prepare_insider_frequency.py"
        ),
        candidates=INSIDER_INSERTION_CSV,
    output:
        insertions=INSIDER_FREQUENCY_INSERTIONS,
        outer=INSIDER_FREQUENCY_OUTER_FLANKS,
        inner=INSIDER_FREQUENCY_INNER_FLANKS,
        depth=INSIDER_FREQUENCY_DEPTH_FLANKS,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/insider_prepare_frequency_intervals.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_prepare_frequency_intervals.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_MAPPING_TO_REF_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} {input.candidates:q} {output.insertions:q} \
            {output.outer:q} {output.inner:q} {output.depth:q} \
            > {log:q} 2>&1
        """


rule subset_insider_frequency_outer_bam:
    input:
        bam=MAPPING_BAM,
        bai=MAPPING_BAM + ".bai",
        regions=INSIDER_FREQUENCY_OUTER_FLANKS,
    output:
        bam=temp(INSIDER_FREQUENCY_OUTER_BAM),
        bai=temp(INSIDER_FREQUENCY_OUTER_BAM + ".bai"),
    params:
        view=SAMTOOLS_VIEW_OPTIONS,
        exclude=INSIDER_FREQUENCY_EXCLUDE_FLAG,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_subset_outer_flanks.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_subset_outer_flanks.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_MAPPING_TO_REF_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        if test -s {input.regions:q}; then
            samtools view {params.view} --threads {threads} -h -b \
                {input.bam:q} {params.exclude} -L {input.regions:q} \
                -o {output.bam:q} 2> {log:q}
        else
            samtools view {params.view} --threads {threads} -H -b \
                {input.bam:q} -o {output.bam:q} 2> {log:q}
        fi
        samtools index -@ {threads} {output.bam:q} {output.bai:q} \
            2>> {log:q}
        """


rule retain_insider_frequency_spanning_reads:
    input:
        bam=INSIDER_FREQUENCY_OUTER_BAM,
        bai=INSIDER_FREQUENCY_OUTER_BAM + ".bai",
        regions=INSIDER_FREQUENCY_INNER_FLANKS,
    output:
        bam=temp(INSIDER_FREQUENCY_SPANNING_BAM),
        bai=temp(INSIDER_FREQUENCY_SPANNING_BAM + ".bai"),
    params:
        view=SAMTOOLS_VIEW_OPTIONS,
        exclude=INSIDER_FREQUENCY_EXCLUDE_FLAG,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_retain_spanning_reads.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_retain_spanning_reads.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_MAPPING_TO_REF_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        if test -s {input.regions:q}; then
            samtools view {params.view} --threads {threads} -h -b \
                {input.bam:q} {params.exclude} -L {input.regions:q} \
                -o {output.bam:q} 2> {log:q}
        else
            samtools view {params.view} --threads {threads} -H -b \
                {input.bam:q} -o {output.bam:q} 2> {log:q}
        fi
        samtools index -@ {threads} {output.bam:q} {output.bai:q} \
            2>> {log:q}
        """


rule count_insider_empty_site_reads:
    input:
        script=str(
            PIPELINE_ROOT / "lib/python/workflow/count_insider_empty_sites.py"
        ),
        bam=INSIDER_FREQUENCY_SPANNING_BAM,
        bai=INSIDER_FREQUENCY_SPANNING_BAM + ".bai",
        candidates=INSIDER_FREQUENCY_INSERTIONS,
    output:
        counts=INSIDER_FREQUENCY_EMPTY_COUNTS,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_count_empty_sites.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_count_empty_sites.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_MAPPING_TO_REF_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} --size-window 30 --breakpoint-distance 30 \
            --threads {threads} {input.bam:q} {input.candidates:q} \
            {output.counts:q} > {log:q} 2>&1
        """


rule measure_insider_frequency_flank_depth:
    input:
        bam=INSIDER_FREQUENCY_SPANNING_BAM,
        bai=INSIDER_FREQUENCY_SPANNING_BAM + ".bai",
        regions=INSIDER_FREQUENCY_DEPTH_FLANKS,
    output:
        depth=INSIDER_FREQUENCY_FLANK_DEPTH,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_measure_flank_depth.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_measure_flank_depth.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_MAPPING_TO_REF_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        if test -s {input.regions:q}; then
            samtools depth {input.bam:q} -b {input.regions:q} 2> {log:q} \
                | awk 'BEGIN {{OFS="\t"}} {{print $1":"$2,$3}}' \
                > {output.depth:q}
        else
            : > {output.depth:q}
            : > {log:q}
        fi
        """


rule summarize_insider_frequency:
    input:
        script=str(
            PIPELINE_ROOT / "lib/python/workflow/summarize_insider_frequency.py"
        ),
        candidates=INSIDER_FREQUENCY_INSERTIONS,
        empty=INSIDER_FREQUENCY_EMPTY_COUNTS,
        depth=INSIDER_FREQUENCY_FLANK_DEPTH,
    output:
        frequency=INSIDER_FREQUENCY,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/insider_summarize_frequency.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_summarize_frequency.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_FREQUENCY_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} --depth-margin 30 {input.candidates:q} \
            {input.empty:q} {input.depth:q} {output.frequency:q} \
            > {log:q} 2>&1
        """
