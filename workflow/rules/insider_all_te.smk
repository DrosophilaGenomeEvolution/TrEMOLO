INSIDER_ALL_TE_FIRST_BLAST = f"{INSIDER_TE_DIR}/ALL_TE_1.bln"
INSIDER_ALL_TE_FIRST_BED = f"{INSIDER_TE_DIR}/ALL_TE_1.bed"
INSIDER_ALL_TE_CLUSTERS = f"{INSIDER_TE_DIR}/ALL_TE_1_merge_cluster.bed"
INSIDER_ALL_TE_CLUSTER_FASTA = f"{INSIDER_TE_DIR}/ALL_TE_1_merge_cluster.fasta"
INSIDER_ALL_TE_SECOND_BLAST = f"{INSIDER_TE_DIR}/ALL_TE_2.bln"
INSIDER_ALL_TE_CSV = f"{INSIDER_TE_DIR}/ALL_TE.csv"
INSIDER_ALL_TE_COMBINE = f"{INSIDER_TE_DIR}/ALL_TE_COMBINE_TE.csv"
INSIDER_ALL_TE_BED = f"{INSIDER_TE_DIR}/POSITION_ALL_TE.bed"
INSIDER_ALL_TE_PUBLIC_BED = f"{WORKDIR}/POSITION_ALL_TE.bed"
INSIDER_ALL_TE_FORMATTER = (
    PIPELINE_ROOT / "lib/python/workflow/format_insider_all_te.py"
)
INSIDER_ALL_TE_CHROM_REGEX = str(
    OUTSIDER_PARAMS.get("TE_DETECTION", {}).get("CHROM_KEEP", ".")
).replace(",", "|")

try:
    re.compile(INSIDER_ALL_TE_CHROM_REGEX)
except re.error as error:
    raise ValueError(
        "PARAMS.OUTSIDER_VARIANT.TE_DETECTION.CHROM_KEEP is not a valid "
        f"regular expression: {error}"
    )


rule insider_all_te:
    input:
        csv=INSIDER_ALL_TE_CSV,
        combine=INSIDER_ALL_TE_COMBINE,
        bed=INSIDER_ALL_TE_BED,
        public_bed=INSIDER_ALL_TE_PUBLIC_BED,


rule blast_whole_assembly_against_te:
    input:
        genome=PREPARED_GENOME,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=temp(INSIDER_ALL_TE_FIRST_BLAST),
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_all_te_first_blast.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_all_te_first_blast.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {INSIDER_TE_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        blastn -num_threads {threads} -db {input.database:q} \
            -query {input.genome:q} -outfmt 6 -out {output.blast:q} \
            > {log:q} 2>&1
        """


rule prepare_whole_assembly_te_regions:
    input:
        blast=INSIDER_ALL_TE_FIRST_BLAST,
        genome=PREPARED_GENOME,
        genome_index=PREPARED_GENOME_INDEX,
    output:
        bed=temp(INSIDER_ALL_TE_FIRST_BED),
        clusters=temp(INSIDER_ALL_TE_CLUSTERS),
        fasta=temp(INSIDER_ALL_TE_CLUSTER_FASTA),
    params:
        chrom_regex=INSIDER_ALL_TE_CHROM_REGEX,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_all_te_regions.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_all_te_regions.tsv",
    shell:
        r"""
        set -euo pipefail
        : > {log:q}
        : > {output.bed:q}
        : > {output.clusters:q}
        : > {output.fasta:q}
        awk -F '\t' -v pattern={params.chrom_regex:q} 'BEGIN {{OFS="\t"}}
            NF == 0 {{next}}
            NF != 12 {{
                print "Malformed first-pass BLAST row at line " NR > "/dev/stderr"
                exit 2
            }}
            {{candidate=$1 OFS $7 OFS $8}}
            candidate ~ pattern {{print $1,$7,$8}}
        ' {input.blast:q} > {output.bed:q} 2>> {log:q}
        if test -s {output.bed:q}; then
            bedtools sort -i {output.bed:q} \
                | bedtools merge -d 100 \
                | bedtools cluster > {output.clusters:q} 2>> {log:q}
            bedtools getfasta -fi {input.genome:q} -bed {output.clusters:q} \
                -name+ > {output.fasta:q} 2>> {log:q}
            expected=$(awk 'END {{print NR+0}}' {output.clusters:q})
            observed=$(awk '/^>/ {{count++}} END {{print count+0}}' {output.fasta:q})
            if test "$expected" -ne "$observed"; then
                printf 'Incomplete whole-assembly FASTA extraction: expected %s records, got %s.\n' \
                    "$expected" "$observed" >> {log:q}
                exit 1
            fi
        fi
        """


rule blast_whole_assembly_te_regions:
    input:
        query=INSIDER_ALL_TE_CLUSTER_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=temp(INSIDER_ALL_TE_SECOND_BLAST),
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_all_te_second_blast.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_all_te_second_blast.tsv",
    shell:
        """
        set -euo pipefail
        : > {output.blast:q}
        : > {log:q}
        if test -s {input.query:q}; then
            blastn -num_threads {threads} -db {input.database:q} \
                -query {input.query:q} -outfmt 6 -out {output.blast:q} \
                >> {log:q} 2>&1
        else
            printf 'No whole-assembly TE region; second BLAST skipped.\n' >> {log:q}
        fi
        """


rule classify_insider_all_te:
    input:
        script=str(PIPELINE_ROOT / "lib/python/parsing/global_sv.py"),
        blast=INSIDER_ALL_TE_SECOND_BLAST,
        database=PREPARED_TE_DATABASE,
    output:
        csv=INSIDER_ALL_TE_CSV,
        combine=INSIDER_ALL_TE_COMBINE,
    params:
        options=INSIDER_BLAST_FILTER,
        temporary_csv=INSIDER_ALL_TE_CSV + ".tmp",
        temporary_combine=INSIDER_ALL_TE_COMBINE + ".tmp",
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_all_te_classify.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_all_te_classify.tsv",
    shell:
        """
        set -euo pipefail
        rm -f {params.temporary_csv:q} {params.temporary_combine:q}
        cleanup_insider_all_te_classification() {{
            rm -f {params.temporary_csv:q} {params.temporary_combine:q}
        }}
        trap cleanup_insider_all_te_classification EXIT
        python3 {input.script:q} {params.options} \
            --combine_name {params.temporary_combine:q} \
            {input.blast:q} {input.database:q} {params.temporary_csv:q} \
            > {log:q} 2>&1
        mv -f {params.temporary_csv:q} {output.csv:q}
        mv -f {params.temporary_combine:q} {output.combine:q}
        trap - EXIT
        """


rule format_insider_all_te:
    input:
        script=str(INSIDER_ALL_TE_FORMATTER),
        csv=INSIDER_ALL_TE_CSV,
    output:
        bed=INSIDER_ALL_TE_BED,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/insider_all_te_format.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_all_te_format.tsv",
    shell:
        "python3 {input.script:q} {input.csv:q} {output.bed:q} > {log:q} 2>&1"


rule publish_insider_all_te:
    input:
        bed=INSIDER_ALL_TE_BED,
    output:
        public_bed=INSIDER_ALL_TE_PUBLIC_BED,
    threads: 1
    resources:
        mem_mb=128,
    shell:
        "ln -s -r -f {input.bed:q} {output.public_bed:q}"
