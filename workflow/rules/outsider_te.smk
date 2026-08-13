"""OUTSIDER transposable-element detection rules.

The legacy workflow mixed four evidence sources in a few monolithic rules.
This file keeps their historical formats while exposing each transformation in
the DAG.  Scientific changes can therefore be introduced later behind explicit
regression tests.
"""

OUTSIDER_TREMOLO_DIR = f"{WORKDIR}/OUTSIDER/TrEMOLO_SV_TE"
OUTSIDER_INS_DIR = f"{OUTSIDER_TREMOLO_DIR}/INS"
OUTSIDER_SOFT_DIR = f"{OUTSIDER_TREMOLO_DIR}/SOFT"
OUTSIDER_HARD_DIR = f"{OUTSIDER_TREMOLO_DIR}/HARD"
OUTSIDER_TE_DIR = f"{WORKDIR}/OUTSIDER/TE_DETECTION"
OUTSIDER_MERGE_DIR = f"{OUTSIDER_TE_DIR}/MERGE_TE"

OUTSIDER_INS_RAW = f"{OUTSIDER_INS_DIR}/SV_INS.bed"
OUTSIDER_INS_CLUSTERED = f"{OUTSIDER_INS_DIR}/SV_INS_CLUST.bed"
OUTSIDER_INS_FASTA = f"{OUTSIDER_INS_DIR}/SV_INS_CLUST.fasta"
OUTSIDER_INS_BLAST = f"{OUTSIDER_INS_DIR}/SV_INS_CLUST.bln"
OUTSIDER_INS_CSV = f"{OUTSIDER_INS_DIR}/INS_TREMOLO.csv"
OUTSIDER_INS_COMBINE = f"{OUTSIDER_INS_DIR}/COMBINE_INS_TREMOLO.csv"
OUTSIDER_INS_COUNT = f"{OUTSIDER_INS_DIR}/INS_TREMOLO_COUNT.csv"
OUTSIDER_INS_BED = f"{OUTSIDER_INS_DIR}/INS_TREMOLO.bed"
OUTSIDER_INS_IDS = f"{OUTSIDER_INS_DIR}/ID.txt"
OUTSIDER_INS_READ_COUNTS = f"{OUTSIDER_INS_DIR}/COUNT_TE_IN_RS.txt"
OUTSIDER_INS_TSD = f"{OUTSIDER_INS_DIR}/INS_FOR_TSD.txt"
OUTSIDER_INS_READ_POSITIONS = f"{OUTSIDER_INS_DIR}/RD_NUMBER.txt"
OUTSIDER_INS_SIZE_BASE = f"{OUTSIDER_INS_DIR}/SV_SIZE.direct.tsv"
OUTSIDER_SV_SIZE = f"{OUTSIDER_INS_DIR}/SV_SIZE.tsv"

OUTSIDER_SOFT_VCF = f"{OUTSIDER_SOFT_DIR}/SV_SOFT.vcf"
OUTSIDER_SOFT_VCF_ALL = OUTSIDER_SOFT_VCF + ".bis"
OUTSIDER_HARD_VCF = f"{OUTSIDER_HARD_DIR}/SV_HARD.tr_vcf"
OUTSIDER_SOFT_SV_BED = f"{OUTSIDER_SOFT_DIR}/SV_SOFT.bed"
OUTSIDER_SOFT_FASTA = f"{OUTSIDER_SOFT_DIR}/SV_SOFT.fasta"
OUTSIDER_SOFT_ALL_FASTA = OUTSIDER_SOFT_FASTA + ".bis"
OUTSIDER_SOFT_BLAST = f"{OUTSIDER_SOFT_DIR}/SV_SOFT.bln"
OUTSIDER_SOFT_ALL_BLAST = OUTSIDER_SOFT_BLAST + ".bis"
OUTSIDER_SOFT_CSV = f"{OUTSIDER_SOFT_DIR}/SOFT_TE.csv"
OUTSIDER_SOFT_COMBINE = f"{OUTSIDER_SOFT_DIR}/COMBINE_SV_SOFT.csv"
OUTSIDER_SOFT_COUNT = f"{OUTSIDER_SOFT_DIR}/SOFT_TE_COUNT.csv"
OUTSIDER_SOFT_BED = f"{OUTSIDER_SOFT_DIR}/SOFT_TE.bed"
OUTSIDER_SOFT_MATCH_BED = f"{OUTSIDER_SOFT_DIR}/TE_SOFT.bed"

OUTSIDER_READS_INDEX = PREPARED_SAMPLE + ".fai"
OUTSIDER_HARD_IDS = f"{OUTSIDER_HARD_DIR}/ID.txt"
OUTSIDER_HARD_FASTQ = f"{OUTSIDER_HARD_DIR}/RS_HARD.fastq"
OUTSIDER_HARD_READS_FASTA = f"{OUTSIDER_HARD_DIR}/RS_HARD.fasta"
OUTSIDER_HARD_READS_FASTA_INDEX = OUTSIDER_HARD_READS_FASTA + ".fa_trml_idx"
OUTSIDER_HARD_RAW_BED = f"{OUTSIDER_HARD_DIR}/HARD.bed"
OUTSIDER_HARD_FASTA = f"{OUTSIDER_HARD_DIR}/HARD.fasta"
OUTSIDER_HARD_BLAST = f"{OUTSIDER_HARD_DIR}/HARD.bln"
OUTSIDER_HARD_CSV = f"{OUTSIDER_HARD_DIR}/HARD_TE.csv"
OUTSIDER_HARD_COMBINE = f"{OUTSIDER_HARD_DIR}/COMBINE_SV_HARD.csv"
OUTSIDER_HARD_COUNT = f"{OUTSIDER_HARD_DIR}/HARD_TE_COUNT.csv"
OUTSIDER_HARD_BED = f"{OUTSIDER_HARD_DIR}/HARD_TE.bed"

OUTSIDER_SNIFFLES_RAW_FASTA = f"{OUTSIDER_VARIANT_DIR}/SEQUENCE_INDEL.raw.fasta"
OUTSIDER_SNIFFLES_FASTA = f"{OUTSIDER_VARIANT_DIR}/SEQUENCE_INDEL.fasta"
OUTSIDER_SNIFFLES_VARIANTS = f"{OUTSIDER_TE_DIR}/TE_VR.bed"
OUTSIDER_SNIFFLES_READ_COUNTS = f"{OUTSIDER_TE_DIR}/RD_NUMBER.txt"
OUTSIDER_SNIFFLES_BLAST = f"{OUTSIDER_TE_DIR}/BLAST_SEQUENCE_INDEL_vs_DBTE.bln"
OUTSIDER_SNIFFLES_CSV = f"{OUTSIDER_TE_DIR}/FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv"
OUTSIDER_SNIFFLES_COMBINE = f"{OUTSIDER_TE_DIR}/COMBINE_TE.csv"
OUTSIDER_SNIFFLES_COUNT = f"{OUTSIDER_TE_DIR}/FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE_COUNT.csv"
OUTSIDER_SNIFFLES_BED = f"{OUTSIDER_TE_DIR}/POSITION_START_TE.bed"
OUTSIDER_SNIFFLES_IDS = f"{OUTSIDER_TE_DIR}/ID.txt"
OUTSIDER_SNIFFLES_TE_READ_COUNTS = f"{OUTSIDER_TE_DIR}/COUNT_TE_IN_RS.txt"

OUTSIDER_MERGED_BED = f"{OUTSIDER_MERGE_DIR}/MERGE_TE_ALL.bed"
OUTSIDER_MERGED_BED_COPY = f"{OUTSIDER_MERGE_DIR}/MERGE_TE_ALL2.bed"
OUTSIDER_MERGED_COUNT = f"{OUTSIDER_MERGE_DIR}/MERGE_TE_ALL_COUNT.csv"
OUTSIDER_POSITION_BED = f"{WORKDIR}/POSITION_TE_OUTSIDER.bed"

OUTSIDER_TE_PARAMS = OUTSIDER_PARAMS.get("TE_DETECTION", {})
OUTSIDER_SIZE_FLANK = OUTSIDER_PARAMS.get("TSD", {}).get("SIZE_FLANK", 10)
OUTSIDER_TIME_LIMIT = OUTSIDER_TE_PARAMS.get("TIME_LIMIT", 0)
OUTSIDER_TIME_LIMIT = (
    OUTSIDER_TIME_LIMIT
    if isinstance(OUTSIDER_TIME_LIMIT, int) and OUTSIDER_TIME_LIMIT > 0
    else 0
)
OUTSIDER_CHROM_KEEP = OUTSIDER_TE_PARAMS.get("CHROM_KEEP", ".")
OUTSIDER_GET_SEQ_OPTIONS = OUTSIDER_TE_PARAMS.get("GET_SEQ_REPORT_OPTION", "-m 1000")
OUTSIDER_BLAST_FILTER = OUTSIDER_PARAMS.get(
    "PARS_BLN_OPTION", "--min-size-percent 90 --min-pident 90 -k 'INS|DEL'"
)
OUTSIDER_CLIPPED_FLAG = (
    "" if OUTSIDER_CHOICES.get("CLIPPED_READS", False) else "--no-clipped"
)


rule outsider_primary_te_detection:
    input:
        OUTSIDER_INS_BED,
        OUTSIDER_SNIFFLES_BED,


rule outsider_te_detection:
    input:
        OUTSIDER_POSITION_BED,
        OUTSIDER_MERGED_COUNT,


rule find_outsider_alignment_candidates:
    input:
        bam=MAPPING_BAM,
        bai=MAPPING_BAM + ".bai",
    output:
        insertions=OUTSIDER_INS_RAW,
        tsd=OUTSIDER_INS_TSD,
        soft=OUTSIDER_SOFT_VCF,
        soft_all=OUTSIDER_SOFT_VCF_ALL,
        hard=OUTSIDER_HARD_VCF,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/find_all_type_ins.py"),
        clipped=OUTSIDER_CLIPPED_FLAG,
        flank=OUTSIDER_SIZE_FLANK,
        time_limit=OUTSIDER_TIME_LIMIT,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_find_alignment_candidates.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_find_alignment_candidates.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_INS_DIR} {OUTSIDER_SOFT_DIR} {OUTSIDER_HARD_DIR}
        python3 {params.script:q} {input.bam:q} \
            --output-seq-tsd {output.tsd:q} \
            --flank-size {params.flank} \
            --threads {threads} \
            --time-limit {params.time_limit} \
            --output-soft {output.soft:q} \
            --output-hard {output.hard:q} \
            --output-ins {output.insertions:q} {params.clipped} \
            > {log:q} 2>&1
        test -s {output.insertions:q}
        touch {output.soft:q} {output.soft_all:q} {output.hard:q}
        """


rule cluster_outsider_alignment_insertions:
    input:
        insertions=OUTSIDER_INS_RAW,
    output:
        clustered=OUTSIDER_INS_CLUSTERED,
        fasta=OUTSIDER_INS_FASTA,
        read_positions=OUTSIDER_INS_READ_POSITIONS,
        sizes=OUTSIDER_INS_SIZE_BASE,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_cluster_insertions.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_cluster_insertions.tsv",
    shell:
        """
        set -euo pipefail
        bedtools sort -i {input.insertions:q} | bedtools cluster -d 50 \
            > {output.clustered:q} 2> {log:q}
        : > {output.read_positions:q}
        awk -v FILE_RD_NUM={output.read_positions:q} '
            BEGIN {{clust = 0}}
            {{
                if (clust == 0) {{
                    chrom=$1; start=$2; end=$3; clust=$NF;
                    dic_clust[clust][1]["header"]=">"chrom":<INS>:"start":"end":TrEMOLO.INS."clust;
                    dic_clust[clust][1]["seq"]=$5;
                    dic_clust[clust][1]["RS"]=1;
                    dic_clust[clust][1]["READ"]=$4;
                }} else if (clust != $NF) {{
                    chrom=$1; start=$2; end=$3;
                    for (i=1; i<dic_clust[clust][1]["RS"]+1; i++) {{
                        print dic_clust[clust][i]["header"]":"dic_clust[clust][1]["RS"]":IMPRECISE:"i;
                        print dic_clust[clust][i]["seq"];
                    }}
                    clust=$NF;
                    dic_clust[clust][1]["header"]=">"chrom":<INS>:"start":"end":TrEMOLO.INS."clust;
                    dic_clust[clust][1]["seq"]=$5;
                    dic_clust[clust][1]["RS"]=1;
                    dic_clust[clust][1]["READ"]=$4;
                    print chrom":<INS>:"start":"end":TrEMOLO.INS."clust":1:"dic_clust[clust][1]["READ"]":"$7":"$8":"($7+$8) >> FILE_RD_NUM;
                }} else {{
                    present=0;
                    for (i=1; i<dic_clust[clust][1]["RS"]+1; i++) {{
                        if (dic_clust[clust][i]["READ"] == $4) {{present=1; i=dic_clust[clust][1]["RS"]+1}}
                    }}
                    if (present == 0) {{
                        dic_clust[clust][1]["RS"]+=1;
                        n=dic_clust[clust][1]["RS"];
                        dic_clust[clust][n]["READ"]=$4;
                        dic_clust[clust][n]["header"]=">"chrom":<INS>:"start":"end":TrEMOLO.INS."clust;
                        dic_clust[clust][n]["seq"]=$5;
                        print chrom":<INS>:"start":"end":TrEMOLO.INS."clust":"n":"dic_clust[clust][n]["READ"]":"$7":"$8":"($7+$8) >> FILE_RD_NUM;
                    }}
                }}
            }}' {output.clustered:q} > {output.fasta:q}
        awk '/^>/ {{head=substr($0,2,length($0))}} /^[^>]/ && OFS="\t" {{print head,length($0)}}' \
            {output.fasta:q} > {output.sizes:q}
        test -s {output.fasta:q}
        """


rule blast_outsider_alignment_insertions:
    input:
        query=OUTSIDER_INS_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=OUTSIDER_INS_BLAST,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_blast_alignment_insertions.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_blast_alignment_insertions.tsv",
    shell:
        """
        blastn -num_threads {threads} -db {input.database:q} \
            -query {input.query:q} -outfmt 6 -out {output.blast:q} \
            > {log:q} 2>&1
        """


rule classify_outsider_alignment_insertions:
    input:
        blast=OUTSIDER_INS_BLAST,
        database_index=TE_FASTA_INDEX,
    output:
        calls=OUTSIDER_INS_CSV,
        combined=OUTSIDER_INS_COMBINE,
        counts=OUTSIDER_INS_COUNT,
        bed=OUTSIDER_INS_BED,
        ids=OUTSIDER_INS_IDS,
        read_counts=OUTSIDER_INS_READ_COUNTS,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/parse_blast_main.py"),
        options=OUTSIDER_BLAST_FILTER,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_classify_alignment_insertions.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_classify_alignment_insertions.tsv",
    shell:
        """
        set -euo pipefail
        python3 {params.script:q} {input.blast:q} {input.database_index:q} \
            {output.calls:q} --combine_name {output.combined:q} -c \
            {params.options} > {log:q} 2>&1
        awk 'NR>1 && OFS="\t" {{split($2,a,":"); print a[1],a[3],a[3]+1,$1"|"a[5],$5,a[9]}}' \
            {output.calls:q} > {output.bed:q}
        awk 'NR>1 {{split($2,a,":"); if (a[2]=="<INS>") print a[5]}}' \
            {output.combined:q} > {output.ids:q}
        grep -w -f {output.ids:q} {input.blast:q} \
            | awk '{{print $1":"$2}}' | sort -u \
            | awk -F ':' '{{print $5":"$9}}' | sort | uniq -c \
            | awk 'OFS="\t" {{print $2,$1}}' > {output.read_counts:q}
        """


rule extract_outsider_sniffles_sequences:
    input:
        vcf=SV_VCF,
        bam=MAPPING_BAM,
        bai=MAPPING_BAM + ".bai",
        direct_sizes=OUTSIDER_INS_SIZE_BASE,
    output:
        raw=temp(OUTSIDER_SNIFFLES_RAW_FASTA),
        variants=OUTSIDER_SNIFFLES_VARIANTS,
        fasta=OUTSIDER_SNIFFLES_FASTA,
        read_counts=OUTSIDER_SNIFFLES_READ_COUNTS,
        sizes=OUTSIDER_SV_SIZE,
    params:
        extract_script=str(PIPELINE_ROOT / "lib/python/parsing/get_seq_vcf.py"),
        reads_script=str(PIPELINE_ROOT / "lib/python/parsing/parse_bam_found_ins.py"),
        options=OUTSIDER_GET_SEQ_OPTIONS,
        chrom=OUTSIDER_CHROM_KEEP,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_extract_sniffles_sequences.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_extract_sniffles_sequences.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {OUTSIDER_TE_DIR}
        python3 {params.extract_script:q} {params.options} -c {params.chrom:q} \
            {input.vcf:q} {output.raw:q} > {log:q} 2>&1
        awk 'BEGIN {{OFS="\t"}} substr($0,1,1)==">" {{split($0,a,":"); header=substr(a[1],2) OFS a[3] OFS a[4] OFS substr($0,2); next}} {{print header,length($0),$0}}' \
            {output.raw:q} | grep -w -E 'INS|DEL' > {output.variants:q}
        python3 {params.reads_script:q} -r {output.read_counts:q} \
            {input.bam:q} {output.variants:q} > {output.fasta:q} 2>> {log:q}
        cp {input.direct_sizes:q} {output.sizes:q}
        awk '/^>/ {{head=substr($0,2,length($0))}} /^[^>]/ && OFS="\t" {{print head,length($0)}}' \
            {output.fasta:q} >> {output.sizes:q}
        test -s {output.fasta:q}
        """


rule blast_outsider_sniffles_sequences:
    input:
        query=OUTSIDER_SNIFFLES_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=OUTSIDER_SNIFFLES_BLAST,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_blast_sniffles.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_blast_sniffles.tsv",
    shell:
        """
        blastn -num_threads {threads} -db {input.database:q} \
            -query {input.query:q} -outfmt 6 -out {output.blast:q} \
            > {log:q} 2>&1
        """


rule classify_outsider_sniffles_te:
    input:
        blast=OUTSIDER_SNIFFLES_BLAST,
        database_index=TE_FASTA_INDEX,
    output:
        calls=OUTSIDER_SNIFFLES_CSV,
        combined=OUTSIDER_SNIFFLES_COMBINE,
        counts=OUTSIDER_SNIFFLES_COUNT,
        bed=OUTSIDER_SNIFFLES_BED,
        ids=OUTSIDER_SNIFFLES_IDS,
        read_counts=OUTSIDER_SNIFFLES_TE_READ_COUNTS,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/parse_blast_main.py"),
        options=OUTSIDER_BLAST_FILTER,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_classify_sniffles.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_classify_sniffles.tsv",
    shell:
        """
        set -euo pipefail
        python3 {params.script:q} {input.blast:q} {input.database_index:q} \
            {output.calls:q} --combine_name {output.combined:q} -c \
            {params.options} > {log:q} 2>&1
        awk 'NR>1 {{print $2":"$1":"$5}}' {output.calls:q} \
            | awk -F ':' 'BEGIN {{OFS="\t"}} {{print $1,$3,($3>=$4 ? $3+1 : $4),$10"|"$5,$11,$9}}' \
            | bedtools sort > {output.bed:q}
        awk 'NR>1 {{split($2,a,":"); if (a[2]=="<INS>") print a[5]}}' \
            {output.combined:q} > {output.ids:q}
        grep -w -f {output.ids:q} {input.blast:q} \
            | awk '{{print $1":"$2}}' | sort -u \
            | awk -F ':' '{{print $5":"$9}}' | sort | uniq -c \
            | awk 'OFS="\t" {{print $2,$1}}' > {output.read_counts:q}
        """


rule prepare_outsider_soft_clipped_sequences:
    input:
        candidates=OUTSIDER_SOFT_VCF,
        all_candidates=OUTSIDER_SOFT_VCF_ALL,
        database=PREPARED_TE_DATABASE,
    output:
        variants=OUTSIDER_SOFT_SV_BED,
        fasta=OUTSIDER_SOFT_FASTA,
        all_fasta=OUTSIDER_SOFT_ALL_FASTA,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_prepare_soft_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_prepare_soft_clipped.tsv",
    shell:
        """
        set -euo pipefail
        max_size=$(awk 'substr($0,1,1)!=">" && length($0)>max {{max=length($0)}} END {{print max+100}}' {input.database:q})
        : > {output.variants:q}
        : > {output.fasta:q}
        : > {output.all_fasta:q}
        if test -s {input.candidates:q}; then
            awk -v max_size="$max_size" -v out_bed={output.variants:q} '
                {{
                    split($4,left,";"); split($5,right,";");
                    split(left[3],left_seq,"="); split(right[3],right_seq,"=");
                    split($6,left_reads,"="); split($7,right_reads,"=");
                    n_left=split(left_reads[2],left_names,",");
                    n_right=split(right_reads[2],right_names,",");
                    if (left_seq[2] != "NONE") {{
                        seq=left_seq[2];
                        if (length(seq)>max_size) seq=substr(seq,length(seq)-max_size+1,max_size);
                        print ">"$1":<SOFT>:"$2":"$2+1":"$3".L:"n_left":PRECISE:0";
                        print seq;
                    }}
                    if (right_seq[2] != "NONE") {{
                        seq=right_seq[2];
                        if (length(seq)>max_size) seq=substr(seq,1,max_size);
                        print ">"$1":<SOFT>:"$2":"$2+1":"$3".R:"n_right":PRECISE:0";
                        print seq;
                    }}
                    print $1"\t"$2"\t"$2+1"\t"$3":"(n_left+n_right) >> out_bed;
                }}' {input.candidates:q} | grep -v '^--' > {output.fasta:q}
        fi
        if test -s {input.all_candidates:q}; then
            awk 'BEGIN {{id=""; count=1}} NR>1 && OFS=":" {{
                if ($4 != id) {{id=$4; count=1}}
                print ">"$1,$2"-"$3,$4,$5,count,$6,$8,$9;
                print $10;
                count+=1;
            }}' {input.all_candidates:q} > {output.all_fasta:q}
        fi
        """


rule blast_outsider_soft_clipped_sequences:
    input:
        fasta=OUTSIDER_SOFT_FASTA,
        all_fasta=OUTSIDER_SOFT_ALL_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=OUTSIDER_SOFT_BLAST,
        all_blast=OUTSIDER_SOFT_ALL_BLAST,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_blast_soft_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_blast_soft_clipped.tsv",
    shell:
        """
        set -euo pipefail
        : > {output.blast:q}
        : > {output.all_blast:q}
        if test -s {input.fasta:q}; then
            blastn -num_threads {threads} -db {input.database:q} \
                -query {input.fasta:q} -outfmt 6 -out {output.blast:q} \
                > {log:q} 2>&1
        fi
        if test -s {input.all_fasta:q}; then
            blastn -num_threads {threads} -db {input.database:q} \
                -query {input.all_fasta:q} -outfmt 6 -out {output.all_blast:q} \
                >> {log:q} 2>&1
        fi
        """


rule classify_outsider_soft_clipped_te:
    input:
        blast=OUTSIDER_SOFT_BLAST,
        database_index=TE_FASTA_INDEX,
    output:
        calls=OUTSIDER_SOFT_CSV,
        combined=OUTSIDER_SOFT_COMBINE,
        counts=OUTSIDER_SOFT_COUNT,
        bed=OUTSIDER_SOFT_BED,
        match_bed=OUTSIDER_SOFT_MATCH_BED,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/parse_blast_main.py"),
        options=OUTSIDER_BLAST_FILTER,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_classify_soft_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_classify_soft_clipped.tsv",
    shell:
        """
        set -euo pipefail
        : > {output.calls:q}
        : > {output.combined:q}
        : > {output.counts:q}
        : > {output.bed:q}
        : > {output.match_bed:q}
        if test -s {input.blast:q}; then
            python3 {params.script:q} {input.blast:q} {input.database_index:q} \
                {output.calls:q} --combine_name {output.combined:q} -c \
                {params.options} -k SOFT > {log:q} 2>&1
            awk 'NR>1 {{print $2":"$1}}' {output.calls:q} \
                | awk -F ':' 'BEGIN {{OFS="\t"}} {{print $1,$3,$4,$10"|"$5}}' \
                > {output.match_bed:q}
            awk 'NR>1 && OFS="\t" {{split($2,a,":"); print a[1],a[3],a[4],$1"|"a[5],$5,a[9]}}' \
                {output.calls:q} > {output.bed:q}
        fi
        """


rule index_outsider_fastq:
    input:
        reads=PREPARED_SAMPLE,
    output:
        index=OUTSIDER_READS_INDEX,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/outsider_fastq_index.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_fastq_index.tsv",
    shell:
        "samtools fqidx {input.reads:q} > {log:q} 2>&1"


rule extract_outsider_hard_clipped_reads:
    input:
        candidates=OUTSIDER_HARD_VCF,
        reads=PREPARED_SAMPLE,
        reads_index=OUTSIDER_READS_INDEX,
    output:
        ids=OUTSIDER_HARD_IDS,
        fastq=OUTSIDER_HARD_FASTQ,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_extract_hard_reads.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_extract_hard_reads.tsv",
    shell:
        """
        set -euo pipefail
        awk 'NR>1 {{print $5}}' {input.candidates:q} | sort -u > {output.ids:q}
        : > {output.fastq:q}
        if test -s {output.ids:q}; then
            samtools fqidx {input.reads:q} -r {output.ids:q} \
                > {output.fastq:q} 2> {log:q}
        fi
        """


rule prepare_outsider_hard_clipped_sequences:
    input:
        candidates=OUTSIDER_HARD_VCF,
        fastq=OUTSIDER_HARD_FASTQ,
        database=PREPARED_TE_DATABASE,
    output:
        reads_fasta=OUTSIDER_HARD_READS_FASTA,
        reads_index=OUTSIDER_HARD_READS_FASTA_INDEX,
        variants=OUTSIDER_HARD_RAW_BED,
        fasta=OUTSIDER_HARD_FASTA,
    params:
        convert_script=str(PIPELINE_ROOT / "lib/python/format_files/fastq_to_fasta.py"),
        extract_script=str(PIPELINE_ROOT / "lib/python/parsing/get_seq_hard.py"),
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_prepare_hard_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_prepare_hard_clipped.tsv",
    shell:
        """
        set -euo pipefail
        max_size=$(awk 'substr($0,1,1)!=">" && length($0)>max {{max=length($0)}} END {{print max+100}}' {input.database:q})
        : > {output.reads_fasta:q}
        : > {output.reads_index:q}
        : > {output.variants:q}
        : > {output.fasta:q}
        if test -s {input.fastq:q}; then
            python3 {params.convert_script:q} {input.fastq:q} {output.reads_fasta:q} \
                > {log:q} 2>&1
            grep -n '>' {output.reads_fasta:q} | tr -d '>' > {output.reads_index:q}
        fi
        if test -s {output.reads_index:q} && test -s {input.candidates:q}; then
            python3 {params.extract_script:q} -s "$max_size" {output.reads_index:q} \
                {output.reads_fasta:q} {input.candidates:q} > {output.variants:q} \
                2>> {log:q}
        fi
        if test -s {output.variants:q}; then
            bedtools sort -i {output.variants:q} | bedtools cluster -d 100 \
                | awk 'BEGIN {{count=0; id=""}} OFS=":" {{
                    if (id != $8) {{
                        if (NR>1) for (i=1; i<count+1; i++) {{print tab_head[1],count,"IMPRECISE",(i-1); print tab_seq[i]}}
                        id=$8; count=0; delete tab_head; delete tab_seq;
                    }} else {{
                        tab_head[count+1]=">"$1":<HARD>:"$2":"$3":HARD."id"."$6;
                        tab_seq[count+1]=$5; count+=1;
                    }}
                }}' > {output.fasta:q} 2>> {log:q}
        fi
        """


rule blast_outsider_hard_clipped_sequences:
    input:
        fasta=OUTSIDER_HARD_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        blast=OUTSIDER_HARD_BLAST,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_blast_hard_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_blast_hard_clipped.tsv",
    shell:
        """
        set -euo pipefail
        : > {output.blast:q}
        if test -s {input.fasta:q}; then
            blastn -num_threads {threads} -db {input.database:q} \
                -query {input.fasta:q} -outfmt 6 -out {output.blast:q} \
                > {log:q} 2>&1
        fi
        """


rule classify_outsider_hard_clipped_te:
    input:
        blast=OUTSIDER_HARD_BLAST,
        database_index=TE_FASTA_INDEX,
    output:
        calls=OUTSIDER_HARD_CSV,
        combined=OUTSIDER_HARD_COMBINE,
        counts=OUTSIDER_HARD_COUNT,
        bed=OUTSIDER_HARD_BED,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/parse_blast_main.py"),
        options=OUTSIDER_BLAST_FILTER,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/outsider_classify_hard_clipped.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_classify_hard_clipped.tsv",
    shell:
        """
        set -euo pipefail
        : > {output.calls:q}
        : > {output.combined:q}
        : > {output.counts:q}
        : > {output.bed:q}
        if test -s {input.blast:q}; then
            python3 {params.script:q} {input.blast:q} {input.database_index:q} \
                {output.calls:q} --combine_name {output.combined:q} -c \
                {params.options} -k HARD > {log:q} 2>&1
            awk 'NR>1 {{print $2":"$1":"$5}}' {output.calls:q} \
                | awk -F ':' 'BEGIN {{OFS="\t"}} {{print $1,$3,$4,$10"|"$5,$11,$9}}' \
                > {output.bed:q}
        fi
        """


rule merge_outsider_te_evidence:
    input:
        sniffles=OUTSIDER_SNIFFLES_BED,
        direct=OUTSIDER_INS_BED,
        soft=OUTSIDER_SOFT_BED,
        hard=OUTSIDER_HARD_BED,
    output:
        sniffles_cluster=f"{OUTSIDER_MERGE_DIR}/tmp_POSITION_TE_OUTSIDER_CLUSTER_100.bed",
        direct_cluster=f"{OUTSIDER_MERGE_DIR}/tmp_TrEMOLO_TE_OUTSIDER_CLUSTER_100.bed",
        direct_ids=f"{OUTSIDER_MERGE_DIR}/tmp_ID_TrEMOLO.txt",
        direct_only=f"{OUTSIDER_MERGE_DIR}/tmp_TE_TrEMOLO_NOT_FOUND_IN_sniffles.bed",
        soft_cluster=f"{OUTSIDER_MERGE_DIR}/tmp_SOFT_CLUSTER_100.bed",
        soft_ids=f"{OUTSIDER_MERGE_DIR}/tmp_ID_SOFT_commun.txt",
        soft_only=f"{OUTSIDER_MERGE_DIR}/tmp_TE_SOFT_NOT_FOUND_IN_sniffles_and_TrEMOLO.bed",
        hard_cluster=f"{OUTSIDER_MERGE_DIR}/tmp_HARD_CLUSTER_100.bed",
        hard_ids=f"{OUTSIDER_MERGE_DIR}/tmp_ID_HARD_commun.txt",
        hard_only=f"{OUTSIDER_MERGE_DIR}/tmp_TE_HARD_NOT_FOUND_IN_SOFT_sniffles_and_TrEMOLO.bed",
        merged=OUTSIDER_MERGED_BED,
        merged_copy=OUTSIDER_MERGED_BED_COPY,
        counts=OUTSIDER_MERGED_COUNT,
        positions=OUTSIDER_POSITION_BED,
    params:
        chrom=OUTSIDER_CHROM_KEEP,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/outsider_merge_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/outsider_merge_te.tsv",
    shell:
        """
        set -eu
        mkdir -p {OUTSIDER_MERGE_DIR}
        chrom_keep=$(printf '%s' {params.chrom:q} | tr ',' '|')
        : > {output.sniffles_cluster:q}
        : > {output.direct_cluster:q}
        : > {output.direct_ids:q}
        : > {output.direct_only:q}
        : > {output.soft_cluster:q}
        : > {output.soft_ids:q}
        : > {output.soft_only:q}
        : > {output.hard_cluster:q}
        : > {output.hard_ids:q}
        : > {output.hard_only:q}

        if test -s {input.sniffles:q}; then
            bedtools sort -i {input.sniffles:q} | grep -E "$chrom_keep" \
                | bedtools cluster -d 100 \
                | awk 'BEGIN {{te=""; cluster=""}} {{split($4,a,"|"); if (cluster==$7 && a[1]==te) {{te=a[1]; cluster=$7}} else {{print; te=a[1]; cluster=$7}}}}' \
                | cut -f 1-6 > {output.sniffles_cluster:q}
        fi
        if test -s {input.direct:q}; then
            bedtools sort -i {input.direct:q} | grep -E "$chrom_keep" \
                | bedtools cluster -d 100 \
                | awk 'BEGIN {{te=""; cluster=""}} {{split($4,a,"|"); if (cluster==$7 && a[1]==te) {{te=a[1]; cluster=$7}} else {{print; te=a[1]; cluster=$7}}}}' \
                | cut -f 1-6 > {output.direct_cluster:q}
        fi
        if test -s {output.direct_cluster:q} && test -s {output.sniffles_cluster:q}; then
            bedtools window -a {output.direct_cluster:q} -b {output.sniffles_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u > {output.direct_ids:q}
        fi
        if test -s {output.direct_cluster:q}; then
            grep -w -v -f {output.direct_ids:q} {output.direct_cluster:q} \
                > {output.direct_only:q} || true
        fi

        if test -s {input.soft:q}; then
            bedtools sort -i {input.soft:q} | grep -E "$chrom_keep" \
                | bedtools cluster -d 100 \
                | awk 'BEGIN {{te=""; cluster=""}} {{split($4,a,"|"); if (cluster==$7 && a[1]==te) {{te=a[1]; cluster=$7}} else {{print; te=a[1]; cluster=$7}}}}' \
                | cut -f 1-6 > {output.soft_cluster:q}
        fi
        if test -s {output.soft_cluster:q} && test -s {output.sniffles_cluster:q}; then
            bedtools window -a {output.soft_cluster:q} -b {output.sniffles_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u > {output.soft_ids:q}
        fi
        if test -s {output.soft_cluster:q} && test -s {output.direct_cluster:q}; then
            bedtools window -a {output.soft_cluster:q} -b {output.direct_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u >> {output.soft_ids:q}
        fi
        if test -s {output.soft_cluster:q}; then
            grep -w -v -f {output.soft_ids:q} {output.soft_cluster:q} \
                > {output.soft_only:q} || true
        fi

        if test -s {input.hard:q}; then
            bedtools sort -i {input.hard:q} | grep -E "$chrom_keep" \
                | bedtools cluster -d 100 \
                | awk 'BEGIN {{te=""; cluster=""}} {{split($4,a,"|"); if (cluster==$7 && a[1]==te) {{te=a[1]; cluster=$7}} else {{print; te=a[1]; cluster=$7}}}}' \
                | cut -f 1-6 > {output.hard_cluster:q}
        fi
        if test -s {output.hard_cluster:q} && test -s {output.sniffles_cluster:q}; then
            bedtools window -a {output.hard_cluster:q} -b {output.sniffles_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u > {output.hard_ids:q}
        fi
        if test -s {output.hard_cluster:q} && test -s {output.direct_cluster:q}; then
            bedtools window -a {output.hard_cluster:q} -b {output.direct_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u >> {output.hard_ids:q}
        fi
        if test -s {output.hard_cluster:q} && test -s {output.soft_cluster:q}; then
            bedtools window -a {output.hard_cluster:q} -b {output.soft_cluster:q} -w 100 \
                | awk '{{split($4,a,"|"); split($10,b,"|"); if (a[1]==b[1]) print $4}}' \
                | sort -u >> {output.hard_ids:q}
        fi
        if test -s {output.hard_cluster:q}; then
            grep -w -v -f {output.hard_ids:q} {output.hard_cluster:q} \
                | grep -E "$chrom_keep" > {output.hard_only:q} || true
        fi

        : > {output.merged:q}
        for evidence in {output.sniffles_cluster:q} {output.direct_only:q} \
            {output.soft_only:q} {output.hard_only:q}; do
            if test -s "$evidence"; then
                grep -E "$chrom_keep" "$evidence" >> {output.merged:q} || true
            fi
        done
        cp {output.merged:q} {output.merged_copy:q}
        cut -f 4 {output.merged:q} | cut -d '|' -f 1 | sort | uniq -c \
            | awk 'BEGIN {{OFS="\t"}} {{print $1,$2}}' | sort -k 1 -n \
            | awk 'BEGIN {{print "x\ty\tz"}} BEGIN {{OFS="\t"}} {{print "",$2,$1}}' \
            > {output.counts:q}
        rm -f {output.positions:q}
        ln -srf {output.merged:q} {output.positions:q}
        printf 'sniffles=%s direct=%s soft=%s hard=%s merged=%s\n' \
            "$(wc -l < {output.sniffles_cluster:q})" \
            "$(wc -l < {output.direct_only:q})" \
            "$(wc -l < {output.soft_only:q})" \
            "$(wc -l < {output.hard_only:q})" \
            "$(wc -l < {output.merged:q})" > {log:q}
        """
