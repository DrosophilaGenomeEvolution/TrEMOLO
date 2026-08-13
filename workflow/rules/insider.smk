INSIDER_VARIANT_DIR = f"{WORKDIR}/INSIDER/VARIANT_CALLING"
INSIDER_PREFIX = f"{INSIDER_VARIANT_DIR}/assemblytics_out"
INSIDER_SAM = f"{INSIDER_VARIANT_DIR}/pm_against_ref.sam"
INSIDER_DELTA = INSIDER_SAM + ".delta"
INSIDER_UNIQUE_DELTA = INSIDER_PREFIX + ".Assemblytics.unique_length_filtered_l20000.delta"
INSIDER_COORDS_TAB = INSIDER_PREFIX + ".coords.tab"
INSIDER_COORDS_CSV = INSIDER_PREFIX + ".coords.csv"
INSIDER_STATS = INSIDER_PREFIX + ".Assemblytics_assembly_stats.txt"
INSIDER_BETWEEN = INSIDER_PREFIX + ".variants_between_alignments.bed"
INSIDER_WITHIN = INSIDER_PREFIX + ".variants_within_alignments.bed"
INSIDER_BED = INSIDER_PREFIX + ".Assemblytics_structural_variants.bed"
INSIDER_TE_DIR = f"{WORKDIR}/INSIDER/TE_DETECTION"
INSIDER_INSERTION_FASTA = f"{INSIDER_TE_DIR}/INSERTION_SEQ.fasta"
INSIDER_DELETION_FASTA = f"{INSIDER_TE_DIR}/DELETION_SEQ.fasta"
INSIDER_INSERTION_BLAST = f"{INSIDER_TE_DIR}/INSERTION.bln"
INSIDER_DELETION_BLAST = f"{INSIDER_TE_DIR}/DELETION.bln"
INSIDER_INSERTION_CSV = f"{INSIDER_TE_DIR}/INSERTION.csv"
INSIDER_DELETION_CSV = f"{INSIDER_TE_DIR}/DELETION.csv"
INSIDER_INSERTION_COMBINE = f"{INSIDER_TE_DIR}/INSERTION_COMBINE_TE.csv"
INSIDER_DELETION_COMBINE = f"{INSIDER_TE_DIR}/DELETION_COMBINE_TE.csv"
INSIDER_INSERTION_TE_BED = f"{INSIDER_TE_DIR}/INSERTION_TE.bed"
INSIDER_DELETION_TE_BED = f"{INSIDER_TE_DIR}/DELETION_TE.bed"

INSIDER_MINIMAP2_PRESET = (
    INSIDER_PARAMS.get("MINIMAP2", {}).get("PRESET_OPTION") or "asm5"
)
INSIDER_MINIMAP2_OPTIONS = option_without_threads(
    INSIDER_PARAMS.get("MINIMAP2", {}).get("OPTION", "")
)
INSIDER_BLAST_FILTER = INSIDER_PARAMS.get(
    "PARS_BLN_OPTION", "--min-size-percent 80 --min-pident 80"
)


rule insider_variant_calling:
    input:
        INSIDER_BED,
        INSIDER_STATS,


rule insider_te_detection:
    input:
        INSIDER_INSERTION_CSV,
        INSIDER_DELETION_CSV,
        INSIDER_INSERTION_TE_BED,
        INSIDER_DELETION_TE_BED,


rule map_assembly_insider:
    input:
        reference=PREPARED_REFERENCE,
        genome=PREPARED_GENOME,
    output:
        sam=INSIDER_SAM,
    params:
        preset=INSIDER_MINIMAP2_PRESET,
        options=INSIDER_MINIMAP2_OPTIONS,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_minimap2.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_minimap2.tsv",
    shell:
        """
        mkdir -p {INSIDER_VARIANT_DIR} {WORKDIR}/log {WORKDIR}/benchmarks
        minimap2 -ax {params.preset} -t {threads} {params.options} \
            {input.reference:q} {input.genome:q} > {output.sam:q} 2> {log:q}
        test -s {output.sam:q}
        """


rule sam_to_delta_insider:
    input:
        sam=INSIDER_SAM,
    output:
        delta=INSIDER_DELTA,
    params:
        script=str(PIPELINE_ROOT / "lib/python/assemblitics/sam2delta.py"),
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_sam2delta.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_sam2delta.tsv",
    shell:
        "python3 {params.script:q} {input.sam:q} > {log:q} 2>&1"


rule assemblytics_unique_anchor_insider:
    input:
        delta=INSIDER_DELTA,
    output:
        unique_delta=INSIDER_UNIQUE_DELTA,
        coords_tab=INSIDER_COORDS_TAB,
        coords_csv=INSIDER_COORDS_CSV,
        stats=INSIDER_STATS,
    params:
        script=str(PIPELINE_ROOT / "lib/python/assemblitics/Assemblytics_uniq_anchor.py"),
        prefix=INSIDER_PREFIX,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_unique_anchor.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_unique_anchor.tsv",
    shell:
        """
        python3 {params.script:q} --delta {input.delta:q} --unique-length 20000 \
            --out {params.prefix:q} --keep-small-uniques > {log:q} 2>&1
        """


rule assemblytics_between_alignments_insider:
    input:
        coords=INSIDER_COORDS_TAB,
    output:
        bed=INSIDER_BETWEEN,
    params:
        script=str(PIPELINE_ROOT / "lib/perl/Assemblytics_between_alignments.pl"),
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_between_alignments.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_between_alignments.tsv",
    shell:
        """
        perl {params.script:q} {input.coords:q} 50 20000 all-chromosomes \
            exclude-longrange bed > {output.bed:q} 2> {log:q}
        """


rule assemblytics_within_alignment_insider:
    input:
        delta=INSIDER_UNIQUE_DELTA,
    output:
        bed=INSIDER_WITHIN,
    params:
        script=str(PIPELINE_ROOT / "lib/python/assemblitics/Assemblytics_within_alignment.py"),
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_within_alignment.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_within_alignment.tsv",
    shell:
        "python3 {params.script:q} --delta {input.delta:q} --min 50 > {output.bed:q} 2> {log:q}"


rule merge_and_annotate_insider_variants:
    input:
        between=INSIDER_BETWEEN,
        within=INSIDER_WITHIN,
        reference=PREPARED_REFERENCE,
        genome=PREPARED_GENOME,
        stats=INSIDER_STATS,
    output:
        bed=INSIDER_BED,
        stats=touch(f"{INSIDER_VARIANT_DIR}/.stats_renamed"),
    params:
        script=str(PIPELINE_ROOT / "lib/python/assemblitics/filter_gap_SVs.py"),
        variant_dir=INSIDER_VARIANT_DIR,
        reference_name=Path(REFERENCE).name,
        genome_name=Path(GENOME).name,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_filter_gaps.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_filter_gaps.tsv",
    shell:
        """
        set -euo pipefail
        cat {input.between:q} {input.within:q} > {output.bed:q}
        test -s {output.bed:q}
        old_directory=$PWD
        log_path=$(readlink -m {log:q})
        reference=$(readlink -f {input.reference:q})
        genome=$(readlink -f {input.genome:q})
        cd {params.variant_dir:q}
        python3 {params.script:q} "$reference" "$genome" > "$log_path" 2>&1
        cd "$old_directory"
        sed -i "s/^Reference: file1$/Reference: {params.reference_name}/" {input.stats:q}
        sed -i "s/^Query: file2$/Query: {params.genome_name}/" {input.stats:q}
        """


rule extract_insider_sv_sequences:
    input:
        variants=INSIDER_BED,
        reference=PREPARED_REFERENCE,
        genome=PREPARED_GENOME,
    output:
        insertion=INSIDER_INSERTION_FASTA,
        deletion=INSIDER_DELETION_FASTA,
    params:
        directory=INSIDER_TE_DIR,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_extract_sv_sequences.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_extract_sv_sequences.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {params.directory:q}
        : > {output.insertion:q}
        : > {output.deletion:q}
        for type in INSERTION Repeat_expansion Tandem_expansion; do
            bed={params.directory:q}/$type.bed
            fasta={params.directory:q}/$type.seq.fasta
            grep -i "$type" {input.variants:q} \
                | awk '{{print $10":"$4}}' \
                | awk -F ':' -v type="$type" 'BEGIN {{OFS="\t"}} {{split($2,a,"-"); print $1,a[1],a[2],$4":"$3":"type}}' \
                > "$bed" || true
            if test -s "$bed"; then
                bedtools getfasta -fi {input.genome:q} -bed "$bed" -name+ \
                    > "$fasta" 2>> {log:q}
                cat "$fasta" >> {output.insertion:q}
            fi
        done
        for type in DELETION Repeat_contraction Tandem_contraction; do
            bed={params.directory:q}/$type.bed
            fasta={params.directory:q}/$type.seq.fasta
            grep -i "$type" {input.variants:q} \
                | awk -v type="$type" 'BEGIN {{OFS="\t"}} {{print $1,$2,$3,$4":"$6":"type}}' \
                > "$bed" || true
            if test -s "$bed"; then
                bedtools getfasta -fi {input.reference:q} -bed "$bed" -name+ \
                    > "$fasta" 2>> {log:q}
                cat "$fasta" >> {output.deletion:q}
            fi
        done
        test -s {output.insertion:q}
        test -s {output.deletion:q}
        """


rule blast_insider_sv_against_te:
    input:
        insertion=INSIDER_INSERTION_FASTA,
        deletion=INSIDER_DELETION_FASTA,
        database=PREPARED_TE_DATABASE,
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    output:
        insertion=INSIDER_INSERTION_BLAST,
        deletion=INSIDER_DELETION_BLAST,
    threads: THREADS
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/insider_blast_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_blast_te.tsv",
    shell:
        """
        set -euo pipefail
        blastn -num_threads {threads} -db {input.database:q} \
            -query {input.insertion:q} -outfmt 6 -out {output.insertion:q} \
            2> {log:q}
        blastn -num_threads {threads} -db {input.database:q} \
            -query {input.deletion:q} -outfmt 6 -out {output.deletion:q} \
            2>> {log:q}
        """


rule classify_insider_te:
    input:
        insertion=INSIDER_INSERTION_BLAST,
        deletion=INSIDER_DELETION_BLAST,
        database=PREPARED_TE_DATABASE,
    output:
        insertion=INSIDER_INSERTION_CSV,
        deletion=INSIDER_DELETION_CSV,
        insertion_combine=INSIDER_INSERTION_COMBINE,
        deletion_combine=INSIDER_DELETION_COMBINE,
    params:
        script=str(PIPELINE_ROOT / "lib/python/parsing/global_sv.py"),
        options=INSIDER_BLAST_FILTER,
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/insider_classify_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_classify_te.tsv",
    shell:
        """
        set -euo pipefail
        python3 {params.script:q} {params.options} \
            --combine_name {output.insertion_combine:q} \
            {input.insertion:q} {input.database:q} {output.insertion:q} \
            > {log:q} 2>&1
        python3 {params.script:q} {params.options} \
            --combine_name {output.deletion_combine:q} \
            {input.deletion:q} {input.database:q} {output.deletion:q} \
            >> {log:q} 2>&1
        """


rule format_insider_te_calls:
    input:
        variants=INSIDER_BED,
        insertion=INSIDER_INSERTION_CSV,
        deletion=INSIDER_DELETION_CSV,
        insertion_combine=INSIDER_INSERTION_COMBINE,
    output:
        insertion_bed=INSIDER_INSERTION_TE_BED,
        deletion_bed=INSIDER_DELETION_TE_BED,
        insertion_ref=f"{INSIDER_TE_DIR}/INSERTION_TE_ON_REF.bed",
        deletion_ref=f"{INSIDER_TE_DIR}/DELETION_TE_ON_REF.bed",
        positions=f"{WORKDIR}/POSITION_TE_INSIDER.bed",
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/insider_format_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/insider_format_te.tsv",
    shell:
        """
        set -euo pipefail
        awk 'NR>1 {{split($2,a,":"); split(a[6],b,"-"); print a[5],b[1]+$6,b[2]-((b[2]-b[1])-$7),$1"|"a[1],a[7]}}' \
            OFS='\t' {input.insertion_combine:q} > {output.insertion_bed:q}
        : > {output.insertion_ref:q}
        awk 'NR>1 {{split($2,a,":"); print a[1]}}' {input.insertion:q} | while read -r id; do
            awk -v id="$id" '$4==id {{print $1,$2,$3; exit}}' OFS='\t' {input.variants:q} \
                | while read -r chrom start end; do
                    te=$(awk -F '\t' -v id="$id" 'NR>1 && $2 ~ ("^"id":") {{print $1; exit}}' {input.insertion:q})
                    printf '%s\t%s\t%s\t%s|%s\n' "$chrom" "$start" "$end" "$te" "$id"
                done
        done > {output.insertion_ref:q}
        awk 'NR>1 {{print $1":"$2}}' {input.deletion:q} \
            | awk -F ':' 'BEGIN {{OFS="\t"}} {{split($7,p,"-"); print $6,p[1],p[2],$1"|"$2}}' \
            > {output.deletion_ref:q}
        : > {output.deletion_bed:q}
        awk 'NR>1 {{split($2,a,":"); print a[1]}}' {input.deletion:q} | while read -r id; do
            location=$(awk -F '\t' -v id="$id" '$4==id {{print $10; exit}}' {input.variants:q})
            te=$(awk -F '\t' -v id="$id" 'NR>1 && $2 ~ ("^"id":") {{print $1; exit}}' {input.deletion:q})
            printf '%s\n' "$location" | awk -F ':' -v te="$te" -v id="$id" \
                'BEGIN {{OFS="\t"}} {{split($2,p,"-"); print $1,p[1],p[2],te"|"id,$3}}'
        done > {output.deletion_bed:q}
        bedtools sort -i {output.insertion_bed:q} \
            | awk 'BEGIN {{OFS="\t"}} {{print $1,$2,$3,$4,$3-$2,$5}}' \
            > {output.positions:q}
        """
