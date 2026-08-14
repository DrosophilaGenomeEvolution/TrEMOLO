PREPARED_TE_DATABASE = f"{WORKDIR}/INPUT/te_database.fasta"
TE_FASTA_INDEX = PREPARED_TE_DATABASE + ".fai"
TE_BLAST_NHR = PREPARED_TE_DATABASE + ".nhr"
TE_BLAST_NIN = PREPARED_TE_DATABASE + ".nin"
TE_BLAST_NSQ = PREPARED_TE_DATABASE + ".nsq"
TE_SIZE_TABLE = f"{WORKDIR}/1-UTILS/TE_SIZE.tsv"


rule te_database_indexes:
    input:
        TE_FASTA_INDEX,
        TE_BLAST_NHR,
        TE_BLAST_NIN,
        TE_BLAST_NSQ,
        TE_SIZE_TABLE,


rule measure_te_sequences:
    input:
        fasta=PREPARED_TE_DATABASE,
    output:
        sizes=TE_SIZE_TABLE,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/te_sequence_sizes.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_sequence_sizes.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {WORKDIR}/1-UTILS {WORKDIR}/log {WORKDIR}/benchmarks
        awk '/^>/ {{header=substr($1,2)}} /^[^>]/ && OFS="\t" {{print header,length($0)}}' \
            {input.fasta:q} > {output.sizes:q} 2> {log:q}
        """


rule index_te_fasta:
    input:
        fasta=PREPARED_TE_DATABASE,
    output:
        fai=TE_FASTA_INDEX,
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/te_fasta_index.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_fasta_index.tsv",
    shell:
        "samtools faidx {input.fasta:q} > {log:q} 2>&1"


rule make_te_blast_database:
    input:
        fasta=PREPARED_TE_DATABASE,
    output:
        nhr=TE_BLAST_NHR,
        nin=TE_BLAST_NIN,
        nsq=TE_BLAST_NSQ,
    threads: 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/te_blast_database.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_blast_database.tsv",
    shell:
        "makeblastdb -in {input.fasta:q} -dbtype nucl > {log:q} 2>&1"
