"""Optional catalogue of TE evidence in reads with no genomic alignment."""
import json
import shlex

def validate_unmapped_mapping(wildcards):
    if "--sam-hit-only" in shlex.split(MINIMAP2_OPTIONS):
        raise ValueError("UNMAPPED_TE requires unmapped SAM records; remove --sam-hit-only from MINIMAP2.OPTION")
    return ""

UNMAPPED_TE_DIR = f"{WORKDIR}/UNMAPPED_TE"
UNMAPPED_TE_SCRIPT = PIPELINE_ROOT / "lib/python/workflow/unmapped_te.py"
UNMAPPED_TE_FASTA = f"{UNMAPPED_TE_DIR}/unmapped_reads.fasta"
UNMAPPED_TE_READS = f"{UNMAPPED_TE_DIR}/unmapped_reads.tsv"
UNMAPPED_TE_EXTRACTION = f"{UNMAPPED_TE_DIR}/extraction.json"
UNMAPPED_TE_BLAST = f"{UNMAPPED_TE_DIR}/TE_vs_reads.bln"
UNMAPPED_TE_TABLE = f"{UNMAPPED_TE_DIR}/UNMAPPED_TE.tsv"
UNMAPPED_TE_HSPS = f"{UNMAPPED_TE_DIR}/UNMAPPED_TE_HSP.tsv"
UNMAPPED_TE_SUMMARY = f"{UNMAPPED_TE_DIR}/summary.json"
UNMAPPED_TE_README = f"{UNMAPPED_TE_DIR}/README.md"
UNMAPPED_TE_FINAL_OUTPUTS = [UNMAPPED_TE_FASTA, UNMAPPED_TE_READS,
                           UNMAPPED_TE_TABLE, UNMAPPED_TE_HSPS,
                           UNMAPPED_TE_SUMMARY, UNMAPPED_TE_README]


rule unmapped_te:
    input:
        UNMAPPED_TE_FINAL_OUTPUTS,


rule extract_unmapped_reads:
    input:
        sam=MAPPING_SAM,
        fastq=PREPARED_SAMPLE,
        script=str(UNMAPPED_TE_SCRIPT),
        helper=str(PIPELINE_ROOT / "lib/python/workflow/annotate_resident_te.py"),
    output:
        fasta=UNMAPPED_TE_FASTA,
        reads=UNMAPPED_TE_READS,
        stats=UNMAPPED_TE_EXTRACTION,
    params:
        mapping_validation=validate_unmapped_mapping,
        method=UNMAPPED_EXTRACTION_PARAMS["METHOD"],
        name_memory_mb=UNMAPPED_EXTRACTION_PARAMS["NAME_MEMORY_MB"],
        samtools=UNMAPPED_EXTRACTION_PARAMS["SAMTOOLS"],
    threads: THREADS if UNMAPPED_EXTRACTION_PARAMS["METHOD"] == "collate" else 1
    resources:
        mem_mb=1024,
    log:
        f"{WORKDIR}/log/extract_unmapped_reads.log",
    benchmark:
        f"{WORKDIR}/benchmarks/extract_unmapped_reads.tsv",
    shell:
        """
        source {WORKFLOW_LOG_HELPER:q}
        trm_begin extract_unmapped_reads "UNMAPPED: selecting unmapped reads ({params.method}) and recovering FASTQ sequences" {WORKDIR:q} {log:q}
        python3 {input.script:q} extract --sam {input.sam:q} --fastq {input.fastq:q} \
            --fasta {output.fasta:q} --reads {output.reads:q} --stats {output.stats:q} \
            --selection-method {params.method:q} --name-memory-mb {params.name_memory_mb} \
            --samtools {params.samtools:q} --threads {threads} > {log:q} 2>&1
        """


rule blast_te_against_unmapped_reads:
    input:
        fasta=UNMAPPED_TE_FASTA,
        stats=UNMAPPED_TE_EXTRACTION,
        te=PREPARED_TE_DATABASE,
    output:
        blast=UNMAPPED_TE_BLAST,
    params:
        database=f"{UNMAPPED_TE_DIR}/BLAST_DB/reads",
        max_evalue=UNMAPPED_TE_PARAMS["MAX_EVALUE"],
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/blast_te_against_unmapped_reads.log",
    benchmark:
        f"{WORKDIR}/benchmarks/blast_te_against_unmapped_reads.tsv",
    shell:
        """
        source {WORKFLOW_LOG_HELPER:q}
        trm_begin blast_te_against_unmapped_reads "UNMAPPED: matching TE sequences to unmapped reads" {WORKDIR:q} {log:q}
        set -euo pipefail
        if test -s {input.fasta:q}; then
            mkdir -p {UNMAPPED_TE_DIR}/BLAST_DB
            makeblastdb -in {input.fasta:q} -dbtype nucl -out {params.database:q} > {log:q} 2>&1
            count=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["unmapped_reads"])' {input.stats:q})
            blastn -task blastn -num_threads {threads} -evalue {params.max_evalue} \
                -query {input.te:q} -db {params.database:q} -max_target_seqs "$count" \
                -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
                -out {output.blast:q} >> {log:q} 2>&1
        else
            : > {output.blast:q}
            printf 'No unmapped reads; BLAST skipped.\n' > {log:q}
        fi
        """


rule annotate_unmapped_te:
    input:
        blast=UNMAPPED_TE_BLAST,
        reads=UNMAPPED_TE_READS,
        extraction=UNMAPPED_TE_EXTRACTION,
        headers=f"{WORKDIR}/INPUT/te_headers.tsv",
        script=str(UNMAPPED_TE_SCRIPT),
        helper=str(PIPELINE_ROOT / "lib/python/workflow/annotate_resident_te.py"),
        documentation=str(PIPELINE_ROOT / "docs/unmapped-te.md"),
    output:
        table=UNMAPPED_TE_TABLE,
        hsps=UNMAPPED_TE_HSPS,
        stats=UNMAPPED_TE_SUMMARY,
        readme=UNMAPPED_TE_README,
    params:
        parameters=json.dumps(UNMAPPED_TE_PARAMS),
    threads: 1
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/annotate_unmapped_te.log",
    benchmark:
        f"{WORKDIR}/benchmarks/annotate_unmapped_te.tsv",
    shell:
        """
        source {WORKFLOW_LOG_HELPER:q}
        trm_begin annotate_unmapped_te "UNMAPPED: building the TE catalogue" {WORKDIR:q} {log:q}
        set -euo pipefail
        python3 {input.script:q} annotate --blast {input.blast:q} --reads {input.reads:q} \
            --extraction-stats {input.extraction:q} --headers {input.headers:q} \
            --output {output.table:q} --hsps {output.hsps:q} --stats {output.stats:q} \
            --parameters {params.parameters:q} > {log:q} 2>&1
        cp {input.documentation:q} {output.readme:q}
        """
