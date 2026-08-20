TE_GENOME_ROOT = f"{WORKDIR}/TE_GENOME"
TE_GENOME_BLAST_NHR = f"{TE_GENOME_ROOT}/{{target}}/BLAST_DB/target.nhr"
TE_GENOME_BLAST_NIN = f"{TE_GENOME_ROOT}/{{target}}/BLAST_DB/target.nin"
TE_GENOME_BLAST_NSQ = f"{TE_GENOME_ROOT}/{{target}}/BLAST_DB/target.nsq"
TE_GENOME_BLAST = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE_BLAST.tsv"
TE_GENOME_FRAGMENTS = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE_FRAGMENTS.tsv"
TE_GENOME_MATCHES = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE_MATCHES.tsv"
TE_GENOME_COPIES = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE_COPIES.tsv"
TE_GENOME_RELATIONS = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE_RELATIONS.tsv"
TE_GENOME_BED = f"{TE_GENOME_ROOT}/{{target}}/POSITION_ALL_TE.bed"
TE_GENOME_GFF = f"{TE_GENOME_ROOT}/{{target}}/ALL_TE.gff3"
TE_GENOME_PUBLIC_BED = f"{WORKDIR}/POSITION_ALL_TE.bed"
TE_GENOME_ANNOTATOR = PIPELINE_ROOT / "lib/python/workflow/annotate_resident_te.py"

TE_GENOME_TARGET_FASTA = {"GENOME": PREPARED_GENOME}
if REFERENCE:
    TE_GENOME_TARGET_FASTA["REFERENCE"] = PREPARED_REFERENCE


def te_genome_target_fasta(wildcards):
    if wildcards.target not in TE_GENOME_TARGET_FASTA:
        raise ValueError(f"unsupported TE_GENOME target: {wildcards.target}")
    return TE_GENOME_TARGET_FASTA[wildcards.target]


def te_genome_parameter(name, default, value_type, minimum=None, maximum=None):
    value = TE_GENOME_PARAMS.get(name, default)
    if isinstance(value, bool):
        raise ValueError(f"PARAMS.TE_GENOME.{name} must be numeric")
    try:
        converted = value_type(value)
    except (TypeError, ValueError):
        raise ValueError(f"PARAMS.TE_GENOME.{name} must be numeric")
    if not math.isfinite(converted):
        raise ValueError(f"PARAMS.TE_GENOME.{name} must be finite")
    if minimum is not None and converted < minimum:
        raise ValueError(f"PARAMS.TE_GENOME.{name} must be >= {minimum}")
    if maximum is not None and converted > maximum:
        raise ValueError(f"PARAMS.TE_GENOME.{name} must be <= {maximum}")
    return converted


TE_GENOME_CHROM_REGEX = str(TE_GENOME_PARAMS.get("CHROM_KEEP", ".")).replace(",", "|")
try:
    re.compile(TE_GENOME_CHROM_REGEX)
except re.error as error:
    raise ValueError(f"PARAMS.TE_GENOME.CHROM_KEEP is not a valid regular expression: {error}")

TE_GENOME_MIN_PIDENT = te_genome_parameter("MIN_PIDENT", 65.0, float, 0, 100)
TE_GENOME_MIN_ALIGNED_BP = te_genome_parameter("MIN_ALIGNED_BP", 80, int, 1)
TE_GENOME_MAX_EVALUE = te_genome_parameter("MAX_EVALUE", 1e-10, float, 0)
TE_GENOME_MIN_CHAIN_GAP = te_genome_parameter("MIN_CHAIN_GAP", 200, int, 0)
TE_GENOME_MAX_GAP_FRACTION = te_genome_parameter("MAX_GAP_FRACTION", 0.2, float, 0)
TE_GENOME_MAX_OVERLAP = te_genome_parameter("MAX_OVERLAP", 100, int, 0)
TE_GENOME_FULL_LENGTH_COVERAGE = te_genome_parameter(
    "FULL_LENGTH_COVERAGE", 80.0, float, 0, 100
)
TE_GENOME_HIGH_CONFIDENCE_PIDENT = te_genome_parameter(
    "HIGH_CONFIDENCE_PIDENT", 80.0, float, 0, 100
)
TE_GENOME_PARTIAL_COVERAGE = te_genome_parameter(
    "PARTIAL_COVERAGE", 20.0, float, 0, 100
)
TE_GENOME_AMBIGUITY_OVERLAP = te_genome_parameter(
    "AMBIGUITY_OVERLAP", 0.8, float, 0.000001, 1
)
if TE_GENOME_PARTIAL_COVERAGE > TE_GENOME_FULL_LENGTH_COVERAGE:
    raise ValueError(
        "PARAMS.TE_GENOME.PARTIAL_COVERAGE must not exceed "
        "FULL_LENGTH_COVERAGE"
    )

TE_GENOME_FINAL_OUTPUTS = []
for target_name in TE_GENOME_TARGETS:
    TE_GENOME_FINAL_OUTPUTS.extend(
        [
            TE_GENOME_FRAGMENTS.format(target=target_name),
            TE_GENOME_MATCHES.format(target=target_name),
            TE_GENOME_COPIES.format(target=target_name),
            TE_GENOME_RELATIONS.format(target=target_name),
            TE_GENOME_BED.format(target=target_name),
            TE_GENOME_GFF.format(target=target_name),
        ]
    )


rule te_genome:
    input:
        TE_GENOME_FINAL_OUTPUTS,
        TE_GENOME_PUBLIC_BED if "GENOME" in TE_GENOME_TARGETS else [],


# Historical named target retained as a command-line alias during the rename.
rule insider_all_te:
    input:
        TE_GENOME_FINAL_OUTPUTS,
        TE_GENOME_PUBLIC_BED if "GENOME" in TE_GENOME_TARGETS else [],


rule make_te_genome_blast_database:
    input:
        fasta=te_genome_target_fasta,
    output:
        nhr=temp(TE_GENOME_BLAST_NHR),
        nin=temp(TE_GENOME_BLAST_NIN),
        nsq=temp(TE_GENOME_BLAST_NSQ),
    params:
        prefix=f"{TE_GENOME_ROOT}/{{target}}/BLAST_DB/target",
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/te_genome_{{target}}_makeblastdb.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_genome_{{target}}_makeblastdb.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {TE_GENOME_ROOT}/{wildcards.target}/BLAST_DB {WORKDIR}/log {WORKDIR}/benchmarks
        makeblastdb -in {input.fasta:q} -dbtype nucl -out {params.prefix:q} > {log:q} 2>&1
        """


rule blast_te_against_genome_target:
    input:
        query=PREPARED_TE_DATABASE,
        nhr=TE_GENOME_BLAST_NHR,
        nin=TE_GENOME_BLAST_NIN,
        nsq=TE_GENOME_BLAST_NSQ,
    output:
        blast=temp(TE_GENOME_BLAST),
    params:
        database=f"{TE_GENOME_ROOT}/{{target}}/BLAST_DB/target",
        max_evalue=TE_GENOME_MAX_EVALUE,
    threads: THREADS
    resources:
        mem_mb=4096,
    log:
        f"{WORKDIR}/log/te_genome_{{target}}_blast.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_genome_{{target}}_blast.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {TE_GENOME_ROOT}/{wildcards.target} {WORKDIR}/log {WORKDIR}/benchmarks
        blastn -task blastn -num_threads {threads} -evalue {params.max_evalue} \
            -query {input.query:q} -db {params.database:q} \
            -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' \
            -out {output.blast:q} > {log:q} 2>&1
        """


rule annotate_te_genome_target:
    input:
        script=str(TE_GENOME_ANNOTATOR),
        blast=TE_GENOME_BLAST,
        te_headers=f"{WORKDIR}/INPUT/te_headers.tsv",
    output:
        fragments=TE_GENOME_FRAGMENTS,
        matches=TE_GENOME_MATCHES,
        copies=TE_GENOME_COPIES,
        relations=TE_GENOME_RELATIONS,
        bed=TE_GENOME_BED,
        gff=TE_GENOME_GFF,
    params:
        chrom_regex=TE_GENOME_CHROM_REGEX,
        min_pident=TE_GENOME_MIN_PIDENT,
        min_aligned_bp=TE_GENOME_MIN_ALIGNED_BP,
        max_evalue=TE_GENOME_MAX_EVALUE,
        min_chain_gap=TE_GENOME_MIN_CHAIN_GAP,
        max_gap_fraction=TE_GENOME_MAX_GAP_FRACTION,
        max_overlap=TE_GENOME_MAX_OVERLAP,
        full_length_coverage=TE_GENOME_FULL_LENGTH_COVERAGE,
        high_confidence_pident=TE_GENOME_HIGH_CONFIDENCE_PIDENT,
        partial_coverage=TE_GENOME_PARTIAL_COVERAGE,
        ambiguity_overlap=TE_GENOME_AMBIGUITY_OVERLAP,
    threads: 1
    resources:
        mem_mb=8192,
    log:
        f"{WORKDIR}/log/te_genome_{{target}}_annotate.log",
    benchmark:
        f"{WORKDIR}/benchmarks/te_genome_{{target}}_annotate.tsv",
    shell:
        """
        set -euo pipefail
        python3 {input.script:q} {input.blast:q} {input.te_headers:q} \
            --target {wildcards.target:q} --chrom-regex {params.chrom_regex:q} \
            --min-pident {params.min_pident} --min-aligned-bp {params.min_aligned_bp} \
            --max-evalue {params.max_evalue} --min-chain-gap {params.min_chain_gap} \
            --max-gap-fraction {params.max_gap_fraction} --max-overlap {params.max_overlap} \
            --full-length-coverage {params.full_length_coverage} \
            --high-confidence-pident {params.high_confidence_pident} \
            --partial-coverage {params.partial_coverage} \
            --ambiguity-overlap {params.ambiguity_overlap} \
            --fragments {output.fragments:q} --matches {output.matches:q} \
            --copies {output.copies:q} --relations {output.relations:q} \
            --bed {output.bed:q} --gff {output.gff:q} > {log:q} 2>&1
        """


rule publish_te_genome_bed:
    input:
        bed=TE_GENOME_BED.format(target="GENOME"),
    output:
        public_bed=TE_GENOME_PUBLIC_BED,
    threads: 1
    resources:
        mem_mb=128,
    shell:
        "ln -s -r -f {input.bed:q} {output.public_bed:q}"
