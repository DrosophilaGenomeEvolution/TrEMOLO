"""Build the legacy public TE call table from declared workflow outputs."""

TE_INFOS = f"{WORKDIR}/TE_INFOS.bed"
TE_INFOS_EMPTY_INPUT = "/dev/null"
TE_INFOS_WITH_OUTSIDER = bool(CHOICES.get("OUTSIDER_VARIANT"))
TE_INFOS_WITH_INSIDER = bool(CHOICES.get("INSIDER_VARIANT"))
TE_INFOS_WITH_INSIDER_FREQUENCY = (
    TE_INFOS_WITH_INSIDER and TE_INFOS_WITH_OUTSIDER
)


def te_infos_input(path, enabled):
    return path if enabled else TE_INFOS_EMPTY_INPUT


rule te_infos:
    input:
        TE_INFOS,


rule build_te_infos:
    input:
        script=str(PIPELINE_ROOT / "lib/python/workflow/build_te_infos.py"),
        outsider_positions=te_infos_input(
            OUTSIDER_POSITION_BED, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_sniffles_combine=te_infos_input(
            OUTSIDER_SNIFFLES_COMBINE, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_direct_combine=te_infos_input(
            OUTSIDER_INS_COMBINE, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_soft_calls=te_infos_input(
            OUTSIDER_SOFT_CSV, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_hard_calls=te_infos_input(
            OUTSIDER_HARD_CSV, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_tsd=te_infos_input(OUTSIDER_TSD, TE_INFOS_WITH_OUTSIDER),
        outsider_frequency_precise=te_infos_input(
            OUTSIDER_FREQUENCY_PRECISE, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_frequency=te_infos_input(
            OUTSIDER_FREQUENCY, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_sv_sizes=te_infos_input(
            OUTSIDER_FREQUENCY_SV_SIZES, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_sniffles_calls=te_infos_input(
            OUTSIDER_SNIFFLES_CSV, TE_INFOS_WITH_OUTSIDER
        ),
        outsider_direct_calls=te_infos_input(
            OUTSIDER_INS_CSV, TE_INFOS_WITH_OUTSIDER
        ),
        insider_deletion_bed=te_infos_input(
            INSIDER_DELETION_TE_BED, TE_INFOS_WITH_INSIDER
        ),
        insider_positions=te_infos_input(
            INSIDER_POSITION_BED, TE_INFOS_WITH_INSIDER
        ),
        insider_insertion_combine=te_infos_input(
            INSIDER_INSERTION_COMBINE, TE_INFOS_WITH_INSIDER
        ),
        insider_deletion_combine=te_infos_input(
            INSIDER_DELETION_COMBINE, TE_INFOS_WITH_INSIDER
        ),
        insider_tsd=te_infos_input(INSIDER_TSD, TE_INFOS_WITH_INSIDER),
        insider_frequency=te_infos_input(
            INSIDER_FREQUENCY, TE_INFOS_WITH_INSIDER_FREQUENCY
        ),
        insider_variants=te_infos_input(INSIDER_BED, TE_INFOS_WITH_INSIDER),
    output:
        infos=TE_INFOS,
    threads: 1
    resources:
        mem_mb=2048,
    log:
        f"{WORKDIR}/log/build_te_infos.log",
    benchmark:
        f"{WORKDIR}/benchmarks/build_te_infos.tsv",
    shell:
        """
        set -euo pipefail
        mkdir -p {WORKDIR}/log {WORKDIR}/benchmarks
        python3 {input.script:q} \
            --outsider-positions {input.outsider_positions:q} \
            --outsider-sniffles-combine {input.outsider_sniffles_combine:q} \
            --outsider-direct-combine {input.outsider_direct_combine:q} \
            --outsider-soft-calls {input.outsider_soft_calls:q} \
            --outsider-hard-calls {input.outsider_hard_calls:q} \
            --outsider-tsd {input.outsider_tsd:q} \
            --outsider-frequency-precise {input.outsider_frequency_precise:q} \
            --outsider-frequency {input.outsider_frequency:q} \
            --outsider-sv-sizes {input.outsider_sv_sizes:q} \
            --outsider-sniffles-calls {input.outsider_sniffles_calls:q} \
            --outsider-direct-calls {input.outsider_direct_calls:q} \
            --insider-deletion-bed {input.insider_deletion_bed:q} \
            --insider-positions {input.insider_positions:q} \
            --insider-insertion-combine {input.insider_insertion_combine:q} \
            --insider-deletion-combine {input.insider_deletion_combine:q} \
            --insider-tsd {input.insider_tsd:q} \
            --insider-frequency {input.insider_frequency:q} \
            --insider-variants {input.insider_variants:q} \
            --output {output.infos:q} > {log:q} 2>&1
        """
