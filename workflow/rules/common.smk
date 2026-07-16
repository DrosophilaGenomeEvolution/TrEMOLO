rule prepare_inputs:
    input:
        genome=GENOME,
        te_database=TE_DATABASE,
        reference=REFERENCE if REFERENCE else [],
        sample=SAMPLE if SAMPLE else [],
    output:
        done=f"{WORKDIR}/.inputs.prepared",
        manifest=f"{WORKDIR}/input_manifest.json",
        genome=PREPARED_GENOME,
        reference=PREPARED_REFERENCE if REFERENCE else [],
        sample=PREPARED_SAMPLE if SAMPLE else [],
        te_database=f"{WORKDIR}/INPUT/te_database.fasta",
        te_headers=f"{WORKDIR}/INPUT/te_headers.tsv",
    params:
        output=WORKDIR,
        script=str(PREPARE_INPUTS_SCRIPT),
        reference=(f"--reference {REFERENCE}" if REFERENCE else ""),
        sample=(f"--sample {SAMPLE}" if SAMPLE else ""),
    log:
        f"{WORKDIR}/log/prepare_inputs.log",
    shell:
        """
        mkdir -p {WORKDIR}/log
        python3 {params.script:q} \
            --genome {input.genome:q} \
            --te-database {input.te_database:q} \
            {params.reference} \
            {params.sample} \
            --output {params.output:q} > {log:q} 2>&1
        """
