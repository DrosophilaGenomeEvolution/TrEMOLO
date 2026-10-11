"""Publish original TE labels after internal scientific calculations."""
import json

TE_NAME_EXPORTS = [
    dict(source=INTERNAL_TE_INFOS, output=TE_INFOS, kind="compound:3"),
    dict(source=INTERNAL_TE_CALL_CANDIDATES, output=TE_CALL_CANDIDATES, kind="table"),
]

def export_root(source):
    TE_NAME_EXPORTS.append(dict(source=source, output=f"{WORKDIR}/{Path(source).name}", kind="compound:3"))

def export_table(source, kind="table"):
    relative = Path(source).relative_to(Path(WORKDIR))
    TE_NAME_EXPORTS.append(dict(source=source, output=str(Path(WORKDIR) / "PUBLIC" / relative), kind=kind))

if TE_INFOS_WITH_INSIDER:
    export_root(INSIDER_POSITION_BED)
    for source in [INSIDER_INSERTION_CSV, INSIDER_DELETION_CSV,
                   INSIDER_INSERTION_COMBINE, INSIDER_DELETION_COMBINE]:
        export_table(source)
    export_table(INSIDER_DELETION_TE_BED, "compound:3")
    export_table(INSIDER_INSERTION_TE_BED, "compound:3")
    export_table(INSIDER_TSD, "compound:0")
if TE_INFOS_WITH_OUTSIDER:
    export_root(OUTSIDER_POSITION_BED)
    for source in [OUTSIDER_SNIFFLES_CSV, OUTSIDER_SNIFFLES_COMBINE, OUTSIDER_INS_CSV,
                   OUTSIDER_INS_COMBINE, OUTSIDER_SOFT_CSV, OUTSIDER_HARD_CSV,
                   OUTSIDER_SOFT_COMBINE, OUTSIDER_HARD_COMBINE]:
        export_table(source)
    export_table(OUTSIDER_TSD, "compound:0")
    export_table(OUTSIDER_FREQUENCY, "column:2")
    export_table(OUTSIDER_FREQUENCY_PRECISE, "column:2")
    for source in [OUTSIDER_INS_COUNT, OUTSIDER_SNIFFLES_COUNT, OUTSIDER_SOFT_COUNT,
                   OUTSIDER_HARD_COUNT, OUTSIDER_MERGED_COUNT]:
        export_table(source, "column:1")
if OUTSIDER_INTEGRATION_ENABLED:
    for source in [OUTSIDER_CANONICAL_PUBLIC_BED, OUTSIDER_OBSERVED_PUBLIC_BED]:
        export_root(source)
    export_table(OUTSIDER_INTEGRATION_AUDIT)
if OUTSIDER_LIFT_ENABLED:
    for source in [OUTSIDER_ON_REFERENCE, POSITION_TE_ON_REFERENCE]:
        export_root(source)
    export_table(OUTSIDER_LIFT_AUDIT)

TE_NAME_PUBLIC_OUTPUTS = [item["output"] for item in TE_NAME_EXPORTS]
TE_NAME_MANIFEST = f"{WORKDIR}/PUBLIC/name-restoration.json"

rule restore_te_names:
    input:
        sources=[item["source"] for item in TE_NAME_EXPORTS],
        mapping=f"{WORKDIR}/INPUT/te_headers.tsv",
        script=str(PIPELINE_ROOT / "lib/python/workflow/restore_te_names.py"),
    output:
        files=TE_NAME_PUBLIC_OUTPUTS,
        manifest=TE_NAME_MANIFEST,
    params:
        exports=json.dumps(TE_NAME_EXPORTS),
    threads: 1
    resources:
        mem_mb=512,
    log:
        f"{WORKDIR}/log/restore_te_names.log",
    benchmark:
        f"{WORKDIR}/benchmarks/restore_te_names.tsv",
    shell:
        "python3 {input.script:q} --te-headers {input.mapping:q} --exports {params.exports:q} "
        "--manifest {output.manifest:q} > {log:q} 2>&1"
