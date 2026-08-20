REPORT_DIR = f"{WORKDIR}/REPORT"
REPORT_DATA = f"{REPORT_DIR}/report-data.json"
REPORT_SOURCE = f"{REPORT_DIR}/report.qmd"
REPORT_HTML = f"{REPORT_DIR}/report.html"
REPORT_OPTIONS = config.get("REPORT", {})
REPORT_TITLE = REPORT_OPTIONS.get("TITLE", "TrEMOLO analysis report")
REPORT_AUTHOR = REPORT_OPTIONS.get("AUTHOR", "")
REPORT_LOCUS_WINDOW = REPORT_OPTIONS.get("LOCUS_WINDOW", 100)
REPORT_QUARTO = config.get("TOOLS", {}).get("QUARTO", "quarto")
REPORT_TEMPLATE = PIPELINE_ROOT / "report/quarto/report.qmd"
REPORT_STYLE = PIPELINE_ROOT / "report/quarto/report.css"
REPORT_SCRIPT = PIPELINE_ROOT / "report/quarto/report.js"
REPORT_BUILDER = PIPELINE_ROOT / "lib/python/workflow/build_quarto_report.py"
REPORT_INPUT_MANIFEST = f"{WORKDIR}/input_manifest.json"
REPORT_LOG_DIR = f"{WORKDIR}/log"

if not isinstance(REPORT_LOCUS_WINDOW, int) or REPORT_LOCUS_WINDOW < 0:
    raise ValueError("REPORT.LOCUS_WINDOW must be a non-negative integer")


def report_optional_input(path, enabled):
    return path if enabled else "/dev/null"


REPORT_ENABLED_SOURCES = []
if TE_INFOS_WITH_INSIDER:
    REPORT_ENABLED_SOURCES.append("INSIDER")
if TE_INFOS_WITH_OUTSIDER:
    REPORT_ENABLED_SOURCES.append("OUTSIDER")
REPORT_SOURCE_ARGS = " ".join(
    "--enabled-source {}".format(source) for source in REPORT_ENABLED_SOURCES
)


rule report:
    input:
        REPORT_HTML,


rule prepare_quarto_report:
    input:
        te_infos=TE_INFOS,
        genome_index=PREPARED_GENOME_INDEX,
        manifest=REPORT_INPUT_MANIFEST,
        mapping_stats=report_optional_input(MAPPING_STATS, TE_INFOS_WITH_OUTSIDER),
        sv_vcf=report_optional_input(SV_VCF, TE_INFOS_WITH_OUTSIDER),
        builder=str(REPORT_BUILDER),
        template=str(REPORT_TEMPLATE),
        style=str(REPORT_STYLE),
        script=str(REPORT_SCRIPT),
    output:
        data=REPORT_DATA,
        source=REPORT_SOURCE,
    params:
        title=REPORT_TITLE,
        author=REPORT_AUTHOR,
        work_directory=WORKDIR,
        source_args=REPORT_SOURCE_ARGS,
        locus_window=REPORT_LOCUS_WINDOW,
        log_dir=REPORT_LOG_DIR,
    log:
        f"{WORKDIR}/log/prepare_quarto_report.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.work_directory:q}/REPORT {params.log_dir:q}
        python3 {input.builder:q} \
            --te-infos {input.te_infos:q} \
            --genome-index {input.genome_index:q} \
            --input-manifest {input.manifest:q} \
            --mapping-stats {input.mapping_stats:q} \
            --sv-vcf {input.sv_vcf:q} \
            --template {input.template:q} \
            --style {input.style:q} \
            --script {input.script:q} \
            --output-qmd {output.source:q} \
            --output-json {output.data:q} \
            --title {params.title:q} \
            --author {params.author:q} \
            --work-directory {params.work_directory:q} \
            --locus-window {params.locus_window} \
            {params.source_args} > {log:q} 2>&1
        """


rule render_quarto_report:
    input:
        source=REPORT_SOURCE,
        data=REPORT_DATA,
    output:
        html=REPORT_HTML,
    params:
        report_dir=REPORT_DIR,
        quarto=REPORT_QUARTO,
        log_dir=REPORT_LOG_DIR,
    log:
        f"{WORKDIR}/log/render_quarto_report.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.report_dir:q}/.quarto-cache {params.log_dir:q}
        if ! command -v {params.quarto:q} >/dev/null 2>&1; then
            echo "Quarto CLI was not found: {params.quarto}" >&2
            echo "Rebuild the TrEMOLO container or render REPORT/report.qmd with Quarto on the host." >&2
            exit 127
        fi
        (
            cd {params.report_dir:q}
            rm -f .report.tmp.html
            XDG_CACHE_HOME=.quarto-cache \
            DENO_DIR=.quarto-cache/deno \
                {params.quarto:q} render report.qmd --output .report.tmp.html
            test -s .report.tmp.html
            mv -f .report.tmp.html report.html
        ) > {log:q} 2>&1
        test -s {output.html:q}
        """
