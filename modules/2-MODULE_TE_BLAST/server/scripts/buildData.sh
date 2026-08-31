#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 WORK_DIRECTORY [OUTPUT_DIRECTORY] [structure thresholds]" >&2
    exit 2
fi

WORK_DIRECTORY="$1"
shift
OUTPUT_DIRECTORY=""
if [[ $# -gt 0 && "$1" != -* ]]; then
    OUTPUT_DIRECTORY="$1"
    shift
fi

WORK_DIRECTORY="$(cd "$WORK_DIRECTORY" && pwd)"
if [[ -z "$OUTPUT_DIRECTORY" ]]; then
    OUTPUT_DIRECTORY="$WORK_DIRECTORY/MODULE_TE_BLAST"
fi
mkdir -p "$OUTPUT_DIRECTORY"
OUTPUT_DIRECTORY="$(cd "$OUTPUT_DIRECTORY" && pwd)"

SCRIPT_DIRECTORY="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="$(cd "$SCRIPT_DIRECTORY/../../../.." && pwd)"
BLAST="$WORK_DIRECTORY/OUTSIDER/TrEMOLO_SV_TE/INS/SV_INS_CLUST.bln"
TE_SIZES="$WORK_DIRECTORY/1-UTILS/TE_SIZE.tsv"
QUERY_SIZES="$WORK_DIRECTORY/OUTSIDER/TrEMOLO_SV_TE/INS/SV_SIZE.tsv"

for required in "$BLAST" "$TE_SIZES" "$QUERY_SIZES"; do
    if [[ ! -f "$required" ]]; then
        echo "ERROR: required module input not found: $required" >&2
        exit 2
    fi
done

OPTIONAL_INPUTS=()
if [[ -s "$WORK_DIRECTORY/TE_INFOS.bed" ]]; then
    OPTIONAL_INPUTS+=(--te-infos "$WORK_DIRECTORY/TE_INFOS.bed")
fi
if [[ -s "$WORK_DIRECTORY/TE_CALL_CANDIDATES.tsv" ]]; then
    OPTIONAL_INPUTS+=(--te-call-candidates "$WORK_DIRECTORY/TE_CALL_CANDIDATES.tsv")
fi

python3 "$PIPELINE_ROOT/lib/python/modules/build_insertion_structure_report.py" \
    --blast "$BLAST" \
    --te-sizes "$TE_SIZES" \
    --query-sizes "$QUERY_SIZES" \
    --output "$OUTPUT_DIRECTORY" \
    "${OPTIONAL_INPUTS[@]}" \
    "$@"

(
    cd "$OUTPUT_DIRECTORY"
    env XDG_CACHE_HOME=.quarto-cache DENO_DIR=.quarto-cache/deno \
        quarto render report.qmd --output report.html
)
cp -f "$OUTPUT_DIRECTORY/report.html" "$OUTPUT_DIRECTORY/index.html"

echo "Insertion structure report: $OUTPUT_DIRECTORY/report.html"
