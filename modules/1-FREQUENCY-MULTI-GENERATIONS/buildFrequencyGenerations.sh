#!/usr/bin/env bash
set -euo pipefail

INPUT=""
OUTPUT="SCATTER-FREQ-TE-TrEMOLO"
GENOME=""
CHROM=""
LOCUS_WINDOW=20
TREND_EPSILON=0.001

usage() {
    echo "Usage: $0 -i MANIFEST -g REFERENCE_FASTA [-o OUTPUT] [--locus-window BP] [--trend-epsilon VALUE]"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT="${2:-}"
            shift 2
            ;;
        -o|--output)
            OUTPUT="${2:-}"
            shift 2
            ;;
        -g|--genome)
            GENOME="${2:-}"
            shift 2
            ;;
        -c|--chrom)
            CHROM="${2:-}"
            shift 2
            ;;
        --locus-window)
            LOCUS_WINDOW="${2:-}"
            shift 2
            ;;
        --trend-epsilon)
            TREND_EPSILON="${2:-}"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
    esac
done

if [[ -z "$INPUT" || ! -s "$INPUT" ]]; then
    echo "ERROR: --input must identify a non-empty sample manifest" >&2
    usage >&2
    exit 2
fi
if [[ -z "$GENOME" || ! -s "$GENOME" ]]; then
    echo "ERROR: --genome is required to bind all samples to one reference checksum" >&2
    usage >&2
    exit 2
fi
if [[ -n "$CHROM" ]]; then
    echo "INFO: --chrom is retained for command compatibility; chromosome filtering is now interactive." >&2
fi

MODULE_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_ROOT="$(cd "$MODULE_DIR/../.." && pwd)"
REFERENCE_SHA256="$(sha256sum "$GENOME" | awk '{print $1}')"
REFERENCE_ID="sha256:${REFERENCE_SHA256}"

mkdir -p "$OUTPUT"
OUTPUT="$(cd "$OUTPUT" && pwd)"
python3 "$PIPELINE_ROOT/lib/python/modules/build_population_report.py" \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --reference-id "$REFERENCE_ID" \
    --genome "$GENOME" \
    --locus-window "$LOCUS_WINDOW" \
    --trend-epsilon "$TREND_EPSILON"

(
    cd "$OUTPUT"
    env XDG_CACHE_HOME=.quarto-cache DENO_DIR=.quarto-cache/deno \
        quarto render report.qmd --output report.html
)
cp -f "$OUTPUT/report.html" "$OUTPUT/index.html"

echo "Population report: $OUTPUT/report.html"
