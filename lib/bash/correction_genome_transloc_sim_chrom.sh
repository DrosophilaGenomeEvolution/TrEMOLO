#!/usr/bin/env bash

# Reference-guided correction of inversions and ordering errors occurring
# within one chromosome. Every query base is emitted exactly once.

set -euo pipefail

SCRIPT_DIRECTORY=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=correction_genome_common.sh
source "${SCRIPT_DIRECTORY}/correction_genome_common.sh"

CHROMOSOME_FILE="NONE"
CHROMOSOME_REGEX="."
MINIMUM_SIZE=20000
MINIMUM_MAPQ=20
MINIMUM_IDENTITY=80
MAXIMUM_ANCHOR_OVERLAP=1000
THREADS=1
MINIMAP_PRESET="asm20"
AUDIT_PREFIX=""
KEEP_TEMP=false
REPORT_ONLY=false

help() {
    cat <<'EOF'
Usage: correction_genome_transloc_sim_chrom.sh [options] <reference> <query> <output> [paf-file]

Reorder and reverse query blocks within their original chromosome according
to non-overlapping reference alignments. Unaligned sequence and chromosomes
without accepted alignments are retained.

Options:
  -c, --chrom FILE              Query/reference chromosome pairs (query:reference)
  -r, --chrom_regex REGEX       Only use query chromosomes matching REGEX [.]
  -s, --size_min INT            Minimum query, reference and block span [20000]
      --min-anchor-size INT     Alias for --size_min
      --min-mapq INT            Minimum PAF mapping quality [20]
      --min-identity FLOAT      Minimum alignment identity in percent [80]
      --max-anchor-overlap INT  Maximum overlap between retained anchors [1000]
  -t, --threads INT             minimap2 threads when no PAF is supplied [1]
      --preset NAME             minimap2 preset [asm20]
      --audit-prefix PATH       Prefix for TSV audit files [<output>.correction]
      --report-only             Build and audit the correction plan only
      --keep-temp               Keep the isolated temporary directory
  -h, --help                    Show this help

The optional PAF must describe: reference as target, query as query.
EOF
}

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -c|--chrom)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            CHROMOSOME_FILE=$2
            shift 2
            ;;
        -r|--chrom_regex)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            CHROMOSOME_REGEX=$2
            shift 2
            ;;
        -s|--size_min|--min-anchor-size)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            MINIMUM_SIZE=$2
            shift 2
            ;;
        --min-mapq)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            MINIMUM_MAPQ=$2
            shift 2
            ;;
        --min-identity)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            MINIMUM_IDENTITY=$2
            shift 2
            ;;
        --max-anchor-overlap)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            MAXIMUM_ANCHOR_OVERLAP=$2
            shift 2
            ;;
        -t|--threads)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            THREADS=$2
            shift 2
            ;;
        --preset)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            MINIMAP_PRESET=$2
            shift 2
            ;;
        --audit-prefix)
            [[ $# -ge 2 ]] || cg_die "missing value for $1"
            AUDIT_PREFIX=$2
            shift 2
            ;;
        --report-only)
            REPORT_ONLY=true
            shift
            ;;
        --keep-temp)
            KEEP_TEMP=true
            shift
            ;;
        -h|--help)
            help
            exit 0
            ;;
        --)
            shift
            POSITIONAL_ARGS+=("$@")
            break
            ;;
        -*)
            cg_die "unknown option: $1"
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done

set -- "${POSITIONAL_ARGS[@]}"
[[ $# -ge 3 && $# -le 4 ]] || {
    help >&2
    cg_die "expected 3 or 4 positional arguments"
}

REFERENCE_FASTA=$1
QUERY_FASTA=$2
OUTPUT_FASTA=$3
PAF_FILE=${4:-}
[[ "$OUTPUT_FASTA" != "$QUERY_FASTA" ]] \
    || cg_die "output FASTA must differ from query FASTA"

cg_require_file "$REFERENCE_FASTA"
cg_require_file "$QUERY_FASTA"
[[ -z "$PAF_FILE" ]] || cg_require_file "$PAF_FILE"
cg_validate_positive_integer "minimum anchor size" "$MINIMUM_SIZE"
cg_validate_nonnegative_integer "minimum MAPQ" "$MINIMUM_MAPQ"
cg_validate_percentage "minimum identity" "$MINIMUM_IDENTITY"
cg_validate_nonnegative_integer "maximum anchor overlap" "$MAXIMUM_ANCHOR_OVERLAP"
cg_validate_positive_integer "threads" "$THREADS"

for command_name in awk bedtools cmp cut ln mktemp readlink samtools sort; do
    cg_require_command "$command_name"
done
[[ -n "$PAF_FILE" ]] || cg_require_command minimap2

mkdir -p -- "$(dirname -- "$OUTPUT_FASTA")"
if [[ -z "$AUDIT_PREFIX" ]]; then
    AUDIT_PREFIX="${OUTPUT_FASTA}.correction"
fi
mkdir -p -- "$(dirname -- "$AUDIT_PREFIX")"

TEMPORARY_DIRECTORY=$(mktemp -d "${TMPDIR:-/tmp}/tremolo-same-chrom.XXXXXX")
cleanup() {
    if [[ "$KEEP_TEMP" == true ]]; then
        printf 'Temporary files retained in %s\n' "$TEMPORARY_DIRECTORY" >&2
    else
        rm -rf -- "$TEMPORARY_DIRECTORY"
    fi
}
trap cleanup EXIT

QUERY_LINK="${TEMPORARY_DIRECTORY}/query.fasta"
REFERENCE_LINK="${TEMPORARY_DIRECTORY}/reference.fasta"
RAW_PAF="${TEMPORARY_DIRECTORY}/mapping.paf"
FILTERED_PAF="${TEMPORARY_DIRECTORY}/mapping.filtered.paf"
PAIR_FILE="${AUDIT_PREFIX}.chromosomes.tsv"
ANCHOR_AUDIT="${AUDIT_PREFIX}.anchors.tsv"
SELECTED_ANCHORS="${TEMPORARY_DIRECTORY}/anchors.selected.tsv"
SELECTION_AUDIT="${AUDIT_PREFIX}.selected_anchors.tsv"
PLAN_FILE="${TEMPORARY_DIRECTORY}/reconstruction.plan.tsv"
CORRECTION_AUDIT="${AUDIT_PREFIX}.corrections.tsv"

cg_link_and_index_fasta "$QUERY_FASTA" "$QUERY_LINK"
cg_link_and_index_fasta "$REFERENCE_FASTA" "$REFERENCE_LINK"

if [[ -n "$PAF_FILE" ]]; then
    cp -- "$PAF_FILE" "$RAW_PAF"
else
    printf 'Aligning query against reference with minimap2...\n' >&2
    minimap2 -x "$MINIMAP_PRESET" --secondary=no -t "$THREADS" \
        "$REFERENCE_LINK" "$QUERY_LINK" > "$RAW_PAF" \
        2> "${AUDIT_PREFIX}.minimap2.log"
fi

cg_filter_paf "$RAW_PAF" "$FILTERED_PAF" "$ANCHOR_AUDIT" \
    "$MINIMUM_SIZE" "$MINIMUM_MAPQ" "$MINIMUM_IDENTITY" "$CHROMOSOME_REGEX"
cg_build_chromosome_pairs "$FILTERED_PAF" "${QUERY_LINK}.fai" \
    "${REFERENCE_LINK}.fai" "$CHROMOSOME_FILE" "$PAIR_FILE" \
    "$TEMPORARY_DIRECTORY"
cg_select_nonoverlapping_anchors "$FILTERED_PAF" "$PAIR_FILE" same \
    "$MAXIMUM_ANCHOR_OVERLAP" "$SELECTED_ANCHORS" "$SELECTION_AUDIT" \
    "$TEMPORARY_DIRECTORY"

# Partition every source chromosome at midpoints between adjacent anchors.
# This assigns each unaligned base to exactly one neighbouring reference block.
awk '
    BEGIN {FS=OFS="\t"}
    FILENAME == ARGV[1] {
        order[++chromosome_count]=$1
        chromosome_size[$1]=$2
        next
    }
    FILENAME == ARGV[2] {home[$1]=$2; next}
    FILENAME == ARGV[3] {
        chrom=$1
        idx=++anchor_count[chrom]
        query_start[chrom,idx]=$2
        query_end[chrom,idx]=$3
        strand[chrom,idx]=$4
        reference_chrom[chrom,idx]=$5
        reference_start[chrom,idx]=$6
        reference_end[chrom,idx]=$7
        matches[chrom,idx]=$8
        alignment_block[chrom,idx]=$9
        mapq[chrom,idx]=$10
        identity[chrom,idx]=$11
        next
    }
    END {
        for (chromosome_index=1; chromosome_index<=chromosome_count; chromosome_index++) {
            chrom=order[chromosome_index]
            count=anchor_count[chrom]
            if (count == 0) {
                print chrom,0,chromosome_size[chrom],chrom,".",0,0,"+",1,1,\
                    "KEEP",0,".","."
                continue
            }

            for (idx=1; idx<=count; idx++) {
                rank=1
                for (other=1; other<=count; other++) {
                    if (reference_start[chrom,other] < reference_start[chrom,idx] || \
                        (reference_start[chrom,other] == reference_start[chrom,idx] && \
                         query_start[chrom,other] < query_start[chrom,idx])) rank++
                }
                reference_rank[idx]=rank
            }

            for (idx=1; idx<=count; idx++) {
                left=(idx == 1 ? 0 : int((query_end[chrom,idx-1] + \
                    query_start[chrom,idx])/2))
                right=(idx == count ? chromosome_size[chrom] : \
                    int((query_end[chrom,idx] + query_start[chrom,idx+1])/2))
                if (right <= left) {
                    print "invalid partition for " chrom ": " left "-" right \
                        > "/dev/stderr"
                    exit 2
                }

                moved=(reference_rank[idx] != idx)
                reversed=(strand[chrom,idx] == "-")
                if (moved && reversed) action="MOVE_AND_REVERSE"
                else if (moved) action="MOVE"
                else if (reversed) action="REVERSE"
                else action="KEEP"

                print chrom,left,right,chrom,reference_chrom[chrom,idx],\
                    reference_start[chrom,idx],reference_end[chrom,idx],\
                    strand[chrom,idx],idx,reference_rank[idx],action,\
                    alignment_block[chrom,idx],identity[chrom,idx],mapq[chrom,idx]
            }
        }
    }
' "${QUERY_LINK}.fai" "$PAIR_FILE" "$SELECTED_ANCHORS" > "$PLAN_FILE" \
    || cg_die "failed to build intrachromosomal reconstruction plan"

cg_validate_plan "$PLAN_FILE" "${QUERY_LINK}.fai"
cg_write_correction_audit "$PLAN_FILE" "$CORRECTION_AUDIT"

if [[ "$REPORT_ONLY" == true ]]; then
    printf 'Correction plan written to %s (report-only mode).\n' \
        "$CORRECTION_AUDIT"
    exit 0
fi

cg_reconstruct_fasta "$QUERY_LINK" "${QUERY_LINK}.fai" "$PLAN_FILE" \
    "$OUTPUT_FASTA" "$TEMPORARY_DIRECTORY"
cg_validate_output_fasta same "${QUERY_LINK}.fai" "$OUTPUT_FASTA" \
    "$TEMPORARY_DIRECTORY"

EDITED_BLOCKS=$(awk 'BEGIN {FS="\t"} $11 != "KEEP" {count++} END {print count+0}' \
    "$PLAN_FILE")
printf 'Corrected FASTA: %s (%s edited block(s)).\n' "$OUTPUT_FASTA" \
    "$EDITED_BLOCKS"
printf 'Audit files: %s.{anchors,selected_anchors,chromosomes,corrections}.tsv\n' \
    "$AUDIT_PREFIX"
