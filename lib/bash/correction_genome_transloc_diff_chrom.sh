#!/usr/bin/env bash

# Reference-guided correction of blocks assigned to the wrong chromosome.
# The same plan can also reorder and reverse blocks inside their destination.
# Every query base is emitted exactly once.

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
THREADS=8
MINIMAP_PRESET="asm20"
AUDIT_PREFIX=""
KEEP_TEMP=false
REPORT_ONLY=false

help() {
    cat <<'EOF'
Usage: correction_genome_transloc_diff_chrom.sh [options] <reference> <query> <output> [paf-file]

Transfer query blocks to the chromosome paired with their reference target,
then order and orient them by non-overlapping reference alignments. Unaligned
sequence and chromosomes without accepted alignments are retained.

Options:
  -c, --chrom FILE              Query/reference chromosome pairs (query:reference)
  -r, --chrom_regex REGEX       Only use query chromosomes matching REGEX [.]
  -s, --size_min INT            Minimum query, reference and block span [20000]
      --min-anchor-size INT     Alias for --size_min
      --min-mapq INT            Minimum PAF mapping quality [20]
      --min-identity FLOAT      Minimum alignment identity in percent [80]
      --max-anchor-overlap INT  Maximum overlap between retained anchors [1000]
  -g, --threads INT             minimap2 threads when no PAF is supplied [8]
  -t INT                        Legacy alias for -g
      --preset NAME             minimap2 preset [asm20]
      --audit-prefix PATH       Prefix for TSV audit files [<output>.correction]
      --report-only             Build and audit the correction plan only
      --keep-temp               Keep the isolated temporary directory
  -h, --help                    Show this help

The optional PAF must describe: reference as target, query as query.
Chromosome correspondence must be one-to-one. Without -c it is inferred from
the strongest cumulative alignment scores; unpaired query contigs are kept.
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
        -g|-t|--threads)
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

for command_name in awk bedtools cmp cut fold ln mktemp readlink samtools sort; do
    cg_require_command "$command_name"
done
[[ -n "$PAF_FILE" ]] || cg_require_command minimap2

mkdir -p -- "$(dirname -- "$OUTPUT_FASTA")"
if [[ -z "$AUDIT_PREFIX" ]]; then
    AUDIT_PREFIX="${OUTPUT_FASTA}.correction"
fi
mkdir -p -- "$(dirname -- "$AUDIT_PREFIX")"

TEMPORARY_DIRECTORY=$(mktemp -d "${TMPDIR:-/tmp}/tremolo-diff-chrom.XXXXXX")
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
GENERATED_PAF="${AUDIT_PREFIX}.mapping.paf"
FILTERED_PAF="${TEMPORARY_DIRECTORY}/mapping.filtered.paf"
PAIR_FILE="${AUDIT_PREFIX}.chromosomes.tsv"
ANCHOR_AUDIT="${AUDIT_PREFIX}.anchors.tsv"
SELECTED_ANCHORS="${TEMPORARY_DIRECTORY}/anchors.selected.tsv"
SELECTION_AUDIT="${AUDIT_PREFIX}.selected_anchors.tsv"
RAW_PLAN="${TEMPORARY_DIRECTORY}/reconstruction.raw.tsv"
ORDER_TABLE="${TEMPORARY_DIRECTORY}/output_order.tsv"
PLAN_FILE="${TEMPORARY_DIRECTORY}/reconstruction.plan.tsv"
CORRECTION_AUDIT="${AUDIT_PREFIX}.corrections.tsv"

cg_link_and_index_fasta "$QUERY_FASTA" "$QUERY_LINK"
cg_link_and_index_fasta "$REFERENCE_FASTA" "$REFERENCE_LINK"

if [[ -n "$PAF_FILE" ]]; then
    cp -- "$PAF_FILE" "$RAW_PAF"
else
    printf 'Aligning query against reference with minimap2...\n' >&2
    minimap2 -x "$MINIMAP_PRESET" --secondary=no -t "$THREADS" \
        "$REFERENCE_LINK" "$QUERY_LINK" > "$GENERATED_PAF" \
        2> "${AUDIT_PREFIX}.minimap2.log"
    RAW_PAF=$GENERATED_PAF
fi

cg_filter_paf "$RAW_PAF" "$FILTERED_PAF" "$ANCHOR_AUDIT" \
    "$MINIMUM_SIZE" "$MINIMUM_MAPQ" "$MINIMUM_IDENTITY" "$CHROMOSOME_REGEX"
cg_build_chromosome_pairs "$FILTERED_PAF" "${QUERY_LINK}.fai" \
    "${REFERENCE_LINK}.fai" "$CHROMOSOME_FILE" "$PAIR_FILE" \
    "$TEMPORARY_DIRECTORY"
cg_select_nonoverlapping_anchors "$FILTERED_PAF" "$PAIR_FILE" all \
    "$MAXIMUM_ANCHOR_OVERLAP" "$SELECTED_ANCHORS" "$SELECTION_AUDIT" \
    "$TEMPORARY_DIRECTORY"

# Partition each source chromosome at midpoints between adjacent selected
# anchors. The reference target determines the destination query chromosome.
awk '
    BEGIN {FS=OFS="\t"}
    FILENAME == ARGV[1] {
        order[++chromosome_count]=$1
        chromosome_size[$1]=$2
        next
    }
    FILENAME == ARGV[2] {
        home[$1]=$2
        owner[$2]=$1
        next
    }
    FILENAME == ARGV[3] {
        chrom=$1
        idx=++anchor_count[chrom]
        query_start[chrom,idx]=$2
        query_end[chrom,idx]=$3
        strand[chrom,idx]=$4
        reference_chrom[chrom,idx]=$5
        reference_start[chrom,idx]=$6
        reference_end[chrom,idx]=$7
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
                print chrom,0,chromosome_size[chrom],chrom,".",0,0,"+",1,0,\
                    "KEEP",0,".","."
                continue
            }

            for (idx=1; idx<=count; idx++) {
                rank=1
                for (other=1; other<=count; other++) {
                    if (reference_chrom[chrom,other] == reference_chrom[chrom,idx] && \
                        (reference_start[chrom,other] < reference_start[chrom,idx] || \
                         (reference_start[chrom,other] == reference_start[chrom,idx] && \
                          query_start[chrom,other] < query_start[chrom,idx]))) rank++
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

                target=reference_chrom[chrom,idx]
                destination=owner[target]
                if (destination == "") {
                    print "no destination query chromosome for reference " target \
                        > "/dev/stderr"
                    exit 2
                }
                transferred=(destination != chrom)
                moved=(!transferred && reference_rank[idx] != idx)
                reversed=(strand[chrom,idx] == "-")
                if (transferred && reversed) action="TRANSFER_AND_REVERSE"
                else if (transferred) action="TRANSFER"
                else if (moved && reversed) action="MOVE_AND_REVERSE"
                else if (moved) action="MOVE"
                else if (reversed) action="REVERSE"
                else action="KEEP"

                print chrom,left,right,destination,target,\
                    reference_start[chrom,idx],reference_end[chrom,idx],\
                    strand[chrom,idx],idx,0,action,alignment_block[chrom,idx],\
                    identity[chrom,idx],mapq[chrom,idx]
            }
        }
    }
' "${QUERY_LINK}.fai" "$PAIR_FILE" "$SELECTED_ANCHORS" > "$RAW_PLAN" \
    || cg_die "failed to build interchromosomal reconstruction plan"

# Assign output order globally for every destination: anchored blocks first in
# reference order, then any chromosome that had no accepted anchor.
awk 'BEGIN {FS=OFS="\t"} {
    anchored=($5 == "." ? 1 : 0)
    reference_position=($5 == "." ? $2 : $6)
    print $4,anchored,reference_position,$1,$2,$3
}' "$RAW_PLAN" | sort -t $'\t' -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n | \
    awk 'BEGIN {FS=OFS="\t"} {
        if ($1 != previous_destination) rank=0
        previous_destination=$1
        rank++
        print $4,$5,$6,rank
    }' > "$ORDER_TABLE"

awk '
    BEGIN {FS=OFS="\t"}
    FILENAME == ARGV[1] {output_order[$1 FS $2 FS $3]=$4; next}
    {
        key=$1 FS $2 FS $3
        if (!(key in output_order)) {
            print "block absent from output-order table: " key > "/dev/stderr"
            errors++
            next
        }
        $10=output_order[key]
        print
    }
    END {if (errors) exit 2}
' "$ORDER_TABLE" "$RAW_PLAN" > "$PLAN_FILE" \
    || cg_die "failed to assign output block order"

cg_validate_plan "$PLAN_FILE" "${QUERY_LINK}.fai"
cg_write_correction_audit "$PLAN_FILE" "$CORRECTION_AUDIT"

if [[ "$REPORT_ONLY" == true ]]; then
    printf 'Correction plan written to %s (report-only mode).\n' \
        "$CORRECTION_AUDIT"
    exit 0
fi

cg_reconstruct_fasta "$QUERY_LINK" "${QUERY_LINK}.fai" "$PLAN_FILE" \
    "$OUTPUT_FASTA" "$TEMPORARY_DIRECTORY"
cg_validate_output_fasta different "${QUERY_LINK}.fai" "$OUTPUT_FASTA" \
    "$TEMPORARY_DIRECTORY"

EDITED_BLOCKS=$(awk 'BEGIN {FS="\t"} $11 != "KEEP" {count++} END {print count+0}' \
    "$PLAN_FILE")
TRANSFERRED_BLOCKS=$(awk 'BEGIN {FS="\t"} $11 ~ /^TRANSFER/ {count++} END {print count+0}' \
    "$PLAN_FILE")
printf 'Corrected FASTA: %s (%s edited block(s), %s transferred).\n' \
    "$OUTPUT_FASTA" "$EDITED_BLOCKS" "$TRANSFERRED_BLOCKS"
printf 'Audit files: %s.{anchors,selected_anchors,chromosomes,corrections}.tsv\n' \
    "$AUDIT_PREFIX"
