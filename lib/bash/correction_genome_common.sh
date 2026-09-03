#!/usr/bin/env bash

# Shared implementation for the two reference-guided assembly correction
# scripts.  The public entry points remain Bash scripts; awk is used for the
# tabular PAF and block planning steps.

cg_die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

cg_warn() {
    printf 'WARNING: %s\n' "$*" >&2
}

cg_require_command() {
    command -v "$1" >/dev/null 2>&1 || cg_die "required command not found: $1"
}

cg_require_file() {
    [[ -s "$1" ]] || cg_die "file not found or empty: $1"
}

cg_validate_positive_integer() {
    [[ "$2" =~ ^[0-9]+$ && "$2" -gt 0 ]] \
        || cg_die "$1 must be a positive integer (got: $2)"
}

cg_validate_nonnegative_integer() {
    [[ "$2" =~ ^[0-9]+$ ]] \
        || cg_die "$1 must be a non-negative integer (got: $2)"
}

cg_validate_percentage() {
    awk -v value="$2" 'BEGIN {exit !(value ~ /^[0-9]+([.][0-9]+)?$/ && value >= 0 && value <= 100)}' \
        || cg_die "$1 must be a number between 0 and 100 (got: $2)"
}

cg_link_and_index_fasta() {
    local source_fasta=$1
    local temporary_fasta=$2
    local absolute_source

    absolute_source=$(readlink -f -- "$source_fasta")
    [[ -n "$absolute_source" ]] || cg_die "cannot resolve FASTA path: $source_fasta"
    ln -s -- "$absolute_source" "$temporary_fasta"
    samtools faidx "$temporary_fasta"
    [[ -s "${temporary_fasta}.fai" ]] \
        || cg_die "failed to index FASTA: $source_fasta"
}

cg_filter_paf() {
    local input_paf=$1
    local filtered_paf=$2
    local audit_file=$3
    local minimum_size=$4
    local minimum_mapq=$5
    local minimum_identity=$6
    local chromosome_regex=$7

    awk \
        -v accepted_file="$filtered_paf" \
        -v minimum_size="$minimum_size" \
        -v minimum_mapq="$minimum_mapq" \
        -v minimum_identity="$minimum_identity" \
        -v chromosome_regex="$chromosome_regex" '
        BEGIN {
            FS=OFS="\t"
            print "query_chrom", "query_start", "query_end", "strand", \
                "reference_chrom", "reference_start", "reference_end", \
                "matches", "alignment_block", "mapq", "identity", \
                "status", "reason"
        }
        {
            status="accepted"
            reason="accepted"
            identity=0
            query_span=0
            reference_span=0
            secondary=0
            if (NF < 12) {
                status="rejected"
                reason="malformed_paf"
            } else {
                query_span=$4-$3
                reference_span=$9-$8
                # PAF matches/alignment-block is strongly depressed by large
                # structural gaps and is therefore not a sequence-identity
                # estimate for assembly alignments. Prefer the minimap2 `dv`
                # divergence tag (then `de`) and keep the column ratio only as
                # a compatibility fallback for tag-less PAF producers.
                identity=($11 > 0 ? 100*$10/$11 : 0)
                dv=""
                de=""
                for (field=13; field<=NF; field++) {
                    if ($field == "tp:A:S") secondary=1
                    if ($field ~ /^dv:f:/) dv=substr($field,6)
                    if ($field ~ /^de:f:/) de=substr($field,6)
                }
                number_regex="^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"
                if (dv ~ number_regex && dv >= 0 && dv <= 1) {
                    identity=100*(1-dv)
                } else if (de ~ number_regex && de >= 0 && de <= 1) {
                    identity=100*(1-de)
                }
                if ($1 !~ chromosome_regex) {
                    status="rejected"
                    reason="chromosome_filter"
                } else if (secondary) {
                    status="rejected"
                    reason="secondary_alignment"
                } else if (query_span < minimum_size) {
                    status="rejected"
                    reason="query_span_below_minimum"
                } else if (reference_span < minimum_size) {
                    status="rejected"
                    reason="reference_span_below_minimum"
                } else if ($11 < minimum_size) {
                    status="rejected"
                    reason="alignment_block_below_minimum"
                } else if ($12 < minimum_mapq) {
                    status="rejected"
                    reason="mapq_below_minimum"
                } else if (identity < minimum_identity) {
                    status="rejected"
                    reason="identity_below_minimum"
                }
            }
            print (NF >= 1 ? $1 : "."), (NF >= 3 ? $3 : "."), \
                (NF >= 4 ? $4 : "."), (NF >= 5 ? $5 : "."), \
                (NF >= 6 ? $6 : "."), (NF >= 8 ? $8 : "."), \
                (NF >= 9 ? $9 : "."), (NF >= 10 ? $10 : "."), \
                (NF >= 11 ? $11 : "."), (NF >= 12 ? $12 : "."), \
                sprintf("%.4f", identity), status, reason
            if (status == "accepted") print $0 > accepted_file
        }
    ' "$input_paf" > "$audit_file"

    touch "$filtered_paf"
}

cg_build_chromosome_pairs() {
    local filtered_paf=$1
    local query_index=$2
    local reference_index=$3
    local chromosome_file=$4
    local output_pairs=$5
    local temporary_directory=$6
    local scores_file="${temporary_directory}/chromosome_pair_scores.tsv"

    awk 'BEGIN {FS=OFS="\t"} {score[$1 FS $6]+=$11} END {
        for (key in score) print key, score[key]
    }' "$filtered_paf" | sort -t $'\t' -k1,1 -k3,3nr -k2,2 > "$scores_file"

    if [[ "$chromosome_file" != "NONE" ]]; then
        cg_require_file "$chromosome_file"
        awk '
            BEGIN {FS=OFS="\t"; errors=0}
            FILENAME == ARGV[1] {query[$1]=1; next}
            FILENAME == ARGV[2] {reference[$1]=1; next}
            FILENAME == ARGV[3] {score[$1 FS $2]=$3; next}
            FILENAME == ARGV[4] {
                line=$0
                sub(/^[[:space:]]+/, "", line)
                sub(/[[:space:]]+$/, "", line)
                if (line == "" || line ~ /^#/) next
                count=split(line, values, /[\t:]/)
                if (count != 2 || values[1] == "" || values[2] == "") {
                    print "invalid chromosome pair: " $0 > "/dev/stderr"
                    errors++
                    next
                }
                if (query[values[1]] && reference[values[2]]) {
                    q=values[1]; r=values[2]
                } else if (query[values[2]] && reference[values[1]]) {
                    q=values[2]; r=values[1]
                } else {
                    print "chromosome pair does not match query/reference FASTA: " \
                        $0 > "/dev/stderr"
                    errors++
                    next
                }
                if (seen[q] && pair[q] != r) {
                    print "query chromosome has two reference pairs: " q \
                        > "/dev/stderr"
                    errors++
                    next
                }
                if (reference_seen[r] && reference_pair[r] != q) {
                    print "reference chromosome has two query pairs: " r \
                        > "/dev/stderr"
                    errors++
                    next
                }
                seen[q]=1
                pair[q]=r
                reference_seen[r]=1
                reference_pair[r]=q
                order[++number]=q
            }
            END {
                for (idx=1; idx<=number; idx++) {
                    q=order[idx]
                    if (printed[q]++) continue
                    value=score[q FS pair[q]]
                    print q, pair[q], (value == "" ? 0 : value)
                }
                if (errors) exit 2
            }
        ' "$query_index" "$reference_index" "$scores_file" \
            "$chromosome_file" > "$output_pairs" \
            || cg_die "invalid chromosome correspondence file: $chromosome_file"
    else
        # Infer a deterministic one-to-one correspondence from the strongest
        # cumulative alignment scores. Queries left without a unique target
        # are deliberately not corrected and remain unchanged in the output.
        sort -t $'\t' -k3,3nr -k1,1 -k2,2 "$scores_file" | \
            awk 'BEGIN {FS=OFS="\t"}
                !query_seen[$1] && !reference_seen[$2] {
                    query_seen[$1]=1
                    reference_seen[$2]=1
                    print $1,$2,$3
                }
            ' | sort -t $'\t' -k1,1 > "$output_pairs"
    fi
}

cg_select_nonoverlapping_anchors() {
    local filtered_paf=$1
    local chromosome_pairs=$2
    local mode=$3
    local maximum_overlap=$4
    local selected_file=$5
    local selection_audit=$6
    local temporary_directory=$7
    local candidates="${temporary_directory}/anchor_candidates.tsv"
    local guarded_home_anchors="${temporary_directory}/anchor_guarded_home.tsv"
    local prioritized_candidates="${temporary_directory}/anchor_candidates.prioritized.tsv"
    local ranked="${temporary_directory}/anchor_candidates.ranked.tsv"
    local eligible_queries="${temporary_directory}/anchor_eligible_queries.tsv"

    # A query is eligible only if at least one accepted alignment supports its
    # own query/reference chromosome correspondence. This guarantees that an
    # interchromosomal plan never empties an output chromosome merely because
    # it was supported only by a cross-chromosome alignment.
    awk '
        BEGIN {FS=OFS="\t"}
        FILENAME == ARGV[1] {home[$1]=$2; next}
        home[$1] == $6 {eligible[$1]=1}
        END {for (query in eligible) print query}
    ' "$chromosome_pairs" "$filtered_paf" | sort > "$eligible_queries"

    awk -v mode="$mode" '
        BEGIN {FS=OFS="\t"}
        FILENAME == ARGV[1] {home[$1]=$2; owner[$2]=$1; next}
        FILENAME == ARGV[2] {eligible[$1]=1; next}
        {
            if (!($1 in eligible)) next
            if (mode == "same" && home[$1] != $6) next
            if (mode != "same" && !($6 in owner)) next
            identity=($11 > 0 ? 100*$10/$11 : 0)
            dv=""
            de=""
            for (field=13; field<=NF; field++) {
                if ($field ~ /^dv:f:/) dv=substr($field,6)
                if ($field ~ /^de:f:/) de=substr($field,6)
            }
            number_regex="^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$"
            if (dv ~ number_regex && dv >= 0 && dv <= 1) {
                identity=100*(1-dv)
            } else if (de ~ number_regex && de >= 0 && de <= 1) {
                identity=100*(1-de)
            }
            print $1,$3,$4,$5,$6,$8,$9,$10,$11,$12,\
                sprintf("%.4f",identity),(home[$1] == $6 ? 1 : 0)
        }
    ' "$chromosome_pairs" "$eligible_queries" "$filtered_paf" > "$candidates"

    # Protect only the strongest home anchor for each query. Protecting every
    # home alignment would allow a short repetitive hit to reject a much
    # larger cross-chromosome anchor, hiding the translocation this script is
    # meant to correct. One guarded home anchor per one-to-one chromosome pair
    # is enough to keep every destination represented.
    awk 'BEGIN {FS=OFS="\t"} $12 == 1 {print}' "$candidates" | \
        sort -t $'\t' -k1,1 -k9,9nr -k10,10nr -k11,11nr -k2,2n | \
        awk 'BEGIN {FS=OFS="\t"} !seen[$1]++ {
            print $1,$2,$3,$4,$5,$6,$7
        }' > "$guarded_home_anchors"

    awk '
        BEGIN {FS=OFS="\t"}
        FILENAME == ARGV[1] {
            guarded[$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7]=1
            next
        }
        {
            key=$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7
            print $0,(key in guarded ? 1 : 0)
        }
    ' "$guarded_home_anchors" "$candidates" > "$prioritized_candidates"

    sort -t $'\t' -k13,13nr -k9,9nr -k10,10nr -k11,11nr -k1,1 -k2,2n \
        "$prioritized_candidates" > "$ranked"

    awk -v selected_file="$selected_file" -v maximum_overlap="$maximum_overlap" '
        function max(a,b) {return a>b?a:b}
        function min(a,b) {return a<b?a:b}
        BEGIN {
            FS=OFS="\t"
            print "query_chrom", "query_start", "query_end", "strand", \
                "reference_chrom", "reference_start", "reference_end", \
                "alignment_block", "mapq", "identity", "status", "reason"
        }
        {
            rejected=0
            reason="selected"
            for (idx=1; idx<=accepted; idx++) {
                query_overlap=0
                reference_overlap=0
                if ($1 == aq[idx]) {
                    query_overlap=max(0, min($3,aqe[idx])-max($2,aqs[idx]))
                }
                if ($5 == ar[idx]) {
                    reference_overlap=max(0, min($7,are[idx])-max($6,ars[idx]))
                }
                if (query_overlap > maximum_overlap) {
                    rejected=1
                    reason="overlapping_query_anchor"
                    break
                }
                if (reference_overlap > maximum_overlap) {
                    rejected=1
                    reason="overlapping_reference_anchor"
                    break
                }
            }
            if (!rejected) {
                accepted++
                aq[accepted]=$1; aqs[accepted]=$2; aqe[accepted]=$3
                ar[accepted]=$5; ars[accepted]=$6; are[accepted]=$7
                print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 > selected_file
                print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,"selected","selected"
            } else {
                print $1,$2,$3,$4,$5,$6,$7,$9,$10,$11,"rejected",reason
            }
        }
    ' "$ranked" > "$selection_audit"

    touch "$selected_file"
    sort -t $'\t' -k1,1 -k2,2n -k3,3n "$selected_file" \
        > "${selected_file}.sorted"
    mv "${selected_file}.sorted" "$selected_file"
}

cg_validate_plan() {
    local plan_file=$1
    local query_index=$2

    awk '
        BEGIN {FS=OFS="\t"; errors=0}
        FILENAME == ARGV[1] {size[$1]=$2; next}
        {
            chrom=$1; start=$2; end=$3; destination=$4
            if (!(chrom in size)) {
                print "unknown source chromosome in plan: " chrom > "/dev/stderr"
                errors++
                next
            }
            if (!(destination in size)) {
                print "unknown destination chromosome in plan: " destination \
                    > "/dev/stderr"
                errors++
            }
            if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/ || end <= start) {
                print "invalid source interval in plan: " $0 > "/dev/stderr"
                errors++
                next
            }
            key=chrom SUBSEP start SUBSEP end
            if (seen_interval[key]++) {
                print "duplicate source interval in plan: " $0 > "/dev/stderr"
                errors++
            }
            count[chrom]++
            begins[chrom,count[chrom]]=start
            ends[chrom,count[chrom]]=end
            total[chrom]+=end-start
        }
        END {
            for (chrom in size) {
                if (count[chrom] == 0) {
                    print "source chromosome absent from plan: " chrom > "/dev/stderr"
                    errors++
                    continue
                }
                # Plans are emitted in source-coordinate order for each source.
                previous=0
                for (idx=1; idx<=count[chrom]; idx++) {
                    if (begins[chrom,idx] != previous) {
                        print "gap/overlap in source plan for " chrom \
                            " at " previous " -> " begins[chrom,idx] > "/dev/stderr"
                        errors++
                    }
                    previous=ends[chrom,idx]
                }
                if (previous != size[chrom] || total[chrom] != size[chrom]) {
                    print "incomplete source coverage for " chrom \
                        ": " total[chrom] "/" size[chrom] > "/dev/stderr"
                    errors++
                }
            }
            if (errors) exit 2
        }
    ' "$query_index" "$plan_file" || cg_die "invalid reconstruction plan"
}

cg_write_correction_audit() {
    local plan_file=$1
    local output_file=$2

    {
        printf '%s\n' \
            $'source_chrom\tsource_start\tsource_end\tdestination_chrom\treference_chrom\treference_start\treference_end\tstrand\toriginal_order\toutput_order\taction\tanchor_length\tidentity\tmapq'
        cat "$plan_file"
    } > "$output_file"
}

cg_reconstruct_fasta() {
    local query_fasta=$1
    local query_index=$2
    local plan_file=$3
    local output_fasta=$4
    local temporary_directory=$5
    local extraction_bed="${temporary_directory}/extract.bed"
    local extraction_fasta="${temporary_directory}/extract.fasta"

    : > "$output_fasta"
    while IFS=$'\t' read -r chromosome chromosome_size _rest; do
        printf '>%s\n' "$chromosome" >> "$output_fasta"
        awk -v destination="$chromosome" 'BEGIN {FS=OFS="\t"}
            $4 == destination {print $1,$2,$3,$1":"$2"-"$3,"1",$8,$10}
        ' "$plan_file" | sort -t $'\t' -k7,7n | cut -f1-6 > "$extraction_bed"
        [[ -s "$extraction_bed" ]] \
            || cg_die "no sequence block assigned to output chromosome: $chromosome"
        bedtools getfasta -fi "$query_fasta" -bed "$extraction_bed" -s \
            > "$extraction_fasta"
        # bedtools writes every extracted block as an independent FASTA
        # sequence. Concatenating those lines directly would create variable
        # line widths inside one reconstructed chromosome, which samtools
        # faidx correctly rejects. Remove block boundaries as a stream, then
        # let the native fold utility wrap it in linear time and constant
        # memory; a character-by-character awk loop becomes quadratic on
        # chromosome-sized lines emitted by bedtools.
        awk '/^>/ || !NF {next} {printf "%s",$0} END {printf "\n"}' \
            "$extraction_fasta" | fold -w 60 >> "$output_fasta"
    done < "$query_index"
}

cg_validate_output_fasta() {
    local mode=$1
    local query_index=$2
    local output_fasta=$3
    local temporary_directory=$4
    local output_link="${temporary_directory}/corrected.fasta"
    local input_total
    local output_total

    ln -s -- "$(readlink -f -- "$output_fasta")" "$output_link"
    samtools faidx "$output_link"
    input_total=$(awk '{sum+=$2} END {print sum+0}' "$query_index")
    output_total=$(awk '{sum+=$2} END {print sum+0}' "${output_link}.fai")
    [[ "$input_total" -eq "$output_total" ]] \
        || cg_die "sequence length changed during reconstruction: ${input_total} -> ${output_total}"

    if [[ "$mode" == "same" ]]; then
        cut -f1-2 "$query_index" > "${temporary_directory}/input_sizes.tsv"
        cut -f1-2 "${output_link}.fai" > "${temporary_directory}/output_sizes.tsv"
        cmp -s "${temporary_directory}/input_sizes.tsv" \
            "${temporary_directory}/output_sizes.tsv" \
            || cg_die "intrachromosomal correction changed chromosome names or sizes"
    else
        cut -f1 "$query_index" > "${temporary_directory}/input_names.txt"
        cut -f1 "${output_link}.fai" > "${temporary_directory}/output_names.txt"
        cmp -s "${temporary_directory}/input_names.txt" \
            "${temporary_directory}/output_names.txt" \
            || cg_die "interchromosomal correction changed the chromosome set or order"
    fi
}
