#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

WORKDIR="$1"
PIPELINE_DIR="$2"
GENOME="$3"
CHROM_KEEP="$4"
LOG_PATH="$5"

trap 'status=$?; printf "%s\n" "[ERROR][build_report_assets] line=$LINENO cmd=$BASH_COMMAND" | tee -a "${LOG_PATH}.err"; exit $status' ERR

build_tsd_chart() {
    local values_file="$1"
    local js_name="$2"
    local template_name="$3"

    if test -s "$values_file"; then
        python3 "$PIPELINE_DIR/lib/python/graphs/conv_js_histo_grouped_ggplot.py" \
            "$values_file" data -c "$WORKDIR/.tmp_orderc.txt" \
            > "$WORKDIR/REPORT/mini_report/js/$js_name"
        sed "s/%%VAR_JS%%/$js_name/g" \
            "$WORKDIR/REPORT/mini_report/lib/template_TSD_CHART.html" \
            > "$WORKDIR/REPORT/mini_report/lib/$template_name"
    fi
}

build_groups_js() {
    local scope="$1"
    local output_name="$2"

    if [ "$scope" = "INOUTSIDER" ]; then
        awk 'NR>1 {print $15}' "$WORKDIR/TE_INFOS.bed" | sort -u | \
            awk 'BEGIN{chain="var typesSV = [";} { if(NR==1){chain=chain"\"" $1 "\""}else{chain=chain", \"" $1 "\""} } END{chain=chain"]"; print chain; }' \
            > "$WORKDIR/js/$output_name"
    else
        awk -v scope="$scope" '$0 ~ scope {print $15}' "$WORKDIR/TE_INFOS.bed" | sort -u | \
            awk 'BEGIN{chain="var typesSV = [";} { if(NR==1){chain=chain"\"" $1 "\""}else{chain=chain", \"" $1 "\""} } END{chain=chain"]"; print chain; }' \
            > "$WORKDIR/js/$output_name"
    fi
}

build_strand_js() {
    local scope="$1"
    local output_name="$2"
    local template_name="$3"
    local presence_file="$4"

    if [ "$scope" = "INOUTSIDER" ]; then
        test -s "$WORKDIR/TE_INFOS.bed" && \
            awk '{
                split($4, sp, "|");

                if($5 == "-"){
                    dico["-"][sp[1]]+=1;
                    if(! dico["+"][sp[1]]){
                        dico["+"][sp[1]] = 0;
                    }
                }
                else{
                    dico["+"][sp[1]] += 1;
                    if(! dico["-"][sp[1]]){
                        dico["-"][sp[1]] = 0;
                    }
                }
            }

            END{
                chain = "var data = ["
                for(TE in dico["-"]){
                    chain = chain"\n{TEname: \""TE"\", reverse: -"dico["-"][TE]", forward: "dico["+"][TE]"},"
                }
                print substr(chain, 1, length(chain)-1)"\n]";
            }' "$WORKDIR/TE_INFOS.bed" > "$WORKDIR/REPORT/mini_report/js/$output_name"
    else
        test -s "$WORKDIR/TE_INFOS.bed" && test -s "$presence_file" && \
            awk -v scope="$scope" '$0 ~ scope {
                split($4, sp, "|");

                if($5 == "-"){
                    dico["-"][sp[1]]+=1;
                    if(! dico["+"][sp[1]]){
                        dico["+"][sp[1]] = 0;
                    }
                }
                else{
                    dico["+"][sp[1]] += 1;
                    if(! dico["-"][sp[1]]){
                        dico["-"][sp[1]] = 0;
                    }
                }
            }

            END{
                chain = "var data = ["
                for(TE in dico["-"]){
                    chain = chain"\n{TEname: \""TE"\", reverse: -"dico["-"][TE]", forward: "dico["+"][TE]"},"
                }
                print substr(chain, 1, length(chain)-1)"\n]";
            }' "$WORKDIR/TE_INFOS.bed" > "$WORKDIR/REPORT/mini_report/js/$output_name"
    fi

    sed "s/%%VAR_JS%%/$output_name/g" \
        "$WORKDIR/REPORT/mini_report/lib/template_COUNT_TE_SENS_ANTISENS.html" \
        > "$WORKDIR/REPORT/mini_report/lib/$template_name"
}

echo "BUILD REPORT ASSETS..."
cp -r "$PIPELINE_DIR/report/"* "$WORKDIR/REPORT/"

echo "  BUILD TSD CHART..."
printf "TSD OK\nTSD KO\n" > "$WORKDIR/.tmp_orderc.txt"

build_tsd_chart "$WORKDIR/VALUES_TSD_GROUP_OUTSIDER.csv" "TSD_OUTSIDER.js" "template_TSD_OUTSIDER_CHART.html" &
pid_tsd_outsider=$!
build_tsd_chart "$WORKDIR/VALUES_TSD_INSIDER_GROUP.csv" "TSD_INSIDER.js" "template_TSD_INSIDER_CHART.html" &
pid_tsd_insider=$!
build_tsd_chart "$WORKDIR/VALUES_TSD_ALL_GROUP.csv" "TSD_INOUTSIDER.js" "template_TSD_INOUTSIDER_CHART.html" &
pid_tsd_inout=$!

wait "$pid_tsd_outsider"
wait "$pid_tsd_insider"
wait "$pid_tsd_inout"

if ! test -s "$WORKDIR/REPORT/mini_report/lib/template_TSD_INSIDER_CHART.html"; then
    sed -i '/TSD_INSIDER_histo/d' "$WORKDIR/REPORT/mini_report/insider.Rmd"
fi

if ! test -s "$WORKDIR/REPORT/mini_report/lib/template_TSD_INOUTSIDER_CHART.html"; then
    sed -i '/TSD_ALL_histo/d' "$WORKDIR/REPORT/mini_report/inoutsider.Rmd"
fi

rm -f "$WORKDIR/.tmp_orderc.txt"

samtools faidx "$GENOME"
mkdir -p "$WORKDIR/js"

echo "  BUILD chroms for scatter FREQUENCY js ..."
cut -f 1 "$WORKDIR/TE_INFOS.bed" | sort -u > "$WORKDIR/tmp_chrom_with_TE.txt"
grep -E "$(echo "$CHROM_KEEP" | tr "," "|")" "$GENOME.fai" | \
    grep -w -f "$WORKDIR/tmp_chrom_with_TE.txt" 2> /dev/null | \
    awk 'BEGIN{print "var chroms = {"} {print "   \""$1"\" : "$2","} END{print "}"}' \
    > "$WORKDIR/js/chroms.js"

echo "  BUILD group for scatter FREQUENCY js ..."
build_groups_js "OUTSIDER" "groups_OUTSIDER.js" &
pid_groups_out=$!
build_groups_js "INSIDER" "groups_INSIDER.js" &
pid_groups_in=$!
build_groups_js "INOUTSIDER" "groups_INOUTSIDER.js" &
pid_groups_all=$!

wait "$pid_groups_out"
wait "$pid_groups_in"
wait "$pid_groups_all"

echo "  STRAND OUTSIDER"
build_strand_js "OUTSIDER" "TE_COUNT_SENS_ANTISENS_OUTSIDER.js" "template_COUNT_TE_SENS_ANTISENS_OUTSIDER.html" "$WORKDIR/POSITION_TE_OUTSIDER.bed" &
pid_strand_out=$!
echo "  STRAND INSIDER"
build_strand_js "INSIDER" "TE_COUNT_SENS_ANTISENS_INSIDER.js" "template_COUNT_TE_SENS_ANTISENS_INSIDER.html" "$WORKDIR/POSITION_TE_INSIDER.bed" &
pid_strand_in=$!
echo "  STRAND INOUTSIDER"
build_strand_js "INOUTSIDER" "TE_COUNT_SENS_ANTISENS_INOUTSIDER.js" "template_COUNT_TE_SENS_ANTISENS_INOUTSIDER.html" "" &
pid_strand_all=$!

wait "$pid_strand_out"
wait "$pid_strand_in"
wait "$pid_strand_all"

cp "$WORKDIR"/js/* "$WORKDIR/REPORT/mini_report/js/"
rm -fr "$WORKDIR/js"
