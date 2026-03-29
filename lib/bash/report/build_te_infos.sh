#!/usr/bin/env bash

set -Eeuo pipefail
IFS=$'\n\t'

WORKDIR="$1"
PIPELINE_DIR="$2"
LOG_PATH="$3"

trap 'status=$?; printf "%s\n" "[ERROR][build_te_infos] line=$LINENO cmd=$BASH_COMMAND" | tee -a "${LOG_PATH}.err"; exit $status' ERR

warn_msg() {
    printf "%s" "$*"
}

echo "BUILD BED ALL INFOS..."

test -s "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_TE.bed" && \
    cat "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_TE.bed" | bedtools sort | awk 'OFS="\t"{ print $1, $2, $3, $4, $3-$2, $5 }' > "$WORKDIR/POSITION_TE_INSIDER.bed"

test -s "$WORKDIR/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed" && \
    ln -sfr "$WORKDIR/OUTSIDER/TE_DETECTION/MERGE_TE/MERGE_TE_ALL.bed" "$WORKDIR/POSITION_TE_OUTSIDER.bed"

NB_TE_OUT=`(test -s "$WORKDIR/POSITION_TE_OUTSIDER.bed" && cat "$WORKDIR/POSITION_TE_OUTSIDER.bed" | wc -l) || echo 0`
NB_TE_IN=`(test -s "$WORKDIR/POSITION_TE_INSIDER.bed" && cat "$WORKDIR/POSITION_TE_INSIDER.bed" | wc -l) || echo 0`
NB_TE=$(($NB_TE_OUT + $NB_TE_IN))

echo "CHECK INTERSECT TO DELETION INSIDER..."
test -s "$WORKDIR/POSITION_TE_OUTSIDER.bed" && test -s "$WORKDIR/DELETION_TE.bed" && \
    bedtools window -a "$WORKDIR/POSITION_TE_OUTSIDER.bed" -b "$WORKDIR/DELETION_TE.bed" -w 30 | \
    awk '{ split($4, sp1, "|"); split($10, sp2, "|"); if(sp1[1] == sp2[1]){ print sp1[2] } }' > "$WORKDIR/ID_OUTSIDER_INTERSECT_TO_DEL.txt"

test -s "$WORKDIR/POSITION_TE_OUTSIDER.bed" && echo "BUILD TE_INFO OUTSIDER..."

i=0
show_step="$i/$NB_TE TE"
number_of_char=`echo "$show_step" | wc -c`
number_of_char=$(($number_of_char-1))
printf "\b%.0s" `seq 1 $number_of_char`
printf "$show_step"

rm -f "$WORKDIR/TE_INFOS.bed"
echo -e "#chrom\tstart\tend\tTE|ID\tstrand\tTSD\tpident\tpsize_TE\tSIZE_TE\tNEW_POS\tFREQ\tFREQ_WITH_CLIPPED\tSV_SIZE\tID_TrEMOLO\tTYPE" > "$WORKDIR/TE_INFOS.bed"

test -s "$WORKDIR/POSITION_TE_OUTSIDER.bed" && \
    while read line;
do
    read -r chrom start end TE ID strand <<< $(echo $line | awk 'OFS="\t"{split($4, TE, "|"); print $1, $2, $3, TE[1], TE[2], $6}')
    pident=`(grep -w "$ID" "$WORKDIR/OUTSIDER/TE_DETECTION/COMBINE_TE.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/INS/COMBINE_INS_TREMOLO.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/SOFT/SOFT_TE.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.csv") 2> /dev/null | cut -f 3 || echo "NONE"`
    size_per=`(grep -w "$ID" "$WORKDIR/OUTSIDER/TE_DETECTION/COMBINE_TE.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/INS/COMBINE_INS_TREMOLO.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/SOFT/SOFT_TE.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/HARD/HARD_TE.csv") 2> /dev/null | cut -f 4 || echo "NONE"`

    test -s "$WORKDIR/ID_OUTSIDER_INTERSECT_TO_DEL.txt" && \
        intersect_to_del=`grep -w "$ID" "$WORKDIR/ID_OUTSIDER_INTERSECT_TO_DEL.txt" 2> /dev/null || echo ""` || intersect_to_del=""

    TYPE_TrEMOLO=`echo "$ID" | grep -E -o "SOFT|HARD|DEL|INS" || echo "UNDEFINED"`
    if [ -n "$intersect_to_del" ]; then
        TYPE_TrEMOLO="${TYPE_TrEMOLO}_DEL"
    fi

    SIZE_TE=""
    TSD_OK=`grep -w "$ID" "$WORKDIR/OUTSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 2 || echo ""`
    if [ ! -n "$TSD_OK" ]; then
        TSD_OK="NONE"
        SIZE_TE=`grep ":$ID:[0-9]*:[IP]" "$WORKDIR/OUTSIDER/TE_DETECTION/COMBINE_TE.csv" 2> /dev/null | cut -f 7 || echo "NONE"`
    else
        SIZE_TE=`grep -w "$ID" "$WORKDIR/OUTSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 4 | cut -d ":" -f 1 || echo "NONE"`
    fi

    if [ ! -n "$SIZE_TE" ]; then
        SIZE_TE="NONE"
    fi

    NEW_POS=`(test -s "$WORKDIR/OUTSIDER/TSD/TSD_TE.tsv" && grep -w $ID "$WORKDIR/OUTSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 3) || echo "$start"`

    FREQ_BEFOR=""
    FREQ_AFTER=""

    test -s "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv" && \
        FREQ_BEFOR=`grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv" 2> /dev/null | awk '$9<=100{ print $9; } $9>100{ print "NONE"; }' || echo "NONE"` && \
        FREQ_AFTER=`grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS_PRECISE.tsv" 2> /dev/null | awk '$10<=100{ print $10; } $10>100{ print "NONE"; }' || echo "NONE"`

    if [ "$FREQ_BEFOR" == "" ] || [ "$FREQ_BEFOR" == "NONE" ] || [ "$FREQ_AFTER" == "" ] || [ "$FREQ_AFTER" == "NONE" ]; then
        FREQ_BEFOR=`grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv" 2> /dev/null | cut -f 9 || echo "NONE"`
        FREQ_AFTER=`grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv" 2> /dev/null | cut -f 10 || echo "NONE"`
    fi

    ! test -n "$FREQ_BEFOR" && FREQ_BEFOR="NONE-FREQB"
    ! test -n "$FREQ_AFTER" && FREQ_AFTER="NONE-FREQA"

    SV_SIZE=`grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/SV_SIZE.tsv" 2> /dev/null | awk 'max<$2{ max=$2; } END{ if(max!=""){ print max } }' || echo "NONE"`
    ! test -n "$SV_SIZE" && SV_SIZE="NONE"

    IS_CLIPPED=`echo "$ID" | grep -E -o "SOFT|HARD" || echo ""`

    if ! test -n "$IS_CLIPPED"; then
        ID_POS=`echo $line | awk '{print $1":"$2}' | "$PIPELINE_DIR/lib/C++/bin/chain_to_id" 2>/dev/null || echo "0"`
        TYPE=`(grep -w "$ID" "$WORKDIR/OUTSIDER/TE_DETECTION/FILTER_BLAST_SEQUENCE_INDEL_vs_DBTE.csv" 2> /dev/null || grep -w "$ID" "$WORKDIR/OUTSIDER/TrEMOLO_SV_TE/INS/INS_TREMOLO.csv" 2> /dev/null ) | cut -f 2  | cut -d ":" -f 2 | sed 's/[<>]//g'`
        RS=`(test -s "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv" && grep -w "$ID" "$WORKDIR/OUTSIDER/FREQUENCY/FREQUENCY_TE_INS.tsv" 2> /dev/null | cut -f 4 | "$PIPELINE_DIR/lib/C++/bin/chain_to_id" 2>/dev/null) || echo "ERROR_RS_$ID" >> "${LOG_PATH}.err"`
        ! test -n "$RS" && RS=0
        ID_TrEMOLO="TE_ID_OUTSIDER.$ID_POS.$TYPE.$RS"
    else
        ID_TrEMOLO="TE_ID_OUTSIDER.$ID"
    fi

    echo -e "$chrom\t$start\t$end\t$TE|$ID\t$strand\t$TSD_OK\t$pident\t$size_per\t$SIZE_TE\t$NEW_POS\t$FREQ_BEFOR\t$FREQ_AFTER\t$SV_SIZE\t$ID_TrEMOLO\t$TYPE_TrEMOLO" >> "$WORKDIR/TE_INFOS.bed"

    i=$(($i + 1))
    show_step="$i/$NB_TE TE "
    number_of_char=`echo "$show_step" | wc -c`
    number_of_char=$(($number_of_char-1))
    printf "\b%.0s" `seq 1 $number_of_char`
    printf "$show_step"
done < "$WORKDIR/POSITION_TE_OUTSIDER.bed" 2>&2

rm -f "$WORKDIR/tmp.txt"
rm -f "$WORKDIR/ID_OUTSIDER_INTERSECT_TO_DEL.txt"
printf "\b%.0s" `seq 1 10`

awk '
    BEGIN{
        OFS="\t";
        print "x", "y", "condition";
    }

    NR>1 && $6!="NONE" {
        split($4, sp, "|");
        FAMILY[sp[1]]=1;
        TSDOK[sp[1]]+=1;
    }

    NR>1 && $6=="NONE" {
        split($4, sp, "|");
        FAMILY[sp[1]];
        TSDKO[sp[1]]+=1;
    }

    END{
        for(key in FAMILY){
            print key, (TSDOK[key]!="" ? TSDOK[key] : 0), "TSD OK";
            print key, (TSDKO[key]!="" ? TSDKO[key] : 0 ), "TSD KO";
        }
    }' "$WORKDIR/TE_INFOS.bed" > "$WORKDIR/VALUES_TSD_GROUP_OUTSIDER.csv"

echo "..."
test -s "$WORKDIR/POSITION_TE_INSIDER.bed" && echo "BUILD TE_INFO INSIDER..."

num_TE_INSIDER=0
test -s "$WORKDIR/POSITION_TE_INSIDER.bed" && while read line;
do
    ID=`echo "$line" | cut -f 4 | cut -d "|" -f 2`
    TE=`echo "$line" | cut -f 4 | cut -d "|" -f 1`
    chrom=`echo "$line" | cut -f 1`
    start=`echo "$line" | cut -f 2`
    end=`echo "$line" | cut -f 3`

    strand=`echo "$line" | cut -f 6`
    pident=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_COMBINE_TE.csv" 2> /dev/null | cut -f 3 || echo "NONE"`
    size_per=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_COMBINE_TE.csv" 2> /dev/null | cut -f 4 || echo "NONE"`

    SIZE_TE=""
    TSD_OK=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 2 || echo ""`
    if [ ! -n "$TSD_OK" ]; then
        TSD_OK="NONE"
        SIZE_TE=`grep "$ID:[+-]:" "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_COMBINE_TE.csv" 2> /dev/null | cut -f 5 || echo "NONE"`
    else
        SIZE_TE=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 4 || echo ""`
    fi

    if [ ! -n "$SIZE_TE" ]; then
        SIZE_TE="NONE"
    fi

    NEW_POS=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 3 || echo "$start"`
    FREQ_BEFOR=`grep -w "$TE|$ID" "$WORKDIR/INSIDER/FREQ_INSIDER/DEPTH_TE_INSIDER.csv" 2>/dev/null | cut -f 6 || echo "NONE"`
    FREQ_AFTER="INSIDER"
    ID_POS=`echo $line | awk '{print $1":"$2}' | "$PIPELINE_DIR/lib/C++/bin/chain_to_id" 2>/dev/null || echo "0"`
    TYPE=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/INSERTION_COMBINE_TE.csv" 2> /dev/null | cut -f 2 | cut -d ":" -f 3 || echo "UNKNOWN"`

    SV_SIZE=`grep -w "$ID" "$WORKDIR/INSIDER/VARIANT_CALLING/assemblytics_out.Assemblytics_structural_variants.bed" 2> /dev/null | awk '
        $7=="Deletion" || $7=="Repeat_contraction" || $7=="Tandem_contraction" { print $5 }
        $7=="Insertion" || $7=="Repeat_expansion" || $7=="Tandem_expansion" { print $9 }
        ' || echo "NONE"`
    ! test -n "$SV_SIZE" && SV_SIZE="NONE"

    ID_TrEMOLO="TE_ID_INSIDER.$ID_POS.$TYPE"

    echo -e "$chrom\t$start\t$end\t$TE|$ID\t$strand\t$TSD_OK\t$pident\t$size_per\t$SIZE_TE\t$NEW_POS\t$FREQ_BEFOR\t$FREQ_AFTER\t$SV_SIZE\t$ID_TrEMOLO\t$TYPE" >> "$WORKDIR/TE_INFOS.bed"

    i=$(($i + 1))
    show_step="$i/$NB_TE TE"
    number_of_char=`echo "$show_step" | wc -c`
    number_of_char=$(($number_of_char-1))
    printf "\b%.0s" `seq 1 $number_of_char`
    printf "$show_step"

    num_TE_INSIDER=$(($num_TE_INSIDER+1))
done < "$WORKDIR/POSITION_TE_INSIDER.bed"

rm -f "$WORKDIR/ID_INSIDER_INTERSECT_TO_DEL.txt"

num_TE_INSIDER=0
test -s "$WORKDIR/DELETION_TE.bed" && while read line;
do
    IFS=$'\t' read -r chrom start end rest strand <<< "$line"
    IFS='|' read -r TE ID <<< "$rest"

    pident=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/DELETION_COMBINE_TE.csv" 2> /dev/null | cut -f 3 || echo "NONE"`
    size_per=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/DELETION_COMBINE_TE.csv" 2> /dev/null | cut -f 4 || echo "NONE"`

    SIZE_TE=""
    TSD_OK=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 2 || echo ""`
    if [ ! -n "$TSD_OK" ]; then
        TSD_OK="NONE"
        SIZE_TE=`grep "$ID:[+-]:" "$WORKDIR/INSIDER/TE_DETECTION/DELETION_COMBINE_TE.csv" 2> /dev/null | cut -f 5 || echo "NONE"`
    else
        SIZE_TE=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 4 || echo ""`
    fi

    if [ ! -n "$SIZE_TE" ]; then
        SIZE_TE="NONE"
    fi

    NEW_POS=`grep -w "$ID" "$WORKDIR/INSIDER/TSD/TSD_TE.tsv" 2> /dev/null | cut -f 3 || echo "$start"`
    FREQ_BEFOR=`grep -w "$TE|$ID" "$WORKDIR/INSIDER/FREQ_INSIDER/DEPTH_TE_INSIDER.csv" 2>/dev/null | cut -f 6 || echo "NONE"`
    FREQ_AFTER="INSIDER"
    ID_POS=`echo $line | awk '{print $1":"$2}' | "$PIPELINE_DIR/lib/C++/bin/chain_to_id" 2>/dev/null || echo "0"`
    TYPE=`grep -w "$ID" "$WORKDIR/INSIDER/TE_DETECTION/DELETION_COMBINE_TE.csv" 2> /dev/null | cut -f 2 | cut -d ":" -f 3 || echo "UNKNOWN"`

    SV_SIZE=`grep -w "$ID" "$WORKDIR/INSIDER/VARIANT_CALLING/assemblytics_out.Assemblytics_structural_variants.bed" 2> /dev/null | awk '
        $7=="Deletion" || $7=="Repeat_contraction" || $7=="Tandem_contraction" { print $5 }
        $7=="Insertion" || $7=="Repeat_expansion" || $7=="Tandem_expansion" { print $9 }
        ' || echo "NONE"`
    ! test -n "$SV_SIZE" && SV_SIZE="NONE"

    ID_TrEMOLO="TE_ID_INSIDER.$ID_POS.$TYPE"

    echo -e "$chrom\t$start\t$end\t$TE|$ID\t$strand\t$TSD_OK\t$pident\t$size_per\t$SIZE_TE\t$NEW_POS\t$FREQ_BEFOR\t$FREQ_AFTER\t$SV_SIZE\t$ID_TrEMOLO\t$TYPE" >> "$WORKDIR/TE_INFOS.bed"

    i=$(($i + 1))
    show_step="$i/$NB_TE TE"
    number_of_char=`echo "$show_step" | wc -c`
    number_of_char=$(($number_of_char-1))
    printf "\b%.0s" `seq 1 $number_of_char`
    printf "$show_step"

    num_TE_INSIDER=$(($num_TE_INSIDER+1))
done < "$WORKDIR/DELETION_TE.bed"

echo

( test -s "$WORKDIR/POSITION_TE_INSIDER.bed" && test -s "$WORKDIR/TE_INFOS.bed" && \
    grep "_INSIDER" "$WORKDIR/TE_INFOS.bed" | awk '
        BEGIN{
            OFS="\t";
            print "x", "y", "condition";
        }

        $6!="NONE" {
            split($4, sp, "|");
            FAMILY[sp[1]]=1;
            TSDOK[sp[1]]+=1;
        }

        $6=="NONE" {
            split($4, sp, "|");
            FAMILY[sp[1]];
            TSDKO[sp[1]]+=1;
        }

        END{
            for(key in FAMILY){
                print key, (TSDOK[key]!="" ? TSDOK[key] : 0), "TSD OK";
                print key, (TSDKO[key]!="" ? TSDKO[key] : 0 ), "TSD KO";
            }
        }' > "$WORKDIR/VALUES_TSD_INSIDER_GROUP.csv" ) || \
    ( test -s "$WORKDIR/POSITION_TE_INSIDER.bed" && warn_msg "Warning : NO INSIDER IN REPORT\n" ) || echo "..."

test -s "$WORKDIR/TE_INFOS.bed" && \
    awk '
        BEGIN{
            OFS="\t";
            print "x", "y", "condition";
        }

        NR>1 && $6!="NONE" {
            split($4, sp, "|");
            FAMILY[sp[1]]=1;
            TSDOK[sp[1]]+=1;
        }

        NR>1 && $6=="NONE" {
            split($4, sp, "|");
            FAMILY[sp[1]];
            TSDKO[sp[1]]+=1;
        }

        END{
            for(key in FAMILY){
                print key, (TSDOK[key]!="" ? TSDOK[key] : 0), "TSD OK";
                print key, (TSDKO[key]!="" ? TSDKO[key] : 0 ), "TSD KO";
            }
        }' "$WORKDIR/TE_INFOS.bed" > "$WORKDIR/VALUES_TSD_ALL_GROUP.csv" || \
            warn_msg "WARN : NO TE DETECTED..."
