#!/bin/bash

#REQUIRE 
#bedtools v2.30.0
#minimap2 2.17-r941


## bash correction_genome_transloc_diff_chrom.sh -c chrom.txt -r "^[23][LR]_|^[4]_|^[X]_" G0-F100_chrom.fasta corrected_genome_${e}.fasta corrected_genome_${i}.fasta
#DEFAULT VALUE
CHROM_COMMUN_FILE="NONE"
CHROM_KEPPING_REGEX="."
SIZE_MIN_TO_REPLACE=20000

#TODO a spprimer
#CHROM_KEPPING_REGEX="^[23][LR]_|^[4]_|^[X]_"


function help(){
    echo "Usage : $0 [options] <reference> <query> <output>"
    echo "Options:"
    echo "  -c, --chrom        file list chromosome commun between reference and query, chrom must be separate by \":\" exemple chr1_Homo:chr1_Sapiens "
    echo "  -r, --chrom_regex  regex chrom to keep [.]"
    echo "  -s, --size_min     size minimum contig to replace [20000]"
    echo "  -h, --help         print help message (this)"
    echo -e "\nWarning : if the names of the chromosomes are not the same between the reference and the genome we recomended put option -c"
}


POSITIONAL_ARGS=()

if [[ ! $# -gt 0 ]]; then
    echo "No parameter found"
    help
    exit 1
fi

while [[ $# -gt 0 ]]; do
  case $1 in
    -c|--chrom)
      CHROM_COMMUN_FILE="$2"
      shift
      shift
      ;;
    -r|--chrom_regex)
      CHROM_KEPPING_REGEX="$2"
      shift
      shift
      ;;
    -s|--size_min)
      SIZE_MIN_TO_REPLACE="$2"
      shift
      shift 
      ;;
    -h|--help)
      help
      exit 0
      ;;
    # --default)
    #   DEFAULT=YES
    #   shift
    #   ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

echo "ARGS:"
echo "  CHROM_COMMUN_FILE    = ${CHROM_COMMUN_FILE}"
echo "  CHROM_KEPPING_REGEX  = ${CHROM_KEPPING_REGEX}"
echo "  SIZE_MIN_TO_REPLACE  = ${SIZE_MIN_TO_REPLACE}"
echo


if [[ ! -n "$1" ]] || [[ ! -n "$2" ]] || [[ ! -n "$3" ]]; then
    echo "Argument required"
    help
    exit 1
fi;


REF=$1
GENOME_FILE=$2
OUTPUT_GENOME=$3



rm -f *.fa*.fai

PAF_FILE=$4
if [[ $PAF_FILE = "" ]]; then
    echo "NO PAF FILE";
    echo "MAPPING...";
    minimap2 -x asm20 -d ${REF}.mmi ${REF}
    minimap2 -x asm20 ${REF} ${GENOME_FILE} > mapping.paf 2> log.err ;
    PAF_FILE="mapping.paf";
    rm -f ${REF}.mmi;
fi;

if [[ "${CHROM_COMMUN_FILE}" == "NONE" ]]; then
    echo "DETERMINATION CHROM COMMUN";
    for i in `cat tmp_chrom.txt`; do 
        c1=`echo $i | cut -d ":" -f 1`; 
        c2=`echo $i | cut -d ":" -f 2`; 
        awk -v c1="$c1" -v c2="$c2" ' $1 == c1 && c2 == $6 {print $0}' mapping.paf | awk 'OFS="\t" { if($2<$3){print $1, $2, $3, $4;}else{print $1, $3, $2, $4;} }' | bedtools sort | bedtools merge | \
            awk -v c1="$c1" -v c2="$c2" '{s+=$3-$2} END{print $1"\t"s"\t"c1":"c2}'; \
    done | sort -k 2 -n -r | awk '{if(! dico[$1]){dico[$1]=1; print $3}}' > chrom.txt;

    CHROM_COMMUN_FILE="chrom.txt";
fi

#GET CHROM QUERY SIZE $1 $2
awk '{print $1":"$2}' ${PAF_FILE} | sort | uniq | tr ":" "\t" | \
    awk 'OFS="\t"{print $1, 0, $2}' | bedtools sort | awk 'OFS="\t"{print $1, $3}' | grep -E $CHROM_KEPPING_REGEX > size_chrom_query.txt

#GET CHROM REF SIZE $6 $7
awk '{print $6":"$7}' ${PAF_FILE} | sort | uniq | tr ":" "\t" | \
    awk 'OFS="\t"{print $1, 0, $2}' | bedtools sort | awk 'OFS="\t"{print $1, $3}' | grep -E $CHROM_KEPPING_REGEX > size_chrom_reference.txt

#ref chrom similary
#REF -> QUERY
echo "GET CHROM similary REF AND QUERY";

if [[ "${CHROM_COMMUN_FILE}" != "NONE" ]] && test -s ${CHROM_COMMUN_FILE} ; then
    awk 'NR == FNR {  split($0, sp, ":"); chrom_commun[sp[1]]=sp[2]; chrom_commun[sp[2]]=sp[1]; next; } chrom_commun[$1] == $6 && OFS="\t"{print $6, $8, $9, $1, $3, $4}' \
        ${CHROM_COMMUN_FILE} ${PAF_FILE}  | bedtools sort | grep -E $CHROM_KEPPING_REGEX  > refSim_sorted.bed
else
    awk 'substr($1, 1, 4) == substr($6, 1, 4) && OFS="\t"{print $6, $8, $9, $1, $3, $4}' \
        ${PAF_FILE}  | bedtools sort | grep -E $CHROM_KEPPING_REGEX  > refSim_sorted.bed
fi;

bedtools merge -i refSim_sorted.bed | grep -E $CHROM_KEPPING_REGEX > refSim_merged_sorted.bed

#DELETION reference
#REF
echo "GET DELETION REF...";
bedtools complement -i refSim_merged_sorted.bed -g size_chrom_reference.txt | awk 'OFS="\t"{print $1, $2, $3, $3-$2}' \
 | grep -E $CHROM_KEPPING_REGEX > deletion_refSim.bed 

#coord ref chrom no similary on query
#REF -> QUERY
echo "GET CHROM DIFFERENT to QUERY...";

if [[ "${CHROM_COMMUN_FILE}" != "NONE" ]] && test -s ${CHROM_COMMUN_FILE} ; then
    awk ' NR == FNR {  split($0, sp, ":"); chrom_commun[sp[1]]=sp[2]; chrom_commun[sp[2]]=sp[1]; next; } chrom_commun[$1] != $6 && OFS="\t"{print $6, $8, $9, $1, $3, $4}' ${CHROM_COMMUN_FILE} ${PAF_FILE} | \
        bedtools sort \
        | grep -E $CHROM_KEPPING_REGEX > ref_sorted_diff_query.bed
else
    awk 'substr($1, 1, 4) != substr($6, 1, 4) && OFS="\t"{print $6, $8, $9, $1, $3, $4}' ${PAF_FILE} | \
        bedtools sort \
        | grep -E $CHROM_KEPPING_REGEX > ref_sorted_diff_query.bed
fi;

#GET QUERY TO REPLACE
#REF -> QUERY
echo "GET QUERY TO REPLACE...";
bedtools intersect -a ref_sorted_diff_query.bed -b deletion_refSim.bed -wo \
    | awk 'OFS="\t"{print $7, $8, $9, $10, $4, $5, $6}' > query_to_replace.bed


cat size_chrom_reference.txt size_chrom_query.txt > size_chrom.txt

#deletion_query_reference sans merge
awk 'BEGIN{
    chrom=""; 
    chrom_2=""
    pred_end_query=0; 
} 

NR == FNR {  
    dico[$1] = $2
    next;
}


OFS="\t" { 

    if(chrom == $1){

        if($2<pred_end_query){
            debq=$2; 
            endq=pred_end_query;
        } 
        else{
            debq=pred_end_query; 
            endq=$2;
        }


        if($5<pred_end_ref){
            debr=$5; 
            endr=pred_end_ref;
        } 
        else{
            debr=pred_end_ref; 
            endr=$5;
        }


        print $1, debq, endq, endq-debq, $4, debr, endr, endr-debr
        pred_end_ref=$6;
        pred_end_query=$3; 
    }
    else if(chrom != $1){
        if(dico[chrom]<pred_end_query){
            debq=dico[chrom]; 
            endq=pred_end_query;
        } 
        else{
            debq=pred_end_query; 
            endq=dico[chrom];
        }


        if(dico[chrom_2]<pred_end_ref){
            debr=dico[chrom_2]; 
            endr=pred_end_ref;
        } 
        else{
            debr=pred_end_ref; 
            endr=dico[chrom_2];
        }


        if(debq < dico[chrom]){
            print chrom, debq, dico[chrom], dico[chrom]-debq, chrom_2, debr, dico[chrom_2], dico[chrom_2]-debr
        }
        pred_end_ref=$6;
        pred_end_query=$3; 
        chrom=$1
        chrom_2=$4
    }
    else{
        pred_end_ref=$6;
        pred_end_query=$3; 
    }
} 

END{

    if(debq < dico[chrom]){
        print chrom, debq, dico[chrom], dico[chrom]-debq, chrom_2, debr, dico[chrom_2], dico[chrom_2]-debr
    }
}' size_chrom.txt refSim_sorted.bed > deletion_query_reference.bed
#REF -> QUERY
#deletion on query des sequence de la ref



echo "GET OVERLAP BETWEEN QUERY TO REPLACE AND DELETION REF MATCH ON QUERY"
bedtools intersect -a query_to_replace.bed -b deletion_query_reference.bed -wo | cut -f 5-7,12-15 \
    | awk 'OFS="\t"{print $1, $2, $3, $3-$2, $4, $5, $6, $7}' > tmp_pos_query_replace.bed



#suprime les multi mappeur
#TODO ? commentaire
echo "REMOVE MULTI MAPPEUR"
cat tmp_pos_query_replace.bed | bedtools sort | bedtools cluster > tmp_cluster.bed

cat tmp_cluster.bed | cut -f 9 | sort | uniq -c | awk '$1==1 {print $2"$"}' > ID_UNIQ.txt
cat tmp_cluster.bed | cut -f 9 | sort | uniq -c | awk '$1>=2 {print $2"$"}' | sort | uniq > ID_OTHER.txt

grep -w -f ID_OTHER.txt tmp_cluster.bed > tmp_pos_query_replace_doublons.bed


#conserve seulement le maximum parmi les multi mappeur

#G31 Probléme

# 2R_RaGOO_RaGOO  23757005    23773028    16023   X_RaGOO_RaGOO   20779166    20783476    4310    99
# 2R_RaGOO_RaGOO  23770012    27022820    3252808 3R_RaGOO_RaGOO  19087982    19088101    119 99
# 2R_RaGOO_RaGOO  23770012    27022820    3252808 3R_RaGOO_RaGOO  19091299    19093267    1968    99

rm -f add_pos.bed
for id in `cat ID_OTHER.txt`;
do
    #Query query
    # chrom1=`grep -w $id tmp_pos_query_replace_doublons.bed | head -n 1 | cut -f 1`
    # chrom2=`grep -w $id tmp_pos_query_replace_doublons.bed | head -n 1 | cut -f 5`
    
    # nb1=`grep -w $id tmp_pos_query_replace_doublons.bed | grep -w "^${chrom1}" -c`
    # nb2=`grep -w $id tmp_pos_query_replace_doublons.bed | cut -f 5 | grep -w "${chrom2}"  -c`

    ID=`grep -w $id tmp_pos_query_replace_doublons.bed | cut -f 9 | head -n 1`
    
    # if [[ $nb1 -eq $nb2 ]]; then
        max1=`grep -w "$id" tmp_pos_query_replace_doublons.bed | awk 'BEGIN{max=0}{if(max<$4){max=$4}}END{print max}'`
        max2=`grep -w "$id" tmp_pos_query_replace_doublons.bed | awk -v max1="$max1" 'BEGIN{ max=0 }{ if(max<$8 && max1==$4){ max=$8 } }END{print max}'`

        nb_max=`grep -w "$id" tmp_pos_query_replace_doublons.bed | awk -v max1="$max1" -v max2="$max2" '$4 == max1 && $8 == max2' | wc -l` 

        if [[ $nb_max -eq 1 ]]; then
            grep -w "$id" tmp_pos_query_replace_doublons.bed | awk -v max1="$max1" -v max2="$max2" '$4 == max1 && $8 == max2' >> add_pos.bed
        else
            nb_max=`grep -w "$id" tmp_pos_query_replace_doublons.bed | awk -v max1="$max1" '$4 == max1' | wc -l` 
            if [[ $nb_max -eq 1 ]]; then
                grep -w "$id" tmp_pos_query_replace_doublons.bed | awk -v max1="$max1" '$4 == max1' >> add_pos.bed
            else
                echo "   ERROR TWO MAX_1 FOUND : id=$id  max1=$max1 ";
            fi;
            echo "ERROR TWO MAX FOUND : id=$id  max1=$max1  max2=$max2";
        fi
    # fi
done;


grep -w -f ID_UNIQ.txt tmp_cluster.bed > pos_query_replace.bed
cat add_pos.bed pos_query_replace.bed | sort -k 9 -n > tmp_pos.bed 
cat tmp_pos.bed | awk -v SIZE_MIN="$SIZE_MIN_TO_REPLACE" '$4>=SIZE_MIN' > pos_query_replace.bed
rm -f tmp_pos.bed

##KEEP GOOD CHROM
grep -E $CHROM_KEPPING_REGEX size_chrom_query.txt > size_chrom_keep.txt
grep -v -E $CHROM_KEPPING_REGEX size_chrom_query.txt > size_chrom_exlude.txt


rm -f ${OUTPUT_GENOME} && touch ${OUTPUT_GENOME}
for chrom_size in `cat size_chrom_keep.txt | tr "\t" ":"`
do
    chrom=`echo $chrom_size | cut -d ":" -f 1` # "_RaGOO_RaGOO_RaGOO"
    #echo $chrom
    size=`echo $chrom_size | cut -d ":" -f 2`

    echo "OLD SIZE CHROMOSOME=$chrom ; SIZE=$size"

    #Exlusion des positions situer sur d'autre chromosome
    # Ancienne position, permet de supprimé les mauvaises insertioin  sur ce chromosome
    # Validé
    grep -w "^${chrom}" pos_query_replace.bed | awk 'OFS="\t"{print $1, $2, $3, $4}' \
        | bedtools sort | bedtools merge > exclude_oldPosition_${chrom}.bed
    #cat tmp_exlude_${chrom}.bed | bedtools sort | bedtools merge > exlude_${chrom}.bed

    #echo "2"
    awk -v chrom="$chrom" '$5==chrom' pos_query_replace.bed > pos_${chrom}.bed
    awk 'OFS="\t"{print $5, $6, $7, $8, $1, $2, $3, $4}' pos_${chrom}.bed | bedtools sort > pos_${chrom}_rev.bed

    #Exlusion des nouvelles positions 
    #echo "3"
    #awk 'OFS="\t"{print $5, $6, $7, $8}' pos_${chrom}.bed > exlude_newPosition_${chrom}.bed
    #cat exlude_newPosition_${chrom}.bed exclude_oldPosition_${chrom}.bed | cut -f 1-3 | bedtools sort | bedtools merge > exlude_${chrom}.bed
    cat exclude_oldPosition_${chrom}.bed | cut -f 1-3 | bedtools sort | bedtools merge > exlude_${chrom}.bed

    bedtools complement -i exlude_${chrom}.bed -g size_chrom_query.txt | grep "${chrom}" | awk 'OFS="\t"{print $0, $3-$2, $0, $3-$2}' | bedtools sort > complement_${chrom}.bed
    
    cat complement_${chrom}.bed pos_${chrom}_rev.bed | bedtools sort | awk 'OFS="\t"{print $5, $6, $7, $8, $1, $2, $3, $4}' > tmp_chrom_${chrom}.bed

    #GET STRAND
    for i in `cat tmp_chrom_${chrom}.bed | tr "\t" ":"`; do 
        
        tmp_chrom=`echo $i | cut -d ":" -f 1`
        start=`echo $i | cut -d ":" -f 2`;
        end=`echo $i | cut -d ":" -f 3`;

        if [[ $tmp_chrom != $chrom ]]; then
            strand=`grep -w "$tmp_chrom" ${PAF_FILE} | grep -w "$end" | cut -f 5 | head -n 1`
            echo -e $tmp_chrom"\t"$start"\t"$end"\tname""\t1""\t$strand"
        else
            echo -e $tmp_chrom"\t"$start"\t"$end"\tname""\t1""\t+"
        fi;

    done > chrom_${chrom}.bed;

    #echo "5"
    #get morceau sequence
    bedtools getfasta -fi ${GENOME_FILE} -bed chrom_${chrom}.bed -s > tmp_chrom_${chrom}.fasta


    awk -v chrom="$chrom" 'BEGIN{print ">"chrom; chain=""} substr($0, 1, 1) != ">" {chain = chain""$0}END{print chain}' tmp_chrom_${chrom}.fasta > chrom_${chrom}.fasta
    
    size_chrom_fa=`grep -v ">" chrom_${chrom}.fasta | wc -c`
    echo "CHROM=$chrom : NEW SIZE CHROMOSOME == $size_chrom_fa"
    cat chrom_${chrom}.fasta >> ${OUTPUT_GENOME}
done;

##Rajoute les chromosome non traité tel qu'il sont
for chrom_size in `cat size_chrom_exlude.txt | tr "\t" ":"`
do
    chrom=`echo $chrom_size | cut -d ":" -f 1`
    size=`echo $chrom_size | cut -d ":" -f 2`
    echo "chrom=$chrom ; size=$size"

    grep -w ">${chrom}" ${GENOME_FILE} -A 1  >> ${OUTPUT_GENOME}
done;

wc -c ${OUTPUT_GENOME} ${GENOME_FILE}



#bash correction_genome_transloc_diff_chrom.sh G0-F100_chrom.fasta G100_on_G0-F100dm6_chrom.fasta corrected_genome_1.fasta
#---------------------------

# rm -f corrected_genome_*.fasta* 
# e=0;
# GENOME_FILE="G100_on_G0-F100dm6_chrom.fasta"
# cat ${GENOME_FILE} > corrected_genome_${e}.fasta
# for (( i = 1; i < 4; i++ )); do
#     bash correction_genome_transloc_diff_chrom.sh G0-F100_chrom.fasta corrected_genome_${e}.fasta corrected_genome_${i}.fasta
#     echo "$e"  
#     e=$(($e+1));
# done;



# rm -f corrected_genome_*.fasta* 
# e=0;
# GENOME_FILE="G31_chrom.fasta"
# cat ${GENOME_FILE} > corrected_genome_${e}.fasta
# for (( i = 1; i < 4; i++ )); do
#     bash correction_genome_transloc_diff_chrom.sh -c chrom.txt -r "^[23][LR]_|^[4]_|^[X]_" G0-F100_chrom.fasta corrected_genome_${e}.fasta corrected_genome_${i}.fasta
#     echo "$e"  
#     e=$(($e+1));
# done;


# bash correction_genome_transloc_diff_chrom.sh G0-F100_chrom.fasta corrected_genome_5.fasta corrected_genome_6.fasta