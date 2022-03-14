#!/bin/bash

#REQUIRE 
#bedtools v2.30.0
#minimap2 2.17-r941


CHROM_COMMUN_FILE="NONE"
CHROM_KEPPING_REGEX="."

#TODO delete
#CHROM_KEPPING_REGEX="^[23][LR]_|^[4]_|^[X]_"


function help(){
    echo "Usage : $0 [options] <reference> <query> <output> [paf-file]"
    echo "Options:"
    echo "  -c, --chrom        file list chromosome commun between reference and query, chrom must be separate by \":\" exemple chr1_Homo:chr1_Sapiens "
    echo "  -r, --chrom_regex  regex chrom to keep [.]"
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
    -h|--help)
      help
      exit 0
      ;;
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

echo "CHROM_COMMUN_FILE    = ${CHROM_COMMUN_FILE}"
echo "CHROM_KEPPING_REGEX  = ${CHROM_KEPPING_REGEX}"

if [[ ! -n "$1" ]] || [[ ! -n "$2" ]] || [[ ! -n "$3" ]]; then
    echo "Argument required"
    help
    exit 1
fi;


REF=$1
GENOME_FILE=$2
OUTPUT_GENOME=$3
PAF_FILE_TMP=$4

rm -f *.fai
if [[ $PAF_FILE_TMP = "" ]]; then
    echo "NO PAF FILE"
    echo "MAPPING..."
    minimap2 -x asm20 -d ${REF}.mmi ${REF}
    minimap2 -x asm20 ${REF} ${GENOME_FILE} > mapping.paf 2> log.err ;

    PAF_FILE_TMP="mapping.paf"
fi

border_diag_size=15
#reduce=$(($border_diag_size/($nb_tour)-2))

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


#PAF_FILE_TMP="map_3R_G31_to_3R_G0-F100.paf"
#PAF_FILE_TMP="test.paf"

awk '{print $1":"$2}' ${PAF_FILE_TMP} | sort | uniq | tr ":" "\t" | \
    awk 'OFS="\t"{print $1, 0, $2}' | bedtools sort | awk 'OFS="\t"{print $1, $3}' >  size_chrom_query.txt

awk '{print $6":"$7}' ${PAF_FILE_TMP} | sort | uniq | tr ":" "\t" | \
    awk 'OFS="\t"{pritn $1, 0, $2}' | bedtools sort | awk 'OFS="\t"{print $1, $3}' >  size_chrom_reference.txt

#REF -> QUERY
if [[ "${CHROM_COMMUN_FILE}" != "NONE" ]] && test -s ${CHROM_COMMUN_FILE} ; then
    awk 'NR == FNR {  split($0, sp, ":"); chrom_commun[sp[1]]=sp[2]; chrom_commun[sp[2]]=sp[1]; next; } chrom_commun[$1] == $6 && OFS="\t"{print $6, $8, $9, $1, $3, $4}' \
        ${CHROM_COMMUN_FILE} ${PAF_FILE_TMP} | bedtools sort > refSim_sorted_tmp.bed
else
    awk 'substr($1, 1, 11) == substr($6, 1, 11) && OFS="\t"{print $6, $8, $9, $1, $3, $4}' \
        ${PAF_FILE_TMP}  | bedtools sort > refSim_sorted_tmp.bed
fi;

#--------------
grep -E ${CHROM_KEPPING_REGEX} size_chrom_query.txt > size_chrom_keep.txt


rm -f ${OUTPUT_GENOME}*
rm -f ${GENOME_FILE}.fai
rm -f *.fa*.fai
for chrom_size in `cat size_chrom_keep.txt | tr "\t" ":"`
do
    
    chrom=`echo $chrom_size | cut -d ":" -f 1`
    size=`echo $chrom_size | cut -d ":" -f 2`
    
    echo "chrom=$chrom ; size=$size"
    #chrom="2L_RaGOO_RaGOO"
    #deplace les coord query vers la gauche et ref vers la droite
    grep -w ${chrom} refSim_sorted_tmp.bed | awk 'OFS="\t"{print $4, $5, $6, $1, $2, $3}' > ${chrom}_sim_norm.bed 

    #ajout des taille est tranformation en cluster
    cat ${chrom}_sim_norm.bed | bedtools sort | bedtools cluster | awk 'OFS="\t"{print $1, $2, $3, $3-$2, $4, $5, $6, $6-$5, $7}' > clust.bed

    #trier sur start ref et ajoute un numéro d'ordre selon le start ref et le clust
    cat clust.bed | sort -k 6 -n | awk 'BEGIN{num=1; clus=0}OFS="\t"{
            if(clus == 0){
                clus=$9;
            }
            else if(clus != $9){
                num+=1;
                clus=$9;
            }

            print $0, num, NR
        }' > all_goups.bed

    cp all_goups.bed ${chrom}_all_goups.bed

    rm -f new_test.bed
    #9 = cluster
    for clust in `cat all_goups.bed | cut -f 9 | sort -n | uniq`
    do
        #echo "clust = $clust"
        cat all_goups.bed | awk -v clust="$clust" '$9==clust' > cl_or.bed
        #nombre de groupe différent
        nb_group=`cat all_goups.bed | awk -v clust="$clust" '$9==clust' | cut -f 10 | sort | uniq | wc -l`
        #nombre de ligne possédant le cluster
        nb_clust=`cat all_goups.bed | awk -v clust="$clust" '$9==clust' | cut -f 9 | wc -l`
        #echo "nb_group = $nb_group"
        rm -f merge.bed && touch merge.bed
        if [[ $nb_clust -ne 1 && $nb_group -eq 1 ]]; then
            #statements
            #echo "here"

            awk -v clust="$clust" '$9==clust' all_goups.bed > t1.bed
            awk -v clust="$clust" '$9==clust' all_goups.bed | bedtools sort > t2.bed

            #merge
            size1=`awk -v clust="$clust" '$9==clust' all_goups.bed | bedtools sort | \
                bedtools merge | awk 'BEGIN{size=0}{size+=$3-$2}END{print size}'`

            #no merge
            size2=`awk -v clust="$clust" '$9==clust' all_goups.bed | bedtools sort | \
                awk 'BEGIN{size=0}{size+=$3-$2}END{print size}'`

            # DEBUG
            # if [[ $clust -eq 7 ]]; then
            #     awk -v clust="$clust" '$9==clust' all_goups.bed
            #     echo "$size1   -- $size2";
            # fi
            diff t1.bed t2.bed > log.txt
            #$clust -eq 7
            #&& $size1 -eq $size2
            # if [[ $? -ne 0 && $clust -eq 7 ]]; then
            #     echo "clust = $clust"
            #     awk -v clust="$clust" '$9==clust' all_goups.bed | bedtools sort \
            #     | awk 'BEGIN{endp=0}OFS="\t"{if(endp==0){print $0;endp=$3;}else{if(endp<$3){$2=endp; print $0; endp=$3}}}' \
            #     | sort -k 11 -n | cut -f 1-3,5-7,9-11 >> new_test.bed
            # else
                for group in `cat cl_or.bed | cut -f 10 | sort | uniq`
                do
                    #chrom ref start ref end ref
                    chrom_tmp=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 5`
                    start=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 6`
                    end=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 7`
                    #chrom ref start ref end ref ----

                    awk -v group="$group" '$10==group' cl_or.bed | bedtools sort | bedtools merge > ttmp_merge.bed
                    bedtools subtract -a ttmp_merge.bed -b merge.bed > tmp_merge.bed
                    cat tmp_merge.bed >> merge.bed
                    nb_line=`cat tmp_merge.bed | wc -l`

                    order=`awk -v group="$group" -v clust="$clust" '$10==group && $9==clust {print $11}' cl_or.bed | head -n 1`

                    for (( line = 1; line <= $nb_line; line++ )); do
                        awk -v line="$line" -v order="$order" -v clust="$clust" -v group="$group" -v chrom_tmp="$chrom_tmp" -v start="$start" -v end="$end" 'NR==line && OFS="\t" {print $0, chrom_tmp, start, end, clust, group, order}' tmp_merge.bed >> new_test.bed
                    done;
                done;
            # fi
            
            #cat cl_or.bed | cut -f 1-3,5-7,9-11 >> new_test.bed
        elif [[ $nb_group -ne 1 ]] || [[ $nb_clust -eq 1 ]]; then
            #statements
            
            for group in `cat cl_or.bed | cut -f 10 | sort | uniq`
            do
                #chrom ref start ref end ref
                chrom_tmp=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 5`
                start=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 6`
                end=`awk -v group="$group" '$10==group' cl_or.bed | head -n 1 | cut -f 7`
                #chrom ref start ref end ref ----

                awk -v group="$group" '$10==group' cl_or.bed | bedtools sort | bedtools merge > ttmp_merge.bed
                bedtools subtract -a ttmp_merge.bed -b merge.bed > tmp_merge.bed
                cat tmp_merge.bed >> merge.bed
                nb_line=`cat tmp_merge.bed | wc -l`

                order=`awk -v group="$group" -v clust="$clust" '$10==group && $9==clust {print $11}' cl_or.bed | head -n 1`

                for (( line = 1; line <= $nb_line; line++ )); do
                    awk -v line="$line" -v order="$order" -v clust="$clust" -v group="$group" -v chrom_tmp="$chrom_tmp" -v start="$start" -v end="$end" 'NR==line && OFS="\t" {print $0, chrom_tmp, start, end, clust, group, order}' tmp_merge.bed >> new_test.bed
                done;
            done;
        fi
    done;


    cat new_test.bed | cut -f 1-9 > tmp_new.bed
    cat tmp_new.bed | sort -k 9 -n > new_test.bed

    awk 'OFS="\t"{print $1, $2, $3, $4, $5, $6, $2-$5}' new_test.bed| awk -v border_diag_size="$border_diag_size" 'OFS="\t"{$7=($7>0)?$7:-1*$7; if($7>border_diag_size){print $0, "I"}else{print $0, "N"}}' > ${chrom}_IN.bed 
    cat ${chrom}_IN.bed | awk 'BEGIN{startn=0}OFS="\t"{if($8=="N"){startn=$2; print $0, $1, $2, $3;}else{startn+=1;print $0, $1, startn, startn+1}}' > ${chrom}_IN_pos.bed
    awk 'OFS="\t"{print $1, $2, $3, $9, $10, $11}' ${chrom}_IN_pos.bed > ${chrom}_IN_pos2.bed
    cat ${chrom}_IN_pos2.bed | bedtools sort | bedtools merge > ${chrom}_IN_merge.bed
    bedtools complement -i ${chrom}_IN_merge.bed -g size_chrom_query.txt | grep ${chrom} | bedtools sort | awk 'OFS="\t"{print $0, $0}' > ${chrom}_IN_complement.bed
    cat ${chrom}_IN_complement.bed ${chrom}_IN_pos2.bed > ${chrom}_IN_new_pos.bed

    for i in `cat ${chrom}_IN_new_pos.bed | tr "\t" ":"`; do 
        
        tmp_chrom=`echo $i | cut -d ":" -f 1`
        start=`echo $i | cut -d ":" -f 2`;
        end=`echo $i | cut -d ":" -f 3`;
        #echo "$tmp_chrom ; start ; end"

        nb=`grep -w "$tmp_chrom" ${PAF_FILE_TMP} | grep -w "$end" | grep -w "$start" | cut -f 5 | head -n 1 | wc -l`
        if [[ $nb -eq 1 ]]; then
            strand=`grep -w "$tmp_chrom" ${PAF_FILE_TMP} | grep -w "$end" | grep -w "$start" | cut -f 5 | head -n 1`
        else
            strand="+"
        fi 
        echo -e $tmp_chrom"\t"$start"\t"$end"\tname""\t1""\t$strand"

    done > chrom_${chrom}.bed;

    #get morceau sequence
    bedtools getfasta -fi ${GENOME_FILE} -bed chrom_${chrom}.bed -s > tmp_chrom_${chrom}.fasta
    awk -v chrom="$chrom" 'BEGIN{print ">"chrom; chain=""} substr($0, 1, 1) != ">" {chain = chain""$0}END{print chain}' tmp_chrom_${chrom}.fasta > chrom_${chrom}.fasta

    size_chrom_fa=`grep -v ">" chrom_${chrom}.fasta | awk '{print length($0)}'`
    echo "$chrom : size_chrom_fa === $size_chrom_fa"

    cat chrom_${chrom}.fasta >> ${OUTPUT_GENOME}
done;

wc -c ${OUTPUT_GENOME} ${GENOME_FILE}




#bash ../correction_genome_transloc_sim_chrom.sh 3R_G0-F100.Fasta 3R_G31.fasta core.fasta
#---------------------------

# rm -f corrected_genome_sim_*.fasta* 
# e=0;
# GENOME_FILE="corrected_genome_3.fasta"
# cat ${GENOME_FILE} > corrected_genome_sim_${e}.fasta
# for (( i = 1; i < 7; i++ )); do
#     bash ../correction_genome_transloc_sim_chrom.sh G0-F100_chrom.fasta corrected_genome_sim_${e}.fasta corrected_genome_sim_${i}.fasta
#     echo "$e"  
#     e=$(($e+1));
# done;


# rm -f corrected_genome_sim_*.fasta* 
# e=0;
# GENOME_FILE="corrected_genome_3.fasta"
# cat ${GENOME_FILE} > corrected_genome_sim_${e}.fasta
# for (( i = 1; i < 7; i++ )); do
#     bash ../correction_genome_transloc_sim_chrom.sh -c chrom.txt -r "^[23][LR]_|^[4]_|^[X]_" G0-F100_chrom.fasta corrected_genome_sim_${e}.fasta corrected_genome_sim_${i}.fasta
#     echo "$e"  
#     e=$(($e+1));
# done;


