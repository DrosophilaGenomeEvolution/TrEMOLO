
OUTPUT=$1

# FILE_INS=$1
# GENOME=$2

# awk 'NR>1 && OFS="\t"{split($2, a, ":"); split(a[6], b, "-"); print a[5], b[1], b[2], $1"|"a[1]}' $FILE_INS > INS.bed
# awk 'NR>1 {print $1}' $FILE_INS | sort -u > TE_name.txt

# for TE in `cat TE_name.txt`; do
#     echo "[$0] TE = $TE"
#     grep "^$TE\s" $FILE_INS | awk 'OFS="\t"{split($2, a, ":"); split(a[6], b, "-"); print a[5], b[1], b[2], $1"|"a[1]}' > INS_$TE.bed
#     rm -f ${OUTPUT}
#     for i in `cat INS_$TE.bed | tr "\t" ":"`; 
#     do 


#         echo $i | tr ":" "\t" > TE.bed;
#         bedtools getfasta -fi $GENOME -bed TE.bed > TE.fasta


#         ID=`cat TE.bed | cut -f 4 | cut -d "|" -f 2`
#         #echo "[$0] ID = $ID"


#         strand=`grep -w $ID $FILE_INS | awk 'OFS="\t"{split($2, a, ":"); print a[7]}'`
#         echo "[$0] ID = $ID; strand = $strand;"


#         cat TE.bed | awk -v size_flank="30" \
#         'OFS="\t"{print $1, $2-size_flank, $2; print $1, $3, $3+size_flank;}' > FLANK.bed;
#         bedtools getfasta -fi $GENOME -bed FLANK.bed > FLANK.fasta

#         python3 `dirname $0`/find_tsd.py FLANK.fasta TE.fasta 10 $ID $strand -1 >> ${OUTPUT} 
#     done;

number_ok=`grep OK ${OUTPUT} -c`
number_ko=`grep KO ${OUTPUT} -c`
number_total=`grep ">" ${OUTPUT} -c`
number_k_o=`grep "K-O" ${OUTPUT} -c`

#RESUME
echo "OK/total : $number_ok/$number_total" >> ${OUTPUT}
echo "KO/total : $number_ko/$number_total" >> ${OUTPUT}
echo "OK+KO/total : $(($number_ok+$number_ko))/$number_total" >> ${OUTPUT}
echo "K-O/total : $number_k_o/$number_total" >> ${OUTPUT}
echo "OK+K-O/total : $(($number_ok+$number_k_o))/$number_total" >> ${OUTPUT}
echo "OK% : $(($number_ok*100/$number_total))%" >> ${OUTPUT}
    #mv ${OUTPUT} TSD_$TE.txt
# done;
# rm -f FLANK.bed TE.bed TE_name.txt

