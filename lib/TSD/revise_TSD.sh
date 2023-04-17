GENOME=$1
#FILE_TSD="total_results_tsd_ZAM_KO.txt"
FILE_TSD=$2
SIZE_FLANK=$3
SIZE_KMER=$4
WORK_TMP_DIR=$5

echo "[$0] <<<< PARAMS >>>> ";
echo "[$0] GENOME : $GENOME ";
echo "[$0] FILE_TSD : $FILE_TSD ";
echo "[$0] SIZE_FLANK : $SIZE_FLANK ";
echo "[$0] SIZE_KMER : $SIZE_KMER ";
echo "[$0] WORK_TMP_DIR : $WORK_TMP_DIR ";
echo "";


NAME_FILE=`basename $FILE_TSD | sed 's/.txt//g'`

echo "[$0] NAME_FILE : $NAME_FILE";

grep -A 2 -B 2 "KO:" $FILE_TSD > ${NAME_FILE}_KO.txt 
grep -B 1 KO ${NAME_FILE}_KO.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.bed

bedtools getfasta -fi $GENOME -bed ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.bed -name+ > ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.fasta

rm -f ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.bed

echo ${NAME_FILE}_KO_corrected.txt;
python3 `dirname $0`/revise_TSD.py ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.fasta ${NAME_FILE}_KO.txt $SIZE_KMER > ${NAME_FILE}_KO_corrected.txt

rm -f ${NAME_FILE}_KO.txt


#GET TSD
grep "OK:" $FILE_TSD | cut -d"," -f 2 | tr -d " " > ${NAME_FILE}_TSD_OK.txt
grep -o "TSM=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> ${NAME_FILE}_TSD_OK.txt

rm -f ${NAME_FILE}_TSD_OK.txt

# grep "OK:" $FILE_TSD | cut -d"," -f 2 | tr -d " " > TSD_OK.txt
# grep -o "TSM=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> TSD_OK.txt

#GET TSM
grep -A 1 OK $FILE_TSD | grep ">" | tr -d ">" | awk -F":" '{print substr($1, length($1)-5, 2)""substr($2, 1, 6)}' > ${NAME_FILE}_TSM_OK.txt
grep -o "SVI=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> ${NAME_FILE}_TSM_OK.txt

rm -f ${NAME_FILE}_TSM_OK.txt

grep -A 1 OK $FILE_TSD | grep ">" | tr -d ">" | awk -F":" '{print substr($1, length($1)-5, 2)""substr($2, 1, 6)}' > ${WORK_TMP_DIR}/TSM_OK_${NAME_FILE}.txt
grep -o "SVI=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> ${WORK_TMP_DIR}/TSM_OK_${NAME_FILE}.txt

rm -f ${NAME_FILE}_KO_corrected.txt


#OK
grep -A 2 -B 2 OK $FILE_TSD > ${NAME_FILE}_OK.txt 
grep -B 1 OK ${NAME_FILE}_OK.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.bed

rm -f ${NAME_FILE}_OK.txt


bedtools getfasta -fi $GENOME -bed ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.bed -name+ > ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.fasta

cat ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.fasta ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.fasta > ${WORK_TMP_DIR}/all_empty_site_${NAME_FILE}.fasta

rm -f ${WORK_TMP_DIR}/empty_site_KO_${NAME_FILE}.fasta
rm -f ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.fasta
rm -f ${WORK_TMP_DIR}/empty_site_OK_${NAME_FILE}.bed
rm -f ${WORK_TMP_DIR}/TSM_OK_${NAME_FILE}.txt
rm -f ${WORK_TMP_DIR}/all_empty_site_${NAME_FILE}.fasta

# echo
# echo ${NAME_FILE}_TSM.txt;
# echo
# python get_TSM.py all_empty_site_${NAME_FILE}.fasta $SIZE_KMER > ${NAME_FILE}_TSM.txt


