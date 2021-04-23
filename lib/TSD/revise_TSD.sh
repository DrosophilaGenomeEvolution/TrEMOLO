GENOME=$1
#FILE_TSD="total_results_tsd_ZAM_KO.txt"
FILE_TSD=$2
SIZE_FLANK=$3
SIZE_KMER=$4

NAME_FILE=`echo $FILE_TSD | sed 's/.txt//g'`

grep -A 2 -B 2 "KO:" $FILE_TSD > ${NAME_FILE}_KO.txt 
grep -B 1 KO ${NAME_FILE}_KO.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > empty_site_KO.bed

bedtools getfasta -fi $GENOME -bed empty_site_KO.bed -name+ > empty_site_KO.fasta

echo ${NAME_FILE}_KO_corrected.txt;
python3 `dirname $0`/revise_TSD.py empty_site_KO.fasta ${NAME_FILE}_KO.txt $SIZE_KMER > ${NAME_FILE}_KO_corrected.txt

#GET TSD
grep "OK:" $FILE_TSD | cut -d"," -f 2 | tr -d " " > ${NAME_FILE}_TSD_OK.txt
grep -o "TSM=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> ${NAME_FILE}_TSD_OK.txt

# grep "OK:" $FILE_TSD | cut -d"," -f 2 | tr -d " " > TSD_OK.txt
# grep -o "TSM=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> TSD_OK.txt

#GET TSM
grep -A 1 OK $FILE_TSD | grep ">" | tr -d ">" | awk -F":" '{print substr($1, length($1)-5, 2)""substr($2, 1, 6)}' > ${NAME_FILE}_TSM_OK.txt
grep -o "SVI=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> ${NAME_FILE}_TSM_OK.txt

grep -A 1 OK $FILE_TSD | grep ">" | tr -d ">" | awk -F":" '{print substr($1, length($1)-5, 2)""substr($2, 1, 6)}' > TSM_OK.txt
grep -o "SVI=[A-Z]*" ${NAME_FILE}_KO_corrected.txt | cut -d"=" -f 2 >> TSM_OK.txt

#OK
grep -A 2 -B 2 OK $FILE_TSD > ${NAME_FILE}_OK.txt 
grep -B 1 OK ${NAME_FILE}_OK.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > empty_site_OK.bed

bedtools getfasta -fi $GENOME -bed empty_site_OK.bed -name+ > empty_site_OK.fasta

cat empty_site_OK.fasta empty_site_KO.fasta > all_empty_site.fasta

rm -f all_empty_site* empty_site*

# echo
# echo ${NAME_FILE}_TSM.txt;
# echo
# python get_TSM.py all_empty_site.fasta $SIZE_KMER > ${NAME_FILE}_TSM.txt


