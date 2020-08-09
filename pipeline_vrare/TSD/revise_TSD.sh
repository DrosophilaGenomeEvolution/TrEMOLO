

GENOME=$1
#FILE_TSD="total_results_tsd_ZAM_KO.txt"
FILE_TSD=$2
SIZE_FLANK=$3
SIZE_KMER=$4

NAME_FILE=`echo $FILE_TSD | cut -d"." -f 1`

grep -A 2 -B 2 KO $FILE_TSD > ${NAME_FILE}_KO.txt 
grep -B 1 KO ${NAME_FILE}_KO.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > empty_site_KO.bed

bedtools getfasta -fi $GENOME -bed empty_site_KO.bed -name > empty_site_KO.fasta

echo ${NAME_FILE}_KO_corrected.txt;
python revise_TSD.py empty_site_KO.fasta ${NAME_FILE}_KO.txt $SIZE_KMER > ${NAME_FILE}_KO_corrected.txt


#OK
grep -A 2 -B 2 OK $FILE_TSD > ${NAME_FILE}_OK.txt 
grep -B 1 OK ${NAME_FILE}_OK.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1, $3-s_fk, $3, $5":FLANK_LEFT\n"$1, $3, $3+s_fk, $5":FLANK_RIGHT"}' > empty_site_OK.bed

bedtools getfasta -fi $GENOME -bed empty_site_OK.bed -name > empty_site_OK.fasta

echo all_empty_site.fasta;
cat empty_site_OK.fasta empty_site_KO.fasta > all_empty_site.fasta

echo ${NAME_FILE}_TSM.txt;
python get_TSM.py all_empty_site.fasta $SIZE_KMER > ${NAME_FILE}_TSM.txt


