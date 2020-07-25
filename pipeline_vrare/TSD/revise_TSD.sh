

GENOME=$1
#FILE_TSD="total_results_tsd_ZAM_KO.txt"
FILE_TSD=$2
SIZE_FLANK=$3
SIZE_KMER=$4

NAME_FILE=`echo $FILE_TSD | cut -d"." -f 1`

grep -A 2 -B 2 KO $FILE_TSD > ${NAME_FILE}_KO.txt 
grep -B 1 KO ${NAME_FILE}_KO.txt | grep ">" | tr -d ">" | awk -F":" -v s_fk="$SIZE_FLANK" 'OFS="\t" {print $1"_RaGOO_RaGOO", $3-s_fk, $3, $5":FLANK_LEFT\n"$1"_RaGOO_RaGOO", $3, $3+s_fk, $5":FLANK_RIGHT"}' > empty_site.bed

bedtools getfasta -fi $GENOME -bed empty_site.bed -name > empty_site.fasta

python revise_TSD.py empty_site.fasta ${NAME_FILE}_KO.txt $SIZE_KMER > ${NAME_FILE}_KO_corrected.txt

