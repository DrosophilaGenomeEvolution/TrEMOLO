#ARGS
FIND_FA=$1
REPO_READS=$2
REPO_READS_FA=$3
DB_TE=$4
FLANK_SIZE=$5
TSD_SIZE=$6

#ERROR
if [ "$#" -ne 6 ]; then
 echo "ERROR : need 6 arguments. You have put $# arguments" ;
 exit 1 ;
fi;


if [ ! -d "$REPO_READS" ]; then
	"$REPO_READS is not a directory"
	exit 1 ;
fi;


if [ ! -d "$REPO_READS_FA" ]; then
 	"$REPO_READS_FA is not a directory"
 	exit 1 ;
fi;

REPO_READS=`echo ${REPO_READS} | sed 's/[/]$//g'`

REPO_READS_FA=`echo ${REPO_READS_FA} | sed 's/[/]$//g'`

#for e in `ls all_fasta_element/ | grep -o "d_.*\." | grep -o "_.*[^.]" | grep -o "[^_]*"`; do
#elem=$e
echo "<<<<<<<<<<<<<<<<<<< BEGIN >>>>>>>>>>>>>>>>>>>>"
echo " " > total_results_tsd.txt
nombre_element=`grep ">" $FIND_FA | grep -o "[0-9]:[0-9]*:[0-9]*:[PI]" | grep -o ":[0-9]*:" | grep -o "[0-9]*" | wc -l`
i=0
for id in `grep ">" ${FIND_FA} | grep -o "[0-9]:[0-9]*:[0-9]*:[PI]" | grep -o ":[0-9]*:" | grep -o "[0-9]*"`; do
		echo "------FIND FILE READS------"
		fr="`ls ${REPO_READS} | grep ":$id:"`"
		echo "id : $id"
		echo "file : $fr"
        name=`echo $fr | grep -o ".*\."`
        echo "name : $name"
        i=$(($i + 1))
        echo "$i/$nombre_element"

        reads=${REPO_READS_FA}/${name}fasta

		echo 
		echo "------TSD------"
		id=`echo $reads | grep -o ":[0-9]*:" | grep -o "[0-9]*"`
		echo "reads : "$reads
		echo "id : "$id
		echo "find_file : ${FIND_FA}"
		head=`grep ".*:.*:$id:[0-9]*:[PI]" ${FIND_FA}`

		echo "head : "$head

		awk -v var=$head 'BEGIN {nb=0} { if( var == $0 || nb == 1 ){print $0; nb = nb + 1;} }' "${FIND_FA}" > sequence_TE.fasta
		if [[ `cat sequence_TE.fasta` = "" ]]; then
			echo "ERROR : can't get sequence TE in ${FIND_FA}" ;
			exit 1 ;
		fi

		#
		echo "*********BLAST 1 TE VS READ**********"
		makeblastdb -in "$reads" -dbtype nucl
		blastn -db "$reads" -query sequence_TE.fasta \
			-perc_identity 100 \
		    -outfmt 6 \
		    -out sequence_TE.bln;

		#get flank bed file
		awk -v var=$FLANK_SIZE '{if($9 - var > 0){ if($9 < $10){ print $2"\t"$9-var-1"\t"$9-1"\n" $2"\t"$10"\t"$10+var }else{ print $2"\t"$10-var-1"\t"$10-1"\n" $2"\t"$9"\t"$9+var  } } }' sequence_TE.bln > flank_TE.bed
		#get TE SEQ
		awk '{ if($9 < $10){ print $2"\t"$9-1"\t"$10"\t" "forward" "\t" "1" "\t" "+" }else{ print $2"\t"$10-1"\t"$9"\t" "reverse" "\t" "1" "\t" "-" } }' sequence_TE.bln > sequence_TE.bed


		bedtools getfasta -fi $reads -bed flank_TE.bed > flank_TE.fasta
		bedtools getfasta -fi $reads -bed sequence_TE.bed -name > sequence_TE.fasta

		echo "*********BLAST 2 TE VS DBTE**********"
		echo ${DB_TE}
		makeblastdb -in ${DB_TE} -dbtype nucl
		blastn -db ${DB_TE} -query sequence_TE.fasta -outfmt 6 -out TE_vs_databaseTE.bln

		head -n 1 TE_vs_databaseTE.bln
		strand=`awk 'NR==1 {if($9 < $10) {print "+"}else{print "-"}}' TE_vs_databaseTE.bln`
		echo $reads >> total_results_tsd.txt
		echo $head  >> total_results_tsd.txt

		#FIND TSD
		python find_tsd.py flank_TE.fasta sequence_TE.fasta $FLANK_SIZE $id $strand $TSD_SIZE >> total_results_tsd.txt
done


rm -f TE_vs_databaseTE.bln
rm -f sequence_TE.fasta sequence_TE.bln sequence_TE.bed
rm -f flank_TE.fasta flank_TE.bed


nombre_ok=`grep OK total_results_tsd.txt | wc -l`
nombre_ko=`grep KO total_results_tsd.txt | wc -l`
nombre_total=`grep reads total_results_tsd.txt | wc -l`
nombre_k_o=`grep "K-O" total_results_tsd.txt | wc -l`

#RESUME
echo "OK/total : $nombre_ok/$nombre_total" >> total_results_tsd.txt
echo "KO/total : $nombre_ko/$nombre_total" >> total_results_tsd.txt
echo "OK+KO/total : $(($nombre_ok+$nombre_ko))/$nombre_total" >> total_results_tsd.txt
echo "K-O/total : $nombre_k_o/$nombre_total" >> total_results_tsd.txt
echo "OK+K-O/total : $(($nombre_ok+$nombre_k_o))/$nombre_total" >> total_results_tsd.txt
echo "OK% : $(($nombre_ok*100/$nombre_total))%" >> total_results_tsd.txt




