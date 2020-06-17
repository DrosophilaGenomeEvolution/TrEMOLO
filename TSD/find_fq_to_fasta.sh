#Convert fastq file of support READS TE to fasta file 
#
#


#ARGS
#File find TE fasta 
FIND_FA=$1
#File contains READS support for find TE
REPO_READS=$2
#Name out directory
OUT_DIR=$3


if [ "$#" -ne 3 ]; then
 echo "ERROR : need 3 arguments. You have put $# arguments" ;
 exit 1 ;
fi

REPO_READS=`echo ${REPO_READS} | sed 's/[/]$//g'`

mkdir -p ${OUT_DIR}
OUT_DIR=`echo ${OUT_DIR} | sed 's/[/]$//g'`


echo " " > total_results_tsd.txt
nombre_element=`grep ">" ${FIND_FA} | grep -o "[0-9]:[0-9]*:[0-9]*:[PI]" | grep -o ":[0-9]*:" | grep -o "[0-9]*" | wc -l`
i=0
echo "begin $nombre_element"
for id in `grep ">" ${FIND_FA} | grep -o "[0-9]:[0-9]*:[0-9]*:[PI]" | grep -o ":[0-9]*:" | grep -o "[0-9]*"`; do
		fr="`ls ${REPO_READS} | grep ":$id:"`"
		echo "id : $id"
		echo "file : $fr"
        name=`echo $fr | grep -o ".*\."`
        echo "name : $name"
		python3 fastq_to_fasta.py ${REPO_READS}/$fr ./${OUT_DIR}/${name}fasta
        i=$(($i + 1))
        echo "$i/$nombre_element"
done

