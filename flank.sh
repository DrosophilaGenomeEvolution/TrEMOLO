reads=$1
size=$2
find_file=$3


id=`echo $reads | grep -o ":[0-9]*:" | grep -o "[0-9]*"`
name=`echo $reads | grep -o -e "_.*\." | grep -o "[^_].*[^.]"`

head=`grep ".*:.*:$id:[0-9]*:[PI]" $find_file`

awk -v var=$head 'BEGIN {nb=0} { if( var == $0 || nb == 1 ){print $0; nb = nb + 1;} }' $find_file > seq_test.fasta
if [[ `cat seq_test.fasta` = "" ]]; then
	exit 1
fi

makeblastdb -in ./$reads -dbtype nucl
blastn -db ./$reads -query seq_test.fasta \
	-perc_identity 100 \
    -outfmt 6 \
    -out seq_test.bln;


awk -v var=$size '{if($9 - var > 0){ if($9 < $10){ print $2"\t"$9-var-1"\t"$9-1"\n" $2"\t"$10"\t"$10+var }else{ print $2"\t"$10-var-1"\t"$10-1"\n" $2"\t"$9"\t"$9+var  } } }' seq_test.bln > flank.bed
awk '{ if($9 < $10){ print $2"\t"$9-1"\t"$10"\t" "forward" "\t" "1" "\t" "+" }else{ print $2"\t"$10-1"\t"$9"\t" "reverse" "\t" "1" "\t" "-" } }' seq_test.bln > seq_test.bed


bedtools getfasta -fi $reads -bed flank.bed > FLANK/Flank_${name}.fasta
