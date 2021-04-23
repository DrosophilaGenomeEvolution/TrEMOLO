###################################################################################################################################
#
# Copyright 2019-2020 IRD-CNRS-Lyon1 University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
# If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Mourdas Mohamed
#                                                                                                                                   
####################################################################################################################################

# FIND TSD
# ==========================
# :author:  Mourdas MOHAMED 
# :contact: mourdas.mohamed@igh.cnrs.fr
# :date: 01/06/2020
# :version: 0.1
# Script description
# ------------------
# tsd_te.sh Convert fastq file of support READS TE to fasta file
# -------
# # tsd_te.sh prefix_find_ZAM.fasta READS_SUPPORT/ READS_FASTA/ cannonical_TE.fa 30 4
# 
# Help Programm
# -------------
# usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size>
# ***WARNING: this program needs the find_tsd.py script in the same directory***


#ARGS
FIND_FA=$1
REPO_READS=$2
REPO_READS_FA=$3 #FASTA FOR BLASTN
DB_TE=$4
FLANK_SIZE=$5
TSD_SIZE=$6
COMBINE=$7

path_this_script=`dirname $0`
echo $path_this_script ;


#ERROR
if [ "$#" -ne 7 ]; then
    echo "ERROR : need 6 arguments. You have put $# arguments" ;
    echo "usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size>"
    exit 1 ;
fi;


if [ ! -d "$REPO_READS" ]; then
    "$REPO_READS is not a directory"
    echo "usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size>"
    exit 1 ;
fi;


if [ ! -d "$REPO_READS_FA" ]; then
    "$REPO_READS_FA is not a directory"
    echo "usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size>"
    exit 1 ;
fi;

REPO_READS=`echo ${REPO_READS} | sed 's/[/]$//g'`

REPO_READS_FA=`echo ${REPO_READS_FA} | sed 's/[/]$//g'`

#for e in `ls all_fasta_element/ | grep -o "d_.*\." | grep -o "_.*[^.]" | grep -o "[^_]*"`; do
#elem=$e
echo "[$0] <<<<<<<<<<<<<<<<<<< TSD >>>>>>>>>>>>>>>>>>>>"
echo " " > total_results_tsd.txt
number_element=`grep ">" $FIND_FA | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*" | wc -l`
i=0
for id in `grep ">" ${FIND_FA} | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*"`; do
        
        if [ ! -n "$id" ]; then
            echo "[$0] **ERROR** ID NOT FOUND : $id";
            exit 1;
        fi

        echo "[$0] ------FIND INFOS READS------"
        fr="`ls ${REPO_READS} | grep ":$id:"`"
        echo "[$0] id : $id--"
        echo "[$0] file : $fr--"
        name=`echo $fr | grep -o ".*\."`
        echo "[$0] name : $name"
        i=$(($i + 1))
        echo "[$0] $i/$number_element"

        if [ ! -n "$fr" ]; then
            echo "[$0] **ERROR** : FILE NOT FOUND FOR ID : $id";
            exit 1;
        fi

        reads=${REPO_READS_FA}/${name}fasta

        echo 
        echo "[$0] ------TSD------"
        id=`echo $reads | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*"`
        echo "reads : "$reads
        echo "id : "$id
        echo "find_file : ${FIND_FA}"
        head=`grep ".*:.*:$id:[0-9]*:[PI]" ${FIND_FA}`

        echo "[$0] head : "$head

        #awk -v var=$head 'BEGIN {nb=0} { if( var == $0 || nb == 1 ){print $0; nb = nb + 1;} }' "${FIND_FA}" > sequence_TE.fasta
        grep -w $head -A 1 "${FIND_FA}" > sequence_TE.fasta 
        # cat sequence_TE.fasta 
        
        if ! test -s sequence_TE.fasta; then
            echo "**ERROR** : can't get sequence TE in ${FIND_FA}" ;
            exit 1 ;
        fi;

        #
        echo "*********BLAST 1 TE VS READ**********"
        makeblastdb -in "$reads" -dbtype nucl
        blastn -db "$reads" -query sequence_TE.fasta \
            -perc_identity 100 \
            -outfmt 6 \
            -out sequence_TE.bln ;


        # get flank bed file
        awk -v var=$FLANK_SIZE '{if($9 - var > 0){ if($9 < $10){ print $2"\t"$9-var-1"\t"$9-1"\n" $2"\t"$10"\t"$10+var }else{ print $2"\t"$10-var-1"\t"$10-1"\n" $2"\t"$9"\t"$9+var  } } }' sequence_TE.bln > flank_TE.bed
        # get TE SEQ
        awk '{ if($9 < $10){ print $2"\t"$9-1"\t"$10"\t" "forward" "\t" "1" "\t" "+" }else{ print $2"\t"$10-1"\t"$9"\t" "reverse" "\t" "1" "\t" "-" } }' sequence_TE.bln > sequence_TE.bed

        mkdir -p DIR_SEQ_TE_READ_POS
        
        cp sequence_TE.bed  DIR_SEQ_TE_READ_POS/sequence_TE_${name}.bed
        cp sequence_TE.bln  DIR_SEQ_TE_READ_POS/sequence_TE_${name}.bln


        bedtools getfasta -fi $reads -bed flank_TE.bed > flank_TE.fasta
        bedtools getfasta -fi $reads -bed sequence_TE.bed -name > sequence_TE.fasta

        echo "[$0] *********BLAST 2 TE VS DBTE**********"
        echo ${DB_TE}
        makeblastdb -in ${DB_TE} -dbtype nucl
        blastn -db ${DB_TE} -query sequence_TE.fasta -outfmt 6 -out TE_vs_databaseTE.bln


#COMBINE AWK
# sstart_global=`awk 'NR==1 {print $7}' test.bln`
# send_global=`awk 'NR==1 {print $8}' test.bln`
# qstart_global=`awk 'NR==1 {print $9}' test.bln`
# qend_global=`awk 'NR==1 {print $10}' test.bln`


# awk -v ssg="$sstart_global" -v seg="$send_global" -v qsg="$qstart_global" -v qeg="$qend_global" '
#     BEGIN{
#         chrom=""
#         sstart_global=ssg; 
#         send_global=seg; 
#         qstart_global=qsg; 
#         qend_global=qeg; 
#     }

#      OFS="\t"  {
#         chrom=$1
#         TE=$2
#         sstart = $7
#         ssend  = $8
#         qstart = $9
#         qend   = $10
#         if (sstart < sstart_global && send < send_global && qstart < qstart_global && qend < qend_global ){
#             sstart_global = sstart
#             qstart_global = qstart
#         }
#         else if ( sstart > sstart_global && send > send_global && qstart > qstart_global && qend > qend_global ){
#             send_global   = send
#             qend_global   = qend
#         }
#     }

#    END{split(chrom, a, ":"); print chrom, sstart_global, send_global, TE"|"a[5]; }' test.bln




        head -n 1 TE_vs_databaseTE.bln

        strand=`awk 'NR==1 {if ($9 < $10){print "+"} else {print "-"}}' TE_vs_databaseTE.bln`
        echo $reads >> total_results_tsd.txt
        echo $head  >> total_results_tsd.txt

        # FIND TSD
        # Warning : path
        echo "[$0] flank_TE.fasta    sequence_TE.fasta     $FLANK_SIZE     $id     $strand     $TSD_SIZE"
        python3 ${path_this_script}/find_tsd.py flank_TE.fasta sequence_TE.fasta $FLANK_SIZE $id $strand $TSD_SIZE >> total_results_tsd.txt
done


rm -f TE_vs_databaseTE.bln
rm -f sequence_TE.fasta sequence_TE.bln sequence_TE.bed
rm -f flank_TE.fasta flank_TE.bed


number_ok=`grep "OK" total_results_tsd.txt -c`
number_ko=`grep "KO" total_results_tsd.txt -c`
number_total=`grep "reads" total_results_tsd.txt -c`
number_k_o=`grep "K-O" total_results_tsd.txt  -c`


#RESUME
echo "OK/total : $number_ok/$number_total" >> total_results_tsd.txt
echo "KO/total : $number_ko/$number_total" >> total_results_tsd.txt
echo "OK+KO/total : $(($number_ok+$number_ko))/$number_total" >> total_results_tsd.txt
echo "K-O/total : $number_k_o/$number_total" >> total_results_tsd.txt
echo "OK+K-O/total : $(($number_ok+$number_k_o))/$number_total" >> total_results_tsd.txt
echo "OK% : $(($number_ok*100/$number_total))%" >> total_results_tsd.txt




