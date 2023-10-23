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
OUTPUT=$8

#NEW
GENOME=$9

echo "[$0] <<<< PARAMS >>>> ";
echo "[$0] FIND_FA : $FIND_FA ";
echo "[$0] REPO_READS : $REPO_READS ";
echo "[$0] REPO_READS_FA : $REPO_READS_FA ";
echo "[$0] DB_TE : $DB_TE ";
echo "[$0] FLANK_SIZE : $FLANK_SIZE ";
echo "[$0] TSD_SIZE : $TSD_SIZE ";
echo "[$0] COMBINE : $COMBINE ";
echo "[$0] OUTPUT : $OUTPUT ";
echo "[$0] GENOME : $GENOME ";
echo "";

path_this_script=`dirname $0`
echo $path_this_script ;

path_out_dir=`dirname ${OUTPUT}`

#ERROR
if [ "$#" -ne 9 ]; then
    echo "[$0] ERROR : need 6 arguments. You have put $# arguments";
    echo "[$0] usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size> <genome>";
    exit 1 ;
fi;


if [ ! -d "$REPO_READS" ]; then
    echo "[$0] $REPO_READS is not a directory";
    echo "[$0] usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size> <genome>";
    exit 2 ;
fi;


if [ ! -d "$REPO_READS_FA" ]; then
    echo "[$0] $REPO_READS_FA is not a directory"
    echo "[$0] usage: tsd_te.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <DIRECTORY_READS_SUPPORT_FASTA> <database_TE_fasta> <flank_size> <tsd_size> <genome>"
    exit 3 ;
fi;

REPO_READS=`echo ${REPO_READS} | sed 's/[/]$//g'`

REPO_READS_FA=`echo ${REPO_READS_FA} | sed 's/[/]$//g'`


echo "[$0] <<<<<<<<<<<<<<<<<<< TSD >>>>>>>>>>>>>>>>>>>>"
echo " " > ${OUTPUT}
number_element=`grep ">" $FIND_FA | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*" | wc -l`
i=0
for id in `grep ">" ${FIND_FA} | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*"`; do
        
        if [ ! -n "$id" ]; then
            echo "[$0] **ERROR** ID NOT FOUND : $id";
            exit 1;
        fi

        echo "[$0] ------FIND INFOS READS------"
        fr="`ls ${REPO_READS} | grep ":$id:"`"
        type=`echo $id | grep -w -o -E "INS|DEL"`
        echo "[$0] id : $id--"
        echo "[$0] file : $fr--"
        echo "[$0] type : $type--"
        name=`echo $fr | grep -o ".*\."`
        echo "[$0] name : $name"
        i=$(($i + 1))
        echo "[$0] $i/$number_element"

        if [ ! -n "$fr" ]; then
            if [ "$type" = "DEL" ]; then
                echo 
                echo "[$0] ------TSD DEL------"
                echo "[$0] reads : NONE"
                echo "[$0] id : "$id
                echo "[$0] find_file : ${FIND_FA}"
                head=`grep ".*:.*:$id:[0-9]*:[PI]" ${FIND_FA}`

                echo "[$0] head : "$head


                grep -w $head -A 1 "${FIND_FA}" > ${path_out_dir}/sequence_TE_${id}.fasta 
                # cat ${path_out_dir}/sequence_TE_${id}.fasta 
                
                if ! test -s ${path_out_dir}/sequence_TE_${id}.fasta; then
                    echo "[$0] **ERROR** : can't get sequence TE in ${FIND_FA}" ;
                    exit 1 ;
                fi;

                grep -w $id ${path_out_dir}/COUNT_READS.txt | awk -v sizeFlank="30" '
                    $19!="." && max<$8 {max=$8; seq=$19; rd=$2; rf_cr=$9; count_read=$15;} 
                
                    END{
                        split(seq, sp, ":"); 
                        print ">"rd":"(rf_cr-sizeFlank)"-"(rf_cr)
                        print sp[1]
                        print ">"rd":"(rf_cr)"-"(rf_cr+sizeFlank)
                        print sp[2]
                    }' > ${path_out_dir}/flank_TE_${id}.fasta

                

                echo "[$0] head ${path_out_dir}/TE_vs_databaseTE_${id}.bln" >&2; 
                head -n 1 ${path_out_dir}/TE_vs_databaseTE_${id}.bln >&2
                echo "" >&2
                
                #Waring strand
                strand="+"
                echo $reads >> ${OUTPUT}
                echo $head  >> ${OUTPUT}

                # FIND TSD
                # Warning : path
                echo "[$0] flank_TE_${id}.fasta:"`test -s "${path_out_dir}/flank_TE_${id}.fasta" && echo "OK" || echo "NONE"`"   sequence_TE.fasta:"`test -s "${path_out_dir}/sequence_TE_${id}.fasta" && echo "OK" || echo "NONE"`"  FK_S:$FLANK_SIZE     id:$id     strand:$strand     TSD_SIZE:$TSD_SIZE"
                echo "";
                echo "";
                echo "[$0] WHY params ${path_out_dir}/flank_TE_${id}.fasta ${path_out_dir}/sequence_TE_${id}.fasta $FLANK_SIZE $id $strand $TSD_SIZE" >&2
                
                if test -s "${path_out_dir}/TE_vs_databaseTE_${id}.bln" && test -n "$FLANK_SIZE" && test -n "$id" && test -n "$strand" && test -n "$TSD_SIZE"; then
                    python3 ${path_this_script}/find_tsd.py ${path_out_dir}/flank_TE_${id}.fasta ${path_out_dir}/sequence_TE_${id}.fasta $FLANK_SIZE $id $strand $TSD_SIZE >> ${OUTPUT}
                    OK=`grep "$head" ${OUTPUT} -A 2 | grep -o "OK" || echo ""`
                else
                    OK="";
                fi;
                
                if [ -n "$OK" ]; then
                    TSD=`grep "$head" ${OUTPUT} -A 4 | grep "++:[A-Z]*:++" -o | head -n 1 | grep -o "[A-Z]*"`
                    echo "$head" | awk -v FLANK_SIZE="$FLANK_SIZE"  -F ":" 'OFS="\t"{print $1, $3-FLANK_SIZE, $3+FLANK_SIZE}' | tr -d ">" >> ${path_out_dir}/tmp_${id}.bed
                    empty_site=`bedtools getfasta -fi ${GENOME} -bed ${path_out_dir}/tmp_${id}.bed | grep -v ">"`
                    position_TE_SV=`echo $head | cut -d ":" -f 3`

                    echo "[$0] TSD:$TSD >> empty_site:$empty_site >> FLANK_SIZE:$FLANK_SIZE >> position_TE_SV:$position_TE_SV"
                    echo "[$0] TSD:$TSD >> empty_site:$empty_site >> FLANK_SIZE:$FLANK_SIZE >> position_TE_SV:$position_TE_SV" >&2

                    ( test -n "$TSD" && test -n "$empty_site" && test -n "$FLANK_SIZE" && test -n "$position_TE_SV" && \
                        python3 ${path_this_script}/find_svi.py $TSD $empty_site $FLANK_SIZE $position_TE_SV >> ${OUTPUT} ) || echo "[$0] ERROR : MISSING ARGS FOR find_svi.py " >&2;
                    rm -f ${path_out_dir}/tmp_${id}.bed
                fi;

                else
                    echo "[$0] **ERROR** : FILE NOT FOUND FOR ID : $id";
                fi;
            # exit 1;
        else
            reads=${REPO_READS_FA}/${name}fasta

            echo 
            echo "[$0] ------TSD------"
            id=`echo $reads | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*"`
            echo "[$0] reads : "$reads
            echo "[$0] id : "$id
            echo "[$0] find_file : ${FIND_FA}"
            head=`grep ".*:.*:$id:[0-9]*:[PI]" ${FIND_FA}`

            echo "[$0] head : "$head

            #awk -v var=$head 'BEGIN {nb=0} { if( var == $0 || nb == 1 ){print $0; nb = nb + 1;} }' "${FIND_FA}" > ${path_out_dir}/sequence_TE_${id}.fasta
            grep -w $head -A 1 "${FIND_FA}" > ${path_out_dir}/sequence_TE_${id}.fasta 
            # cat ${path_out_dir}/sequence_TE_${id}.fasta 
            
            if ! test -s ${path_out_dir}/sequence_TE_${id}.fasta; then
                echo "[$0] **ERROR** : can't get sequence TE in ${FIND_FA}" ;
                exit 1 ;
            fi;

            #Warning : WORK ONLY BECAUSE ERROR IN LONG READS
            echo "[$0] *********BLAST 1 TE VS READ**********"
            echo "[$0] --READS : $reads"
            makeblastdb -in "$reads" -dbtype nucl
            blastn -db "$reads" -query ${path_out_dir}/sequence_TE_${id}.fasta \
                -perc_identity 100 \
                -outfmt 6 \
                -out ${path_out_dir}/sequence_TE_${id}.bln ;

            echo "[$0] --BEGIN sequence_TE_${id}.bln--"
            cat ${path_out_dir}/sequence_TE_${id}.bln
            echo "[$0] --END sequence_TE_${id}.bln--"

            # TODO if SIZE FLANK is too big, 
            # get flank bed file
            awk -v var=$FLANK_SIZE '{if($9 - var > 0){ if($9 < $10){ print $2"\t"$9-var-1"\t"$9-1"\n" $2"\t"$10"\t"$10+var }else{ print $2"\t"$10-var-1"\t"$10-1"\n" $2"\t"$9"\t"$9+var  } } }' ${path_out_dir}/sequence_TE_${id}.bln > ${path_out_dir}/flank_TE_${id}.bed
            # get TE SEQ
            awk '{ if($9 < $10){ print $2"\t"$9-1"\t"$10"\t" "forward" "\t" "1" "\t" "+" }else{ print $2"\t"$10-1"\t"$9"\t" "reverse" "\t" "1" "\t" "-" } }' ${path_out_dir}/sequence_TE_${id}.bln > ${path_out_dir}/sequence_TE_${id}.bed

            bedtools getfasta -fi $reads -bed ${path_out_dir}/flank_TE_${id}.bed > ${path_out_dir}/flank_TE_${id}.fasta
            bedtools getfasta -fi $reads -bed ${path_out_dir}/sequence_TE_${id}.bed -name+ > ${path_out_dir}/sequence_TE_${id}.fasta

            echo "[$0] SHOW FLANK : ${path_out_dir}/flank_TE_${id}.bed"
            cat ${path_out_dir}/flank_TE_${id}.bed
            echo "[$0] ..SHOW FLANK END.."

            echo "[$0] *********BLAST 2 TE VS DBTE**********"
            echo "[$0] ${DB_TE}"
            #makeblastdb -in ${DB_TE} -dbtype nucl ##comment for threads task
            blastn -db ${DB_TE} -query ${path_out_dir}/sequence_TE_${id}.fasta -outfmt 6 -out ${path_out_dir}/TE_vs_databaseTE_${id}.bln

            echo "[$0] head ${path_out_dir}/TE_vs_databaseTE_${id}.bln" >&2; 
            head -n 1 ${path_out_dir}/TE_vs_databaseTE_${id}.bln >&2
            echo "" >&2
            
            strand=`awk 'NR==1 {if ($9 < $10){print "+"} else {print "-"}}' ${path_out_dir}/TE_vs_databaseTE_${id}.bln`
            echo $reads >> ${OUTPUT}
            echo $head  >> ${OUTPUT}

            # FIND TSD
            # Warning : path
            echo "[$0] flank_TE_${id}.fasta:"`test -s "${path_out_dir}/flank_TE_${id}.fasta" && echo "OK" || echo "NONE"`"   sequence_TE.fasta:"`test -s "${path_out_dir}/sequence_TE_${id}.fasta" && echo "OK" || echo "NONE"`"  FK_S:$FLANK_SIZE     id:$id     strand:$strand     TSD_SIZE:$TSD_SIZE"
            echo "";
            echo "";
            echo "[$0] WHY params ${path_out_dir}/flank_TE_${id}.fasta ${path_out_dir}/sequence_TE_${id}.fasta $FLANK_SIZE $id $strand $TSD_SIZE" >&2
            
            if test -s "${path_out_dir}/TE_vs_databaseTE_${id}.bln" && test -n "$FLANK_SIZE" && test -n "$id" && test -n "$strand" && test -n "$TSD_SIZE"; then
                python3 ${path_this_script}/find_tsd.py ${path_out_dir}/flank_TE_${id}.fasta ${path_out_dir}/sequence_TE_${id}.fasta $FLANK_SIZE $id $strand $TSD_SIZE >> ${OUTPUT}
                OK=`grep "$head" ${OUTPUT} -A 2 | grep -o "OK" || echo ""`
            else
                OK="";
            fi;
            
            if [ -n "$OK" ]; then
                TSD=`grep "$head" ${OUTPUT} -A 4 | grep "++:[A-Z]*:++" -o | head -n 1 | grep -o "[A-Z]*"`
                echo "$head" | awk -v FLANK_SIZE="$FLANK_SIZE"  -F ":" 'OFS="\t"{print $1, $3-FLANK_SIZE, $3+FLANK_SIZE}' | tr -d ">" >> ${path_out_dir}/tmp_${id}.bed
                empty_site=`bedtools getfasta -fi ${GENOME} -bed ${path_out_dir}/tmp_${id}.bed | grep -v ">"`
                position_TE_SV=`echo $head | cut -d ":" -f 3`

                echo "[$0] TSD:$TSD >> empty_site:$empty_site >> FLANK_SIZE:$FLANK_SIZE >> position_TE_SV:$position_TE_SV"
                echo "[$0] TSD:$TSD >> empty_site:$empty_site >> FLANK_SIZE:$FLANK_SIZE >> position_TE_SV:$position_TE_SV" >&2

                ( test -n "$TSD" && test -n "$empty_site" && test -n "$FLANK_SIZE" && test -n "$position_TE_SV" && \
                    python3 ${path_this_script}/find_svi.py $TSD $empty_site $FLANK_SIZE $position_TE_SV >> ${OUTPUT} ) || echo "[$0] ERROR : MISSING ARGS FOR find_svi.py " >&2;
                rm -f ${path_out_dir}/tmp_${id}.bed
            fi;

            rm -f ${path_out_dir}/TE_vs_databaseTE_${id}.bln
            rm -f ${path_out_dir}/sequence_TE_${id}.fasta \
                ${path_out_dir}/sequence_TE_${id}.bln \
                ${path_out_dir}/sequence_TE_${id}.bed
            rm -f ${path_out_dir}/flank_TE_${id}.fasta \
                ${path_out_dir}/flank_TE_${id}.bed
        fi;
done;

number_ok=`grep "OK" ${OUTPUT} -c`
number_ko=`grep "KO" ${OUTPUT} -c`
number_total=`grep "reads" ${OUTPUT} -c`
number_k_o=`grep "K-O" ${OUTPUT}  -c`

if [ $number_total -ne 0 ]; then

    #RESUME
    echo "OK/total : $number_ok/$number_total" >> ${OUTPUT}
    echo "KO/total : $number_ko/$number_total" >> ${OUTPUT}
    echo "OK+KO/total : $(($number_ok+$number_ko))/$number_total" >> ${OUTPUT}
    echo "K-O/total : $number_k_o/$number_total" >> ${OUTPUT}
    echo "OK+K-O/total : $(($number_ok+$number_k_o))/$number_total" >> ${OUTPUT}
    echo "OK% : $(($number_ok*100/$number_total))%" >> ${OUTPUT}

else
    echo "ERROR-404: NUMBER TOTAL OF TE EQUAL 0 " >> ${OUTPUT}
fi;