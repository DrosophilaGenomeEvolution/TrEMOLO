#!/bin/bash

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
        OUTPUT="$2"
        shift # past argument
        shift # past value
        ;;
        -i|--input)
        INPUT="$2"
        shift # past argument
        shift # past value
        ;;
        -c|--chrom)
        CHROM="$2"
        shift # past argument
        shift # past value
        ;;
        -g|--genome)
        GENOME="$2"
        shift # past argument
        shift # past value
        ;;
        -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
        *)
        POSITIONAL_ARGS+=("$1") # save positional arg
        shift # past argument
        ;;
    esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters


if ! [[ -n "$INPUT" ]] || ! test -s $INPUT;
then
    echo "[$0] ERROR : need --input"
    echo
    echo "[$0] Usage : buildFrequencyGenerations.sh -i <input-init-file> [-o OUTPUT-NAME-DIRECTORY] [-g GENOME-FASTA-FILE] [-c REGEX-CHROM]"
    exit 1;
fi;

if [[ -n "$GENOME" ]] && ! test -s $GENOME;
then
    echo
    echo "[$0] ERROR : -g|--genome ${GENOME} No sush file or directory"
    echo
    echo "[$0] Usage : buildFrequencyGenerations.sh -i <input-init-file> [-o OUTPUT-NAME-DIRECTORY] [-g GENOME-FASTA-FILE] [-c REGEX-CHROM]"
    exit 1;
fi;


! test -n "${OUTPUT}" && OUTPUT="SCATTER-FREQ-TE-TrEMOLO";
! test -n "$CHROM" && CHROM=".";

echo "[$0] BEGIN"
echo "---------ARGS--------"
echo "INPUT PATH GENERATIONS     = ${INPUT}"
echo "OUTPUT NAME                = ${OUTPUT}"
echo "REGEX CHROMOSOME           = ${CHROM}"
test -n "$GENOME" && \
    echo "GENOME                     = ${GENOME}"

echo

path_to_module=`dirname "$0"`
echo "   PATH : $path_to_module";
mkdir -p "${OUTPUT}"/web ;
mkdir -p "${OUTPUT}"/work;

##FORMAT INPUT REQUIRED
#work_directory_path:G4
#-----------------------
#or
#-----------------------
#work_directory_path
grep ":" "$INPUT" 2> /dev/null
if [ $? -eq 0 ]; then
    grep ":" "$INPUT" | awk -F ":" '{print $0"\t"substr($2, 2, length($2))}' | \
        sort -k 2 -n | cut -f 1 > "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt" 2>/dev/null 
else
    echo "INIT GENERATION" && sed 's/\/$//g' ${INPUT}  | awk '{print $1":G"NR}' \
        > "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt"
fi;

echo "NB GENERATIONS : "$(cat "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt" | wc -l)
echo 

echo "PREPARE DATA SCATTER..."
echo -e "chrom\tx\ty1\ty2\tgroup\tname\tID" > "${OUTPUT}/work/SCATTER.csv"
for work in `cat "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt"`
do
    dir=`echo $work | cut -d ":" -f 1`;
    name=`echo $work | cut -d ":" -f 2`;
    echo "   generation : $name; WORK_DIRECTORY : $dir"
    ! test -s ${dir}/TE_INFOS.bed && echo "[$0] ERROR : FILE NOT FOUND ${dir}/TE_INFOS.bed" && exit 1
    ! test -s ${dir}/TE_FREQUENCY_TrEMOLO.bed && echo "BUILD FILE TE_FREQUENCY_TrEMOLO..." && \
        grep -E "${CHROM}" ${dir}/TE_INFOS.bed | grep -v -E "HARD|SOFT" | grep OUTSIDER | awk '$11!="NONE" && $12!="NONE"' | grep -w -v "DEL" > ${dir}/TE_FREQUENCY_TrEMOLO.bed;

    mkdir -p "${OUTPUT}/work/${name}";
    test -s ${dir}/TE_FREQUENCY_TrEMOLO.bed && echo "   INFOS NB TE : "$(wc -l ${dir}/TE_FREQUENCY_TrEMOLO.bed)
    test -s ${dir}/TE_FREQUENCY_TrEMOLO.bed && \
        awk 'substr($0, 1, 1)!="#" && OFS="\t"{print $1, $2, $3, $4, $11, $12}' ${dir}/TE_FREQUENCY_TrEMOLO.bed > "${OUTPUT}/work/${name}/POS_TE_LIFT.bdg"

    awk -v name="$name" 'OFS="\t"{split($4, sp, "|"); print $1, $2, $5/100, $6/100, name, $1":"$2":"sp[1], sp[2]}' "${OUTPUT}/work/$name/POS_TE_LIFT.bdg" | tr "," "." >> "${OUTPUT}/work/SCATTER.csv";
done;

wc -l ${OUTPUT}/work/SCATTER.csv
echo 

echo -e "chrom\tx\ty1\ty2\tgroup\tname\tIN_OUT\tID\ttype" > "${OUTPUT}/work/FT1_SCATTER.csv"
awk 'NR>1 && OFS="\t"{split($7, sp, "." ); print $1, $2, $3, $4, $5, $6, "OUTSIDER", $7, sp[2]}' "${OUTPUT}/work/SCATTER.csv" >> "${OUTPUT}/work/FT1_SCATTER.csv"

wc -l ${OUTPUT}/work/FT1_SCATTER.csv
echo 

echo "MERGE POSITIONS..."
awk 'NR>1 && OFS="\t"{print $1, $2, $2+1, $3, $4, $5, $6, $7, $8, $9}' "${OUTPUT}"/work/FT1_SCATTER.csv | bedtools sort | bedtools cluster -d 200 | awk '
BEGIN{
    clust="";
    pos="";
    print "chrom\tx\ty1\ty2\tgroup\tname\tIN_OUT\tID\ttype";
} 

OFS="\t" {
    split($7, sp, ":"); 
    TE=sp[3]; 
    if(clust != $11){
        pos=$2; 
        print $1, pos, $4, $5, $6, $1":"pos":"TE, $8, $9, $10;
        clust=$9;
        pos=$2;
    }
    else{
        print $1, pos, $4, $5, $6, $1":"pos":"TE, $8, $9, $10;
    }  
}' > "${OUTPUT}/work/SCATTER_MERGE.csv"

wc -l "${OUTPUT}/work/SCATTER_MERGE.csv"
echo

echo "BUILD METRIC..."
python3 "$path_to_module/metric.py" "${OUTPUT}/work/SCATTER_MERGE.csv" "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt" > "${OUTPUT}/work/METRIC.csv";    
echo "SORTING..."
python3 "$path_to_module/putMetricScatter.py" "${OUTPUT}/work/METRIC.csv" "${OUTPUT}/work/SCATTER_MERGE.csv" "${OUTPUT}/work/SCATTER_SORT.csv"

echo "BUILD DATA JS SCATTER..."
##BUILD JS SCATTER
awk 'NR>1' "${OUTPUT}/work/SCATTER_SORT.csv" | sort -k 6 | awk 'BEGIN{print "chrom\tx\ty1\ty2\tgroup\tname\ttype1\ttype2\tIN_OUT \tID"} {print $0}' | \
  grep -w -v "NONE" | awk '
BEGIN{
    print "allData = {"
    gen="";
    name="";
    chrom="";
    id="";
    type1="";
    type2="";
    IN_OUT="";
} 

NR>1 {
    if( chrom != $1 ){
        if( chrom != "" ){
            print "\"type1\": \""type1"\",";
            print "\"type2\": \""type2"\",";
            print "\"id\": \""id"\",";
            print "\"name\": \""name"\"\n}],\n";
        }

        chrom=$1;

        print "\""chrom"\": [{";
        print "\""$5"x\": "$2",";
        print "\""$5"y1\": "$3",";
        print "\""$5"y2\": "$4",";
        print "\""$5"io\": \""$9"\",";
        gen=$5;
        id=$6;
        type1=$7;
        type2=$8;
        IN_OUT=$9;
        split($6, sp, ":");
        name=sp[3];
    }
    else if( id != $5 ){
        print "\"type1\": \""type1"\",";
        print "\"type2\": \""type2"\",";
        print "\"id\": \""id"\",";
        print "\"name\": \""name"\"";
        print "}, {\n""\""$5"x\": "$2",";
        print "\""$5"y1\": "$3",";
        print "\""$5"y2\": "$4",";
        print "\""$5"io\": \""$9"\",";
        gen=$5;
        id=$6;
        type1=$7;
        type2=$8;
        IN_OUT=$9;
        split($6, sp, ":");
        name=sp[3];
    }
    else {
        print "\""$5"x\": "$2",";
        print "\""$5"y1\": "$3",";
        print "\""$5"y2\": "$4",";
        print "\""$5"io\": \""$9"\",";
    }
}

END{
    print "\"type1\": \""type1"\",";
    print "\"type2\": \""type2"\",";
    print "\"id\": \""id"\",";
    print "\"name\": \""name"\"\n}]};";
}' > "${OUTPUT}/web/data.js"


awk 'NR>1 {print $1}' "${OUTPUT}/work/SCATTER_SORT.csv" | sort -u > "${OUTPUT}/work/ID_CHROM.txt"

pathG1=`awk -F ":" 'NR==1 {print $1}'  "${OUTPUT}/work/INIT_FREQ_TE_TrEMOLO.txt"`;

echo "GETTING SIZE OF CHROMOSOMES..."

if ! test -s "$GENOME"; then
    samtools view -h `realpath ${pathG1}`/OUTSIDER/FREQUENCY/MAPPING_POSTION_TE.bam | awk 'substr($0, 1, 1)=="@"{print $0; next;} {exit;}' | grep "^@SQ" | \
        grep -w -f "${OUTPUT}"/work/ID_CHROM.txt | \
        awk '
        BEGIN {
            print "var chroms = {"
        }

        {
            split($2, sn, ":"); 
            split($3, ln, ":"); 
            print "   \""sn[2]"\":", ln[2]","
        }

        END{
            print "}"
        }' > "${OUTPUT}"/web/chrom.js
else
    echo "   BUILD FAIdx : ${GENOME} ...";
    samtools faidx "$GENOME";
    
    echo "   BUILD JS...";
    cat "${GENOME}.fai" | \
        grep -w -f "${OUTPUT}"/work/ID_CHROM.txt | \
        awk '
        BEGIN {
            print "var chroms = {"
        }

        {
            print "   \""$1"\":", $2","
        }

        END{
            print "}"
        }' > "${OUTPUT}"/web/chrom.js
fi;




echo "GETTING GENERATIONS..."
##GENERATIONS
cat  "${OUTPUT}"/work/INIT_FREQ_TE_TrEMOLO.txt | awk -v OUTPUT="${OUTPUT}" -F ":" '
    BEGIN{
        print "var generations = [" > OUTPUT"/web/generations.js";
    }

    NR==1 {
        print "\""$2"\"" >> OUTPUT"/web/generations.js"
    }

    NR>1 {
        print ",\""$2"\"" >> OUTPUT"/web/generations.js"
    }

    END{
        print "];" >> OUTPUT"/web/generations.js"
    }
'

cp -f $path_to_module/GRAPH.html "${OUTPUT}"/web/
cp -f $path_to_module/index.html "${OUTPUT}"/

echo "[$0] END"
