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

# FIND FQ AND CHANGE TO FASTA
# ==========================
# :author:  Mourdas MOHAMED 
# :contact: mourdas.mohamed@igh.cnrs.fr
# :date: 01/06/2020
# :version: 0.1
# Script description
# ------------------
# find_fq_to_fasta.sh Convert fastq file of all support READS TE to fasta file
# -------
# bash find_fq_to_fasta.sh prefix_find_ZAM.fasta READS_SUPPORT/ OUT_DIR/
#
# Help Programm
# -------------
# usage: find_fq_to_fasta.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <NAME_OUT_DIRECTORY>
# **WARNING: this program needs the fastq_to_fasta.py script in the same directory**


#ARGS
#File find TE fasta 
FIND_FA=$1
#File contains READS support for find TE
REPO_READS=$2
#Name out directory
OUT_DIR=$3

path_this_script=`dirname $0`
echo $path_this_script ;

if [ "$#" -ne 3 ]; then
    echo "[$0] ERROR : need 3 arguments. You have put $# arguments" ;
    echo "[$0] usage : find_fq_to_fasta.sh <file_find_nameTE.fasta> <DIRECTORY_READS_SUPPORT> <NAME_OUT_DIRECTORY>"
    exit 1 ;
fi


! test -s "$FIND_FA" && echo "[$0] ERROR : $FIND_FA d'osent exist" && exit 2; 
! test -d "$REPO_READS" && echo "[$0] ERROR : $REPO_READS is not a directory" && exit 3;

REPO_READS=`echo ${REPO_READS} | sed 's/[/]$//g'`

mkdir -p ${OUT_DIR}
OUT_DIR=`echo ${OUT_DIR} | sed 's/[/]$//g'`


nombre_element=`grep ">" ${FIND_FA} | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*" | wc -l`
i=0
for id in `grep ">" ${FIND_FA} | grep -o "[0-9]:[A-Za-z\.0-9]*:[0-9]*:[PI]" | grep -o ":[A-Za-z\.0-9]*:" | grep -o "[A-Za-z\.0-9]*"`; do
    fr="`ls ${REPO_READS} | grep ":$id:"`"
    echo "[$0] id : $id"
    echo "[$0] file : $fr"
    name=`echo $fr | grep -o ".*\."`
    echo "[$0] name : $name"
    test -f ${REPO_READS}/$fr || echo "[$0] ERROR FILE ${REPO_READS}/$fr NOT FOUND";
    #Warning : path
    test -f ${REPO_READS}/$fr && \
        python3 ${path_this_script}/fastq_to_fasta.py ${REPO_READS}/$fr ${OUT_DIR}/${name}fasta;
    i=$(($i + 1));
    echo "[$0] $i/$nombre_element";
done;

