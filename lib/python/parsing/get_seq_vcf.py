#!/usr/bin python3
# -*- coding: utf-8 -*-

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
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Mourdas Mohamed
#                                                                                                                                   
####################################################################################################################################


"""
    GET SEQUENCE REPORT
    ==========================
    :author:  Mourdas MOHAMED 
    :contact: mourdas.mohamed@igh.cnrs.fr
    :date: 01/06/2020
    :version: 0.1
    Script description
    ------------------
    get_seq_vcf.py get the sequences report for INDEL in vcf file 
    -------
    >>> get_seq_vcf.py file.vcf 
    
    Help Programm
    -------------
    usage: get_seq_vcf.py [-h] [-m MIN_LEN_SEQ] [-t TYPE] [-c CHROM]
                      vcf_file fasta_out

    positional arguments:
      vcf_file              VCF file to parse
      fasta_out             name output fasta file

    optional arguments:
      -h, --help            show this help message and exit
      -m MIN_LEN_SEQ, --min_len_seq MIN_LEN_SEQ
                            minimum size of the sequence to keep [1000]
      -t TYPE, --type TYPE  not keep this type on vcf ["<DEL>"]
      -c CHROM, --chrom CHROM
                            chromosome (or part) to keep [2L,2R,3L,3R,^4_,X_]
"""


import sys
import os
import re
import argparse

parser = argparse.ArgumentParser(description="Get the sequences report for INDEL in vcf file (format VCFv4.3)")

#MAIN ARGS
parser.add_argument("vcf_file", type=str,
                    help="VCF file to parse")
parser.add_argument("fasta_out", type=str,
                    help="name output fasta file")

#OPTION
parser.add_argument("-m", "--min_len_seq", type=int, default=1000,
                    help="minimum size of the sequence to keep (default: [1000])")

parser.add_argument("--max-len-seq", dest='max_len_seq', type=int, default=-1,
                    help="maximum size of the sequence to keep (default: [-1])")

parser.add_argument("-t", "--type", type=str, default='<DEL>',
                    help="not keep this type on vcf (give a list of arguments separate the values with commas \"<DEL>,<INS>\") (default: [\"<DEL>\"])")

#Warning : format of the list of chromosome must be separate by "," it could be regex
parser.add_argument("-c", "--chrom", type=str, default='^2L_,^2R_,^3L_,^3R_,^4_,^X_',
                    help='chromosome (or part/contig) to keep (give a list of arguments separate the values with commas "X,Y") put \".\" for keep all chromosome (default: [2L,2R,3L,3R,^4_,X_])')

parser.add_argument("-k", "--keep", default=False, action='store_true',
                    help='keep the sequence contain only character N (default: False)')

parser.add_argument("-i", "--idrs", default=False, action='store_true',
                    help='for make all id reads support')


# print("usage example : python3 get_seq_vcf.py file.vcf outpout.fasta --keep --chrom \".\" " )
# print("usage example : python3 get_seq_vcf.py file.vcf outpout.fasta --keep --chrom chrX,chrY " )
# print("usage example : python3 get_seq_vcf.py file.vcf outpout.fasta \n" )
args = parser.parse_args()

#file in
file = open(args.vcf_file, "r")

#file out
file_fasta  = open(args.fasta_out, "w")

type_list   = args.type.split(",")#NOT KEEP THIS
chrom_list  = args.chrom.split(",")#KEEP ONLY
min_len_seq = args.min_len_seq
max_len_seq = args.max_len_seq
keep_seq_N  = args.keep
ID_RS       = args.idrs

print("[" + str(sys.argv[0]) + "] : exclude type : ", type_list)
print("[" + str(sys.argv[0]) + "] : regex contig liste : ", chrom_list)


#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False

line = file.readline()

version_vcf =  line.split("=")[1].strip()
print("Version VCF :", version_vcf)
##fileformat=VCFv4.3
#Check format vcf file
if version_vcf != "VCFv4.3" and version_vcf != "VCFv4.2" and version_vcf != "VCFv4.1":
    print("[" + str(sys.argv[0]) + "] : ERROR format vcf must be VCFv4.3 or VCFv4.2 not " + str(line.split("=")[1].strip()) )
    print("[" + str(sys.argv[0]) + "] : please change format of vcf file, or use Snifflesv1.0.10 for genrate the good vcf file")
    exit(1)


contig_vcf = []
contigs_vcf_org = []
while line and line[:2] == "##":
    #line.split("=")[2].split(",")[1] == chrom
    if "contig=" in line:
        contig = line.split("=")[2].split(",")[0]
        contigs_vcf_org.append(contig)

        if regex_in_list(contig, chrom_list):
            contig_vcf.append(contig)
        
    line = file.readline()


print("[" + str(sys.argv[0]) + "] : contig keeping = " + str(contig_vcf))
if len(contig_vcf) == 0:
    print("[" + str(sys.argv[0]) + "] : ERROR your regex chrom d'osnt match with any chromosome in your vcf file")
    print("[" + str(sys.argv[0]) + "] : ERROR please change your regex (-c option) > regex=" + str(args.chrom) + " :: list=" + str(chrom_list))
    print("[" + str(sys.argv[0]) + "] : ERROR your regex must match on " + str(contigs_vcf_org))
    exit(2)

##Sniffles
count = 0
if version_vcf == "VCFv4.3":
    
    while line :
        
        seq = None
        if re.search("^[^#]", line) :
            #FOR HEADER
            spl     = line.split("\t")
            chrom   = spl[0]
            start   = spl[1]
            ID      = spl[2]
            type_v  = spl[4]
            
            if re.search("END=([0-9]+)", spl[7]) != None and re.search("RE=([0-9]+)", spl[7]) != None :
                end     = re.search("END=([0-9]+)", spl[7]).group(1)
                precise = spl[7].split(";")[0]
                read_support = re.search("RE=([0-9]+)", spl[7]).group(1)
            
                if "SEQ=" in spl[7] :
                    seq    = re.search("SEQ=([A-Z][A-Z]+);", spl[7]).group(1)
                    #seq    = spl[7].split(";")[-2].split("=")[1]
                    seq_NN = re.search("^[N]+$", seq)#Not keep sequence contains N (pptional)
                    seq    = seq.replace("N", "A")#replace

                condition = seq != None and (not seq_NN or keep_seq_N) and min_len_seq < len(seq) and type_v not in type_list and ( (max_len_seq != -1 and len(seq) <= max_len_seq) or (max_len_seq == -1) )
                if regex_in_list(chrom, chrom_list) and condition :
                    #EXEMPLE  FORMAT
                    #2R:<INS>:1:190:sniffles.INS.2:1:PRECISE
                    ID = "sniffles." + type_v.replace("<", "").replace(">", "") + "." + ID
                    file_fasta.write(">" + ":".join([chrom, type_v, start, end, ID, read_support, precise]) + "\n" + seq + "\n")
                    count += 1

        line = file.readline()


if version_vcf == "VCFv4.2":
    while line :
        
        if re.search("^[^#]", line) :
            #FOR HEADER
            spl     = line.split("\t")
            chrom   = spl[0]
            ID      = spl[2]
            start   = spl[1]
            type_v  = spl[4]
            
            if re.search("END=([0-9]+)", spl[7]) != None and re.search("SUPPORT=([0-9]+)", spl[7]) != None :
                end     = re.search("END=([0-9]+)", spl[7]).group(1)
                precise = "PRECISE"
                read_support = re.search("SUPPORT=([0-9]+)", spl[7]).group(1)
                
                if "SEQS=" in spl[7] :

                    seqs   = re.search("SEQS=([A-Z][A-Z,]+);", spl[7]).group(1).split(",")

                    for i, seq in enumerate(seqs) :
                        
                        if ID_RS :
                            ID_VR = ID + "." + str(i)
                        else:
                            ID_VR = ID

                        seq_NN = re.search("^[N]+$", seq)#Not keep sequence contains N (pptional)
                        seq    = seq.replace("N", "A")#replace
                    
                        condition = seq != None and (not seq_NN or keep_seq_N) and min_len_seq < len(seq) and type_v not in type_list and ( (max_len_seq != -1 and len(seq) <= max_len_seq) or (max_len_seq == -1) )
                        if regex_in_list(chrom, chrom_list) and condition :
                            #EXEMPLE  FORMAT
                            #2R:<INS>:1:190:svim.INS.2:1:PRECISE
                            file_fasta.write(">" + ":".join([chrom, type_v, start, end, ID_VR, read_support, precise, str(i)]) + "\n" + seq + "\n")
                            count += 1

        line = file.readline()

##Sniffles
if version_vcf == "VCFv4.1":
    while line :
        
        seq = None
        if re.search("^[^#]", line) :
            #FOR HEADER
            spl     = line.split("\t")
            chrom   = spl[0]
            start   = spl[1]
            ID      = spl[2]
            seq     = "" 
            
            if re.search("END=([0-9]+)", spl[7]) != None and re.search("RE=([0-9]+)", spl[7]) != None :
                end     = re.search("END=([0-9]+)", spl[7]).group(1)
                precise = spl[7].split(";")[0]
                read_support = re.search("RE=([0-9]+)", spl[7]).group(1)
                
                type_v  = re.search("SVTYPE=([A-Z]+)", spl[7]).group(1)

                if type_v in "INS":
                    seq     = spl[4]
                elif type_v in "DEL":
                    seq     = spl[3]

                #seq     = spl[4]
                #seq    = spl[7].split(";")[-2].split("=")[1]
                seq_NN = re.search("^[N]+$", seq)#Not keep sequence contains N (pptional)
                seq    = seq.replace("N", "A")#replace

                condition = seq != None and (not seq_NN or keep_seq_N) and min_len_seq < len(seq) and type_v not in type_list and ( (max_len_seq != -1 and len(seq) <= max_len_seq) or (max_len_seq == -1) )
                if regex_in_list(chrom, chrom_list) and condition :
                    #EXEMPLE  FORMAT
                    #2R:<INS>:1:190:sniffles.INS.2:1:PRECISE
                    ID = "sniffles." + type_v.replace("<", "").replace(">", "") + "." + ID
                    file_fasta.write(">" + ":".join([chrom, "<" + type_v + ">", start, end, ID, read_support, precise]) + "\n" + seq + "\n")
                    count += 1

        line = file.readline()

print("[" + str(sys.argv[0]) + "] : Total sequence found=" + str(count))
file.close()
file_fasta.close()
