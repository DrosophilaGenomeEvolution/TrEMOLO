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
                            chromosome (or part) to keep [2L,2R,3L,3R,4,X]
"""


import sys
import os
import re
import argparse

parser = argparse.ArgumentParser(description="get the sequences report for INDEL in vcf file")

#MAIN ARGS
parser.add_argument("vcf_file", type=str,
                    help="VCF file to parse")
parser.add_argument("fasta_out", type=str,
                    help="name output fasta file")

#OPTION
parser.add_argument("-m", "--min_len_seq", type=int, default=1000,
                    help="minimum size of the sequence to keep [1000]")
parser.add_argument("-t", "--type", type=str, default='<DEL>',
                    help="not keep this type on vcf (give a list of arguments separate the values ​​with commas \"<DEL>,<INS>\") [\"<DEL>\"]")
parser.add_argument("-c", "--chrom", type=str, default='2L,2R,3L,3R,4,X',
                    help='chromosome (or part) to keep (give a list of arguments separate the values ​​with commas "X,Y") [2L,2R,3L,3R,4,X]')

args = parser.parse_args()

#file in
file = open(args.vcf_file, "r")

#file out
file_fasta  = open(args.fasta_out, "w")

type_list   = args.type.split(",")#NOT KEEP THIS
chrom_list  = args.chrom.split(",")#KEEP ONLY
min_len_seq = args.min_len_seq

print("exclude type : ", type_list)
print("chromosome liste : ", chrom_list)

#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False

line = file.readline()
while line :
    
    seq = None
    if re.search("^[^#]", line) :
        #FOR HEADER
        spl     = line.split("\t")
        chrom   = spl[0]
        ID      = spl[2]
        type_v  = spl[4]
        start   = spl[1]
        end     = spl[7].split(";")[3].split("=")[1]
        precise = spl[7].split(";")[0]
        read_support = spl[7].split(";")[-1].split("=")[1]
        
        
        if spl[7].split(";")[-2].split("=")[0] == "SEQ":
            seq    = spl[7].split(";")[-2].split("=")[1]
            seq_NN = re.search("^[N]+$", seq)#Not keep sequence contains N (pptional)
            seq    = seq.replace("N", "A")#replace
        

        condition = seq != None and not seq_NN and min_len_seq < len(seq) and type_v not in type_list
        if regex_in_list(chrom, chrom_list) and condition :
            #EXEMPLE  FORMAT
            #2R:<INS>:1;190;2:1:PRECISE
            file_fasta.write(">" + ":".join([chrom, type_v, start, end, ID, read_support, precise]) + "\n" + seq + "\n")

    line = file.readline()

file.close()
file_fasta.close()
