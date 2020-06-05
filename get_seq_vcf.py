import sys
import os
import re
import argparse

parser = argparse.ArgumentParser()

#MAIN ARGS
parser.add_argument("vcf_file", type=str,
                    help="VCF file to parse")
parser.add_argument("fasta_out", type=str,
                    help="name fasta file out")

#OPTION
parser.add_argument("-m", "--min_len_seq", type=int, default=1000,
                    help="size seq minimum to keep [1000]")

args = parser.parse_args()

#file in
file  = open(args.vcf_file , "r")

#file out
file_fasta = open(args.fasta_out, "w")


tab_type  = ["<DEL>"]#NOT KEEP THIS
chrom_tab = ["2L", "2R", "3L", "3R", "4", "X"]
min_len_seq = args.min_len_seq


line = file.readline()
while line :
    
    seq = None
    if re.search("^[^#]", line) :
        #FOR HEADER
        spl     = line.split("\t")
        chrom   = spl[0].split("_")[0]
        ID      = spl[2]
        type_v  = spl[4]
        start   = spl[1]
        end     = spl[7].split(";")[3].split("=")[1]
        precise = spl[7].split(";")[0]
        read_support = spl[7].split(";")[-1].split("=")[1]
        
        
        if spl[7].split(";")[-2].split("=")[0] == "SEQ":
            seq    = spl[7].split(";")[-2].split("=")[1]
            seq_NN = re.search("^[N]+$", seq)#Not keep sequence contains N
            seq    = seq.replace("N", "A")
        

        condition = seq != None and not seq_NN and min_len_seq < len(seq) and type_v not in tab_type
        if chrom.split("_")[0] in chrom_tab and condition :
            #EXEMPLE  FORMAT
            #2R:<INS>:1;190;2:1:PRECISE
            file_fasta.write(">" + ":".join([chrom, type_v, start, end, ID, read_support, precise]) + "\n" + seq + "\n")

    line = file.readline()


file.close()
file_fasta.close()
