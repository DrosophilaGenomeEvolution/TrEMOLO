import pysam
import os
import sys
import re
import argparse
import numpy

#TODO CHECK ALL METHODS
parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion (HARD reads), at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")


#OPTION
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=30,
                    help="minimum size of sequence.")
args = parser.parse_args()

name_bamfile = args.bam_file

min_size    = args.min_size

bamfile  = pysam.AlignmentFile(name_bamfile, "rb")

ID       = 0

tab_doublons  = []
ref_name_pred = ""
dico_clust    = []
for e, read in enumerate(bamfile.fetch()):
    start_query     = read.query_alignment_start
    reference_start = read.reference_start
    seq             = read.seq
    read_name       = read.query_name
    REF             = read.reference_name #chrom

    count_ref  = 0 #number of nucleotides on the ref before reaching the insertion site
    count_read = 0
    count_read_real = 0

    if read_name :

        for tupl in read.cigartuples :
            
            if tupl[0] in [0, 2, 7] : #Check M,D,= CIGAR for position on ref
                count_ref  += tupl[1]

            if tupl[0] in [0, 1, 7, 4] : #Check M,I,=,S CIGAR for postion on reads
                count_read += tupl[1]

            if tupl[0] in [0, 1, 7, 4, 5] : #Check M,I,=,S,H CIGAR for postion on reads real (HARD)
                count_read_real += tupl[1]

            #if we have found HARD to a good position
            if tupl[0] == 5 and tupl[1] >= min_size :

                identify = ":".join([REF, str(count_ref + reference_start), str(read_name), str(count_read_real), "L", str(len(seq))])
                if seq and tupl[1] > min_size and identify not in tab_doublons :

                    side = ""
                    if count_read_real - tupl[1] == 0 :
                        side = "L"
                    else :
                        side = "R"

                    if side == "L" :
                        dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq), "SEQ":seq})
                    else :
                        dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq), "SEQ":seq})
                    
                    tab_doublons.append(identify)

                    ID += 1

print("\t".join(["#REF", "START", "END","ID", "ID_READ", "POS_REAL_READ", "FLAG", "SIDE", "SIZE_SEQ", "SEQ"]))
for indice, dic in enumerate(dico_clust) :
    print("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "HARD." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]))






