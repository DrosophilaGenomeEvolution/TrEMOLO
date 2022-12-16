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
parser.add_argument("-m", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=30,
                    help="minimum size of sequence.")
parser.add_argument("-c", "--chrom", type=str, default='.',
                    help='chromosome (or part/contig) to keep (give a list of arguments separate the values with commas "X,Y") put \".\" for keep all chromosome (default: [2L,2R,3L,3R,^4_,X_])')
args = parser.parse_args()

name_bamfile = args.bam_file

min_size    = args.min_size
window      = args.window
chrom_list  = args.chrom.split(",")#KEEP ONLY

bamfile  = pysam.AlignmentFile(name_bamfile, "rb")

ID       = 0

#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False

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

            #if we have found INS to a good position
            if tupl[0] == 5 and tupl[1] >= min_size :
                
                # if ref_name_pred != REF :
                #     ref_name_pred = REF
                #     for indice, dic in enumerate(dico_clust) :
                #         print("\t".join([dic["REF"], str(dic["POS"]), "HARD." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"])]))

                #     dico_clust = []


                if seq and tupl[1] > min_size :

                    side = ""
                    if count_read_real - tupl[1] == 0 :
                        side = "L"
                    else :
                        side = "R"

                    # find = False
                    # if len(dico_clust) :
                    #     # for indice, dic in enumerate(dico_clust) :

                    #     #     if count_ref + reference_start <= dic["POS"] + window and count_ref + reference_start >= dic["POS"] - window :
                    #     #         find = True
                        
                    #     #         break
                    #     #     elif count_ref + reference_start <= dic["POS"] + window :
                    #     #         find = True
                    #     #         break

                    #     if not find :
                    #         if side == "L" :
                    #             dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real})
                    #         else :
                    #             dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1]})
                            
                    #         ID += 1
                    # else :
                    #     if side == "L" :
                    #         dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real})
                    #     else :
                    #         dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1]})
                        
                    #     ID += 1

                    if side == "L" :
                        dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq), "SEQ":seq})
                    else :
                        dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq), "SEQ":seq})
                    
                    ID += 1

print("\t".join(["#REF", "START", "END","ID", "ID_READ", "POS_REAL_READ", "FLAG", "SIDE", "SIZE_SEQ", "SEQ"]))
for indice, dic in enumerate(dico_clust) :
    print("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "HARD." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]))
# for indice, dic in enumerate(dico_clust) :
#     print("\t".join([dic["REF"], str(dic["POS"]), "HARD." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(dic["RS_LEFT"]), "RS_RIGHT=" + ",".join(dic["RS_RIGHT"])]))







