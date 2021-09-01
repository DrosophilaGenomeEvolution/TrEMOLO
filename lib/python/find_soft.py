import pysam
import os
import sys
import argparse

#TODO CHECK ALL METHODS
parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion (soft reads), at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")


#OPTION
parser.add_argument("-m", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=30,
                    help="minimum size of sequence.")
args = parser.parse_args()

name_bamfile = args.bam_file

min_size = args.min_size
window   = args.window

bamfile  = pysam.AlignmentFile(name_bamfile, "rb")

ID       = 0

ref_name_pred = ""
dico_clust    = []
for e, read in enumerate(bamfile.fetch()):
    start_query     = read.query_alignment_start
    reference_start = read.reference_start
    seq             = read.seq
    read_name       = read.query_name
    REF             = read.reference_name
    
    count_ref  = 0 #number of nucleotides on the ref before reaching the insertion site
    count_read = 0

    if read_name :

        for tupl in read.cigartuples :
            
            if tupl[0] in [0, 2, 7] : #Check M,D,= CIGAR for position on ref
                count_ref  += tupl[1]

            if tupl[0] in [0, 1, 7, 4] : #Check M,I,=,S CIGAR for postion on reads
                count_read += tupl[1]

            #if we have found INS to a good position
            if tupl[0] == 4 and tupl[1] >= min_size :
                
                if ref_name_pred != REF :
                    ref_name_pred = REF
                    for indice, dic in enumerate(dico_clust) :
                        print("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(dic["RS_LEFT"]), "RS_RIGHT=" + ",".join(dic["RS_RIGHT"])]))

                    dico_clust = []


                if seq and tupl[1] > min_size :

                    side = ""
                    if count_read - tupl[1] == 0 :
                        side = "L"
                    else :
                        side = "R"

                    find = False
                    if len(dico_clust) :
                        for indice, dic in enumerate(dico_clust) :

                            if count_ref + reference_start <= dic["POS"] + window and count_ref + reference_start >= dic["POS"] - window :
                                find = True
                                if side == "L":
                                    dico_clust[indice]["RS_LEFT"].append(str(read_name))
                                    
                                    if int(dico_clust[indice]["BEST_LEFT"][1]) < tupl[1] :
                                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ INS in reads
                                        dico_clust[indice]["BEST_LEFT"] = [str(read_name) , str(tupl[1]), seq_vr]
                                else :
                                    dico_clust[indice]["RS_RIGHT"].append(str(read_name))
                                    
                                    if int(dico_clust[indice]["BEST_RIGHT"][1]) < tupl[1] :
                                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ INS in reads
                                        dico_clust[indice]["BEST_RIGHT"] = [str(read_name) , str(tupl[1]), seq_vr]

                                break
                            elif count_ref + reference_start <= dic["POS"] + window :
                                find = True
                                break

                        if not find :
                            seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ INS in reads
                            if side == "L" :
                                dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" ,  str(0), "NONE"], "RS_LEFT":[str(read_name)], "RS_RIGHT":[]})
                            else :
                                dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name)]})
                            
                            ID += 1
                    else :
                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ INS in reads
                        if side == "L" :
                            dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" , str(0), "NONE"], "RS_LEFT":[str(read_name)], "RS_RIGHT":[]})
                        else :
                            dico_clust.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID, "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name)]})
                        
                        ID += 1


for indice, dic in enumerate(dico_clust) :
    print("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(dic["RS_LEFT"]), "RS_RIGHT=" + ",".join(dic["RS_RIGHT"])]))







