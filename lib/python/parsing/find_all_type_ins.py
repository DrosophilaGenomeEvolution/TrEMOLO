import pysam
import argparse
import numpy


parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion, at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")


#OPTION
parser.add_argument("--output-soft", dest="output_soft", type=str, default="SOFT.txt",
                    help="output soft")
parser.add_argument("--output-ins", dest="output_ins", type=str, default="INS.txt",
                    help="output ins")
parser.add_argument("--output-hard", dest="output_hard", type=str, default="HARD.txt",
                    help="output hard")

parser.add_argument("-w", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")
parser.add_argument("-m", "--max_distance", dest='max_distance', type=int, default=30,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=200,
                    help="minimum size sequence.")

args = parser.parse_args()

name_bamfile = args.bam_file

min_size     = args.min_size
max_distance = args.max_distance
window       = args.window

bamfile = pysam.AlignmentFile(name_bamfile, "rb")

ID_HARD         = 0
ID_SOFT         = 0
tab_doublons    = []
ref_name_pred   = ""

output_soft = open(args.output_soft, "w") 
output_hard = open(args.output_hard, "w") 
output_ins  = open(args.output_ins, "w")

dico_clust_ins  = []
dico_clust_soft = []
dico_clust_hard = []
for e, read in enumerate(bamfile.fetch()):
    start_query     = read.query_alignment_start
    reference_start = read.reference_start
    seq             = read.seq
    read_name       = read.query_name
    REF             = read.reference_name #chrom
        
    count_ref       = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
    count_read      = 0
    count_read_real = 0

    if read_name :

        for tupl in read.cigartuples:
            
            if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
                count_ref  += tupl[1]

            if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
                count_read += tupl[1]

            if tupl[0] in [0, 1, 7, 4, 5] : #Check M,I,=,S,H CIGAR for postion on reads real (HARD)
                count_read_real += tupl[1]

            #if we have found INS to a good position
            if tupl[0] == 1 and tupl[1] >= min_size :
                
                if seq :
                    seq_vr = seq[count_read-tupl[1]:count_read]
                    output_ins.write("\t".join([str(REF), str(reference_start + count_ref), str(reference_start + count_ref + 1), str(read_name), str(seq_vr)]) + "\n")

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
                        dico_clust_hard.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_HARD, "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq), "SEQ":seq})
                    else :
                        dico_clust_hard.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_HARD, "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq), "SEQ":seq})
                    
                    tab_doublons.append(identify)

                    ID_HARD += 1

            #if we have found SOFT to a good position
            if tupl[0] == 4 and tupl[1] >= min_size :
                
                if ref_name_pred != REF :
                    ref_name_pred = REF
                    for indice, dic in enumerate(dico_clust_soft) :
                        output_soft.write("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(list(numpy.unique(dic["RS_LEFT"]))), "RS_RIGHT=" + ",".join(list(numpy.unique(dic["RS_RIGHT"])))]) + "\n")
                    
                    dico_clust_soft = []


                if seq and tupl[1] > min_size :

                    side = ""
                    if count_read - tupl[1] == 0 :
                        side = "L"
                    else :
                        side = "R"

                    find = False
                    if len(dico_clust_soft) :
                        for indice, dic in enumerate(dico_clust_soft) :

                            if count_ref + reference_start <= dic["POS"] + window and count_ref + reference_start >= dic["POS"] - window :
                                find = True
                                if side == "L":
                                    dico_clust_soft[indice]["RS_LEFT"].append(str(read_name) + ":" + seq_vr)
                                    
                                    if int(dico_clust_soft[indice]["BEST_LEFT"][1]) < tupl[1] :
                                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                        dico_clust_soft[indice]["BEST_LEFT"] = [str(read_name) , str(tupl[1]), seq_vr]
                                else :
                                    dico_clust_soft[indice]["RS_RIGHT"].append(str(read_name) + ":" + seq_vr)
                                    
                                    if int(dico_clust_soft[indice]["BEST_RIGHT"][1]) < tupl[1] :
                                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                        dico_clust_soft[indice]["BEST_RIGHT"] = [str(read_name) , str(tupl[1]), seq_vr]
                                
                                dico_clust_soft[indice]["NB_RS"] += 1

                                break
                            elif count_ref + reference_start <= dic["POS"] + window :
                                find = True
                                break

                        if not find :
                            seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                            if side == "L" :
                                dico_clust_soft.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_SOFT, "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" ,  str(0), "NONE"], "RS_LEFT":[str(read_name) + ":" + str(seq_vr)], "RS_RIGHT":[], "NB_RS":1})
                            else :
                                dico_clust_soft.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_SOFT, "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name) + ":" + str(seq_vr)], "NB_RS":1})
                            
                            ID_SOFT += 1
                    else :
                        seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                        if side == "L" :
                            dico_clust_soft.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_SOFT, "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" , str(0), "NONE"], "RS_LEFT":[str(read_name) + ":" + str(seq_vr)], "RS_RIGHT":[], "NB_RS":1})
                        else :
                            dico_clust_soft.append({"REF":REF, "POS": count_ref + reference_start, "ID":ID_SOFT, "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name) + ":" + str(seq_vr)], "NB_RS":1})
                        
                        ID_SOFT += 1


#INS


#HARD
output_hard.write("\t".join(["#REF", "START", "END","ID", "ID_READ", "POS_REAL_READ", "FLAG", "SIDE", "SIZE_SEQ", "SEQ"]) + "\n")
for indice, dic in enumerate(dico_clust_hard) :
    output_hard.write("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "HARD." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]) + "\n")

#SOFT
for indice, dic in enumerate(dico_clust_soft) :
    output_soft.write("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(dic["RS_LEFT"]), "RS_RIGHT=" + ",".join(dic["RS_RIGHT"]), "NB_RS=" + str(dic["NB_RS"])]) + "\n")



