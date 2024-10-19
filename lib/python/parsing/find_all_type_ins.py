import pysam
import argparse
import numpy
import os
from multiprocessing import Pool, Event, Manager
import time


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
parser.add_argument("--output-seq-tsd", dest="output_seq_tsd", type=argparse.FileType('w'), default=None,
                    help="output seq for found TSD")

parser.add_argument("--no-clipped", dest="no_clipped", action='store_true',
                    help="d'ont get clipped reads")
parser.add_argument("--time-limit", dest="time_limit", type=int, default=0,
                    help="max time to process")

parser.add_argument("-t", "--threads", dest="threads", type=int, default=4,
                    help="Number of threads")

parser.add_argument("-f", "--flank-size", dest="flank_size", type=int, default=30,
                    help="flanking size sequence for get TSD.")
parser.add_argument("-w", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=30,
                    help="minimum size sequence.")

args = parser.parse_args()


def process_chrom(args):
    pid = os.getpid()
    try:
        i_chrom, chrom, stop_event = args
        print(f"Process PID {pid} : {chrom.ljust(31, ' ')} : {i_chrom} : begining...", flush=True)
        parse_chrom_get_sv(chrom, i_chrom, stop_event)
        print(f"Process PID {pid} : {chrom.ljust(31, ' ')} : {i_chrom} : STOP", flush=True)
    except Exception as e:
        print(f"Raise with process PID {pid} : {e}", flush=True)
        raise

def parse_chrom_get_sv(chrom, i_chrom, stop_event):

    output_soft_bis = open(f"{args.output_soft}.bis.tmp.{i_chrom}", "w")
    output_soft     = open(f"{args.output_soft}.tmp.{i_chrom}", "w") 
    output_hard     = open(f"{args.output_hard}.tmp.{i_chrom}", "w") 
    output_ins      = open(f"{args.output_ins}.tmp.{i_chrom}", "w")

    soft_line = []
    soft_bis_line = []
    hard_line = []

    ID_HARD         = 0
    ID_SOFT         = 0

    bamfile = pysam.AlignmentFile(args.bam_file.name, "rb")
    for e, read in enumerate(bamfile.fetch(chrom)):

        if stop_event.is_set():
            print("STOP:", i_chrom)
            break

        reference_start = read.reference_start
        seq             = read.seq
        read_name       = read.query_name
        REF             = chrom #chrom
        
        count_ref       = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
        count_read      = 0
        count_read_real = 0

        tab_doublons    = set()

        if read_name :

            for tupl in read.cigartuples:
                
                if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
                    count_ref  += tupl[1]

                if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
                    count_read += tupl[1]

                if tupl[0] in [0, 1, 7, 4, 5] : #Check M,I,=,S,H CIGAR for postion on reads real (HARD)
                    count_read_real += tupl[1]

                #if we have found INS to a good position
                if tupl[0] == 1 and tupl[1] >= args.min_size :
                    
                    if seq :
                        pos_min = max((count_read-tupl[1])-args.flank_size, 0)
                        pos_max = min((count_read+args.flank_size), len(seq))
                        # seq_te_tsd = seq[count_read-tupl[1]:count_read]
                        fk_L = seq[pos_min:pos_min+args.flank_size]
                        fk_R = seq[pos_max-args.flank_size:pos_max]
                        seq_full = seq[pos_min:pos_max]
                        #seq_vr = seq[count_read-tupl[1]:count_read]
                        seq_vr = seq_full
                        output_ins.write("\t".join([str(REF), str(reference_start + count_ref), str(reference_start + count_ref + 1), str(read_name), str(seq_vr), str(count_read), str(count_read_real), str(tupl[1])]) + "\n")
                        if args.output_seq_tsd != None:
                            args.output_seq_tsd.write("\t".join([str(REF), str(reference_start + count_ref), str(read_name), str(f'{fk_L}|{seq_vr}|{fk_R}'), str(seq_full), str(count_read), str(count_read_real), str(tupl[1])]) + "\n")

                #if we have found HARD to a good position
                if not args.no_clipped and tupl[0] == 5 and tupl[1] >= args.min_size :

                    identify = ":".join([REF, str(count_ref + reference_start), str(read_name), str(count_read_real), "L", str(len(seq))])
                    if not args.no_clipped and seq and tupl[1] > args.min_size and identify not in tab_doublons :

                        side = ""
                        if count_read_real - tupl[1] == 0 :
                            side = "L"
                        else :
                            side = "R"

                        if side == "L" :
                            hard_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_HARD}", "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq), "SEQ":seq})
                        else :
                            hard_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_HARD}", "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq), "SEQ":seq})
                        
                        tab_doublons.add(identify)

                        ID_HARD += 1

                #if we have found SOFT to a good position
                if not args.no_clipped and tupl[0] == 4 and tupl[1] >= args.min_size :
                    
                    if seq and tupl[1] > args.min_size :

                        side = ""
                        if count_read - tupl[1] == 0 :
                            side = "L"
                        else :
                            side = "R"

                        find = False
                        if len(soft_line) :
                            for indice, dic in enumerate(soft_line) :

                                if count_ref + reference_start <= dic["POS"] + args.window and count_ref + reference_start >= dic["POS"] - args.window :
                                    find = True
                                    if side == "L":
                                        soft_line[indice]["RS_LEFT"].append(str(read_name) + ":" + seq_vr)

                                        seq_vr_bis = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                        soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":soft_line[indice]["ID"], "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq_vr_bis), "SEQ":seq_vr_bis})
                                        
                                        if int(soft_line[indice]["BEST_LEFT"][1]) < tupl[1] :
                                            seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads TODO : Check
                                            soft_line[indice]["BEST_LEFT"] = [str(read_name) , str(tupl[1]), seq_vr]
                                    else :
                                        soft_line[indice]["RS_RIGHT"].append(str(read_name) + ":" + seq_vr)

                                        seq_vr_bis = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                        soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":soft_line[indice]["ID"], "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq_vr_bis), "SEQ":seq_vr_bis})
                                        
                                        if int(soft_line[indice]["BEST_RIGHT"][1]) < tupl[1] :
                                            seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads TODO : Check
                                            soft_line[indice]["BEST_RIGHT"] = [str(read_name) , str(tupl[1]), seq_vr]
                                        
                                    
                                    soft_line[indice]["NB_RS"] += 1

                                    break
                                elif count_ref + reference_start <= dic["POS"] + args.window :
                                    find = True
                                    break

                            if not find :
                                seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                if side == "L" :
                                    soft_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" ,  str(0), "NONE"], "RS_LEFT":[str(read_name) + ":" + str(seq_vr)], "RS_RIGHT":[], "NB_RS":1})
                                    soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq_vr), "SEQ":seq_vr})
                                else :
                                    soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq_vr), "SEQ":seq_vr}) 
                                    soft_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name) + ":" + str(seq_vr)], "NB_RS":1})
                                
                                ID_SOFT += 1
                        else :
                            seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                            if side == "L" :
                                soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "ID_READ":read_name, "POS_REAL_READ":count_read_real, "FLAG":read.flag, "SIDE":"L", "SIZE_SEQ":len(seq_vr), "SEQ":seq_vr})
                                soft_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "BEST_LEFT":[str(read_name) ,  str(tupl[1]), seq_vr],"BEST_RIGHT":["NONE" , str(0), "NONE"], "RS_LEFT":[str(read_name) + ":" + str(seq_vr)], "RS_RIGHT":[], "NB_RS":1})
                            else :
                                soft_bis_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "ID_READ":read_name, "POS_REAL_READ":count_read_real-tupl[1], "FLAG":read.flag, "SIDE":"R", "SIZE_SEQ":len(seq_vr), "SEQ":seq_vr}) 
                                soft_line.append({"REF":REF, "POS": count_ref + reference_start, "ID":f"{i_chrom}_{ID_SOFT}", "BEST_LEFT":["NONE" ,  str(0), "NONE"],"BEST_RIGHT":[str(read_name) , str(tupl[1]), seq_vr], "RS_LEFT":[], "RS_RIGHT":[str(read_name) + ":" + str(seq_vr)], "NB_RS":1})
                            
                            ID_SOFT += 1

                    if len(soft_line) > 5000: 
                        for indice, dic in enumerate(soft_line) :
                            output_soft.write("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(list(numpy.unique(dic["RS_LEFT"]))), "RS_RIGHT=" + ",".join(list(numpy.unique(dic["RS_RIGHT"]))), "NB_RS=" + str(dic["NB_RS"])]) + "\n")
                        
                    ## SOFT BIS
                    if len(soft_bis_line) > 5000: 
                        for indice, dic in enumerate(soft_bis_line) :
                            output_soft_bis.write("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "SOFT." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]) + "\n")

                    soft_line = []
                    soft_bis_line = []

    #HARD
    for indice, dic in enumerate(hard_line) :
        output_hard.write("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "HARD." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]) + "\n")

    #SOFT_BIS : LAST DATA 
    for indice, dic in enumerate(soft_bis_line) :
        output_soft_bis.write("\t".join([dic["REF"], str(dic["POS"]), str(dic["POS"]+1), "SOFT." + str(dic["ID"]), str(dic["ID_READ"]), str(dic["POS_REAL_READ"]), str(dic["FLAG"]), str(dic["SIDE"]), str(dic["SIZE_SEQ"]), str(dic["SEQ"])]) + "\n")

    #SOFT : LAST DATA 
    for indice, dic in enumerate(soft_line) :
        output_soft.write("\t".join([dic["REF"], str(dic["POS"]), "SOFT." + str(dic["ID"]), ";".join(["BEST_L_RS=" + dic["BEST_LEFT"][0], "BEST_L_SIZE=" + dic["BEST_LEFT"][1], "BEST_L_SEQ=" + dic["BEST_LEFT"][2]]), ";".join(["BEST_R_RS=" + dic["BEST_RIGHT"][0], "BEST_R_SIZE=" + dic["BEST_RIGHT"][1], "BEST_R_SEQ=" + dic["BEST_RIGHT"][2]]), "RS_LEFT=" + ",".join(dic["RS_LEFT"]), "RS_RIGHT=" + ",".join(dic["RS_RIGHT"]), "NB_RS=" + str(dic["NB_RS"])]) + "\n")



if __name__ == '__main__':
    import multiprocessing

    bamfile = pysam.AlignmentFile(args.bam_file.name, "rb")
    chromosomes = bamfile.references
    print("chromosomes:", chromosomes)
    print("--no-clipped", args.no_clipped)
    # bamfile.close()

    # event for send signal
    manager = Manager()
    stop_event = manager.Event()

    time_limit = args.time_limit * 3600
    print("Time limit:", time_limit)
    start_time = time.time()

    args_list = [(i_chr, chrom, stop_event) for i_chr, chrom in enumerate(chromosomes)]

    with Pool(processes=args.threads) as pool:
        results = pool.map_async(process_chrom, args_list)

        while not results.ready():
            elapsed_time = time.time() - start_time
            if time_limit > 0 and elapsed_time >= time_limit:
                print("Time limit collapsed, stop processus...")
                stop_event.set()  # send sign
                break
        pool.close()
        pool.join()

    # merge temporary files
    with open(f"{args.output_ins}", "w") as out:
        for i_chrom, chrom in enumerate(chromosomes):
            try:
                with open(f"{args.output_ins}.tmp.{i_chrom}", "r") as infile:
                    for read in infile:
                        out.write(read)
            except FileNotFoundError:
                continue

    if not args.no_clipped:
        with open(f"{args.output_soft}.bis", "w") as out:
            out.write("\t".join(["#REF", "START", "END","ID", "ID_READ", "POS_REAL_READ", "FLAG", "SIDE", "SIZE_SEQ", "SEQ"]) + "\n")
            for i_chrom, chrom in enumerate(chromosomes):
                try:
                    with open(f"{args.output_soft}.bis.tmp.{i_chrom}", "r") as infile:
                        for read in infile:
                            out.write(read)
                except FileNotFoundError:
                    continue
        
        with open(f"{args.output_soft}", "w") as out:
            for i_chrom, chrom in enumerate(chromosomes):
                try:
                    with open(f"{args.output_soft}.tmp.{i_chrom}", "r") as infile:
                        for read in infile:
                            out.write(read)
                except FileNotFoundError:
                    continue
        
        with open(f"{args.output_hard}", "w") as out:
            out.write("\t".join(["#REF", "START", "END","ID", "ID_READ", "POS_REAL_READ", "FLAG", "SIDE", "SIZE_SEQ", "SEQ"]) + "\n")
            for i_chrom, chrom in enumerate(chromosomes):
                try:
                    with open(f"{args.output_hard}.tmp.{i_chrom}", "r") as infile:
                        for read in infile:
                            out.write(read)
                except FileNotFoundError:
                    continue
        
