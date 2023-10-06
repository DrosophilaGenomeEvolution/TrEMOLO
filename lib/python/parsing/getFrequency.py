import pysam
import argparse
import numpy
import pandas as pd


parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion, at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")
parser.add_argument("position", metavar='<position-file>', type=argparse.FileType('r'), help="position of TE")
parser.add_argument("te_size", metavar='<te-size-file>', type=argparse.FileType('r'), help="TE size")
parser.add_argument("sv_size", metavar='<sv-size-file>', type=argparse.FileType('r'), help="SV size")
parser.add_argument("output", metavar='<output-file>', type=argparse.FileType('w'), help="output file")


#OPTION
parser.add_argument("-w", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")
parser.add_argument("-m", "--max_distance", dest='max_distance', type=int, default=30,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=200,
                    help="minimum size sequence.")
parser.add_argument("-p", "--per-size", dest="per_size", type=int, default=80,
                    help="percent size sequence.")

args = parser.parse_args()

name_bamfile = args.bam_file

min_size     = args.min_size
max_distance = args.max_distance
window       = args.window
size_percent = args.per_size

bamfile = pysam.AlignmentFile(name_bamfile, "rb")

ID_HARD         = 0
ID_SOFT         = 0
tab_doublons    = []
ref_name_pred   = ""
position = None
te_size_dic = {}
sv_size_dic = {}



df = pd.read_csv(filepath_or_buffer=args.position, sep="\t")

#TE SIZE
file_te_size = args.te_size
for line in file_te_size.readlines():
    te_size_dic[line.strip().split("\t")[0]] = line.strip().split("\t")[1]

#SV SIZE
file_sv_size = args.sv_size
for line in file_sv_size.readlines():
    sv_size_dic[line.strip().split("\t")[0].split(":")[4]] = line.strip().split("\t")[1]

#sseqid  qseqid  pident  size_per        size_el mismatch        gapopen qstart  qend    sstart  send    evalue  bitscore
#df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


if position != None or True :
    for index, row in enumerate(df.values):
        qseqid = df["qseqid"].values[index]
        ID     = qseqid.split(":")[4]
        sseqid = df["sseqid"].values[index] #TE
        RS     = qseqid.split(":")[5]

        position = {"chrom": qseqid.split(":")[0], "start":qseqid.split(":")[2], "end": qseqid.split(":")[3]}

        for e, read in enumerate(bamfile.fetch(str(position["chrom"]), max(int(position["start"])-window, 1), int(position["end"])+window) ):
            start_query     = read.query_alignment_start
            reference_start = read.reference_start
            seq             = read.seq
            read_name       = read.query_name
            REF             = read.reference_name #chrom
                
            count_ref       = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
            count_read      = 0
            count_read_real = 0

            if read_name :
                found = False

                for tupl in read.cigartuples:
                    
                    if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
                        count_ref  += tupl[1]

                    if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
                        count_read += tupl[1]

                    if tupl[0] in [0, 1, 7, 4, 5] : #Check M,I,=,S,H CIGAR for postion on reads real (HARD)
                        count_read_real += tupl[1]

                    distance = min(abs( (reference_start + count_ref) -  int(position["end"]) ), abs( (reference_start + count_ref) -  int(position["start"]) ))
                    infos = " ".join([read_name, ID, str(sseqid), str(RS), str(sv_size_dic[ID]), str(te_size_dic[sseqid]), str(reference_start + count_ref), str(position["chrom"]+":"+position["start"]+"-"+position["end"]), qseqid, str(distance), str(reference_start), str(count_ref),  str(count_read), str(count_read_real)])

                    if max(int(position["start"])-window, 1) <= reference_start + count_ref and reference_start + count_ref <= int(position["end"])+window:

                        #if we have found INS to a good position
                        if tupl[0] == 1 and tupl[1] >= ((size_percent/100) * int(te_size_dic[sseqid])) :
                            size_seq = tupl[1]
                            if seq :
                                seq_vr = seq[count_read-tupl[1]:count_read]
                                #print("I", infos, size_seq, seq_vr)
                                args.output.write(" ".join(["I", str(infos), str(size_seq), str(seq_vr)]) + "\n")
                            else:
                                args.output.write(" ".join(["I", str(infos), str(size_seq), "."]) + "\n")
                                #print("I", infos, size_seq, ".")
                            found = True

                        #if we have found HARD to a good position
                        if tupl[0] == 5 and tupl[1] >= min_size :

                            identify = ":".join([REF, str(count_ref + reference_start), str(read_name), str(count_read_real), "L", str(len(seq))])
                            if seq and tupl[1] > min_size and identify not in tab_doublons :

                                tab_doublons.append(identify)
                                args.output.write(" ".join(["H", str(infos), ".", "."]) + "\n")
                                #print("H", infos, ".", ".")
                                ID_HARD += 1
                                found = True

                        #if we have found SOFT to a good position
                        if tupl[0] == 4 and tupl[1] >= min_size :
                            
                            if ref_name_pred != REF :
                                ref_name_pred = REF

                            if seq and tupl[1] > min_size :

                                find = False
                                seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
                                args.output.write(" ".join(["S", str(infos), str(tupl[1]), str(seq_vr)]) + "\n")
                                #print("S", infos, tupl[1], seq_vr)
                                ID_SOFT += 1
                                found = True

                distance = min(abs( (reference_start + count_ref) -  int(position["end"]) ), abs( (reference_start + count_ref) -  int(position["start"]) ))
                infos = " ".join([read_name, ID, str(sseqid), str(RS),  str(sv_size_dic[ID]), str(te_size_dic[sseqid]), str(reference_start + count_ref) , str(position["chrom"]+":"+position["start"]+"-"+position["end"]), qseqid, str(distance), str(reference_start), str(count_ref),  str(count_read), str(count_read_real)])
                if not found:
                    args.output.write(" ".join(["E", str(infos), ".", "."]) + "\n")
                    #print("E", infos, ".", ".")


