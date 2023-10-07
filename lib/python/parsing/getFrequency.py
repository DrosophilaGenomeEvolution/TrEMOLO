import pysam
import argparse
import numpy
import pandas as pd
from multiprocessing import Pool, cpu_count, Lock, Manager, current_process
import os

lock = Lock()


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



te_size_dic = {}
sv_size_dic = {}


df = pd.read_csv(filepath_or_buffer=args.position, sep="\t")
bamfile = pysam.AlignmentFile(name_bamfile, "rb")


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


for _, row in df.iterrows():
    qseqid = row["qseqid"]
    ID     = qseqid.split(":")[4]
    sseqid = row["sseqid"] #TE
    RS     = qseqid.split(":")[5]

    position = {"chrom": qseqid.split(":")[0], "start":qseqid.split(":")[2], "end": qseqid.split(":")[3]}

    tab_doublons    = []
    ref_name_pred   = ""

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
                            found = True

            distance = min(abs( (reference_start + count_ref) -  int(position["end"]) ), abs( (reference_start + count_ref) -  int(position["start"]) ))
            infos = " ".join([read_name, ID, str(sseqid), str(RS),  str(sv_size_dic[ID]), str(te_size_dic[sseqid]), str(reference_start + count_ref) , str(position["chrom"]+":"+position["start"]+"-"+position["end"]), qseqid, str(distance), str(reference_start), str(count_ref),  str(count_read), str(count_read_real)])
            if not found:
                args.output.write(" ".join(["E", str(infos), ".", "."]) + "\n")
                #print("E", infos, ".", ".")




# def cout_reads(row, tmp_file):
#     qseqid = row["qseqid"]
#     ID     = qseqid.split(":")[4]
#     sseqid = row["sseqid"] #TE
#     RS     = qseqid.split(":")[5]

#     position = {"chrom": qseqid.split(":")[0], "start":qseqid.split(":")[2], "end": qseqid.split(":")[3]}

#     tab_doublons    = []
#     ref_name_pred   = ""

#     itr = bamfile.fetch(str(position["chrom"]), max(int(position["start"])-window, 1), int(position["end"])+window, multiple_iterators = True)

#     for read in itr:
#         #for e, read in enumerate(fetch):

#         start_query     = read.query_alignment_start
#         reference_start = read.reference_start
#         seq             = read.seq
#         read_name       = read.query_name
#         REF             = read.reference_name #chrom
            
#         count_ref       = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
#         count_read      = 0
#         count_read_real = 0

#         if read_name :
#             found = False

#             for tupl in read.cigartuples:
                
#                 if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
#                     count_ref  += tupl[1]

#                 if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
#                     count_read += tupl[1]

#                 if tupl[0] in [0, 1, 7, 4, 5] : #Check M,I,=,S,H CIGAR for postion on reads real (HARD)
#                     count_read_real += tupl[1]

#                 distance = min(abs( (reference_start + count_ref) -  int(position["end"]) ), abs( (reference_start + count_ref) -  int(position["start"]) ))
#                 infos = " ".join([read_name, ID, str(sseqid), str(RS), str(sv_size_dic[ID]), str(te_size_dic[sseqid]), str(reference_start + count_ref), str(position["chrom"]+":"+position["start"]+"-"+position["end"]), qseqid, str(distance), str(reference_start), str(count_ref),  str(count_read), str(count_read_real)])

#                 if max(int(position["start"])-window, 1) <= reference_start + count_ref and reference_start + count_ref <= int(position["end"])+window:

#                     #if we have found INS to a good position
#                     if tupl[0] == 1 and tupl[1] >= ((size_percent/100) * int(te_size_dic[sseqid])) :
#                         size_seq = tupl[1]
#                         if seq :
#                             seq_vr = seq[count_read-tupl[1]:count_read]
#                             #print("I", infos, size_seq, seq_vr)
#                             tmp_file.write(" ".join(["I", str(infos), str(size_seq), str(seq_vr)]) + "\n")
#                         else:
#                             tmp_file.write(" ".join(["I", str(infos), str(size_seq), "."]) + "\n")
#                             #print("I", infos, size_seq, ".")
#                         found = True

#                     #if we have found HARD to a good position
#                     if tupl[0] == 5 and tupl[1] >= min_size :

#                         identify = ":".join([REF, str(count_ref + reference_start), str(read_name), str(count_read_real), "L", str(len(seq))])
#                         if seq and tupl[1] > min_size and identify not in tab_doublons :

#                             tab_doublons.append(identify)
#                             tmp_file.write(" ".join(["H", str(infos), ".", "."]) + "\n")
#                             #print("H", infos, ".", ".")

#                             found = True

#                     #if we have found SOFT to a good position
#                     if tupl[0] == 4 and tupl[1] >= min_size :
                        
#                         if ref_name_pred != REF :
#                             ref_name_pred = REF

#                         if seq and tupl[1] > min_size :

#                             find = False
#                             seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ SOFT in reads
#                             tmp_file.write(" ".join(["S", str(infos), str(tupl[1]), str(seq_vr)]) + "\n")
#                             #print("S", infos, tupl[1], seq_vr)

#                             found = True

#             distance = min(abs( (reference_start + count_ref) -  int(position["end"]) ), abs( (reference_start + count_ref) -  int(position["start"]) ))
#             infos = " ".join([read_name, ID, str(sseqid), str(RS),  str(sv_size_dic[ID]), str(te_size_dic[sseqid]), str(reference_start + count_ref) , str(position["chrom"]+":"+position["start"]+"-"+position["end"]), qseqid, str(distance), str(reference_start), str(count_ref),  str(count_read), str(count_read_real)])
#             if not found:
#                 tmp_file.write(" ".join(["E", str(infos), ".", "."]) + "\n")
                #print("E", infos, ".", ".")




# def process_rows(chunk, files_temp):
#     pid = current_process().pid
#     temp_filename = f"temp_{pid}.txt"

#     if not temp_filename in files_temp:
#         files_temp.append(temp_filename)

#     with open(temp_filename, 'w') as f:
#         for _, row in chunk.iterrows():
#             cout_reads(row, f)
    


# def main():
#     manager = Manager()
#     temp_files = manager.list()

#     num_processes = cpu_count()
#     num_processes = 5
#     chunks = numpy.array_split(df, num_processes)

#     with Pool(processes=num_processes) as pool:
#         pool.starmap(process_rows, [(chunk, temp_files) for chunk in chunks])


#     for temp_file in temp_files:
#         if temp_file:
#             with open(temp_file, 'r') as f:
#                 args.output.write(f.read())
#             os.remove(temp_file)

# if __name__ == "__main__":
#     main()
