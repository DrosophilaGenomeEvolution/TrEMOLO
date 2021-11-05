import pysam
import os
import sys
import re
import argparse

#TODO CHECK ALL METHODS
parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion (soft reads), at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")
parser.add_argument("position", type=str, help="position of TE")

#OPTION
parser.add_argument("-p", "--per-max-size-diff", dest='percent_size', type=float, default=0.8,
                    help="Maximum size difference.")
parser.add_argument("-d", "--min-dist", dest='min_dist', type=int, default=30,
                    help="Maximum distance of INSERTION to TE.")

args = parser.parse_args()

name_bamfile = args.bam_file
percent_size = args.percent_size
min_dist     = args.min_dist

#format chrom:start-end:size
if args.position != None:
    position = {"chrom": args.position.split(":")[0], "start": int(args.position.split(":")[1].split("-")[0]), "end": int(args.position.split(":")[1].split("-")[1]), "size": int(args.position.split(":")[2])}
    size_TE  = int(position["size"])
    if int(position["start"]) == int(position["end"]):
        print("Warning: start position equal end position, +1 to end")
        position["end"] += 1
    elif int(position["start"]) > int(position["end"]) :
        print("ERROR: start grather than end")
        exit(-1)
else:
    position = None

bamfile  = pysam.AlignmentFile(name_bamfile, "rb")

ens_rd  = set()
ens_ins = set()

doublons_rd = set()
total = []
if position != None :
    
    for e, read in enumerate(bamfile.fetch(str(position["chrom"]), int(position["start"]), int(position["end"])) ):
        start_query     = read.query_alignment_start
        reference_start = read.reference_start
        seq             = read.seq
        read_name       = read.query_name
        REF             = read.reference_name #chrom
        mapping_quality = read.mapping_quality

        count_ref  = 0 #number of nucleotides on the ref before reaching the insertion site
        count_read = 0

        
        
        if read_name and mapping_quality >= 30 :#TODO quality reajust
            if read_name in ens_rd :
                doublons_rd.add(read_name)

            ens_rd.add(read_name)

            total.append(read_name)
            for tupl in read.cigartuples :
                
                if tupl[0] in [0, 2, 7] : #Check M,D,= CIGAR for position on ref
                    count_ref  += tupl[1]

                if tupl[0] in [0, 1, 7, 4] : #Check M,I,=,S CIGAR for postion on reads
                    count_read += tupl[1]

                position_read_start = count_read - tupl[1]
                position_read_end   = count_read
                position_ref        = reference_start + count_ref

                #if we have found INS to a good position
                if tupl[0] == 1 and int(position["start"]) - min_dist <= position_ref and position_ref <= int(position["start"]) + min_dist and tupl[1] >= percent_size * int(position["size"]) :
                    ens_ins.add(read_name)

print("DEPTH:" + str(len(ens_rd)))
print("INS:" + str(len(ens_ins)))
print("DOUBLONS:", doublons_rd)
print("total:", total, "len:", len(total))
print("percent:", (len(ens_ins)/len(ens_rd)))

