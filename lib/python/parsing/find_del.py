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
parser.add_argument("-m", "--max-size-diff", dest='window', type=int, default=20,
                    help="Maximum size difference.")
parser.add_argument("-d", "--min-dist", dest='min_dist', type=int, default=30,
                    help="minimum distance of DELETION to TE.")

args = parser.parse_args()

name_bamfile = args.bam_file
window       = args.window
min_dist     = args.min_dist


#format chrom:start-end
if args.position != None:
    position = {"chrom": args.position.split(":")[0], "start": args.position.split(":")[1].split("-")[0], "end": args.position.split(":")[1].split("-")[1]}
    size_TE  = int(position["end"]) - int(position["start"])
else:
    position = None

bamfile  = pysam.AlignmentFile(name_bamfile, "rb")

NB_DEL        = 0
read_names    = []

if position != None :
    for e, read in enumerate(bamfile.fetch(str(position["chrom"]), int(position["start"]), int(position["end"])) ):
        start_query     = read.query_alignment_start
        reference_start = read.reference_start
        seq             = read.seq
        read_name       = read.query_name
        REF             = read.reference_name #chrom

        count_ref  = 0 #number of nucleotides on the ref before reaching the insertion site
        count_read = 0

        if read_name :

            for tupl in read.cigartuples :
                
                if tupl[0] in [0, 2, 7] : #Check M,D,= CIGAR for position on ref
                    count_ref  += tupl[1]

                if tupl[0] in [0, 1, 7, 4] : #Check M,I,=,S CIGAR for postion on reads
                    count_read += tupl[1]

                #if we have found DEL to a good position
                if tupl[0] == 2 and size_TE-window <= tupl[1] and tupl[1] <= size_TE+window and read_name not in read_names and ( abs((count_ref + reference_start)-int(position["start"])) < min_dist or abs((count_ref + reference_start)-int(position["end"])) < min_dist )  :
                    read_names.append(read_name)
                    NB_DEL += 1

print(NB_DEL)



