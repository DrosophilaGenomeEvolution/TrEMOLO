import pysam
import argparse
import numpy
import pandas as pd


parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion, at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")
parser.add_argument("position", type=str, help="position of TE")


#OPTION

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
position = None

format chrom:start-end:size
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


if position != None or True :
    for e, read in enumerate(bamfile.fetch(str(position["chrom"]), max(int(position["start"])-window, 1), int(position["end"])+window) ):
    #for e, read in enumerate(bamfile.fetch()):
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

                #if tupl[0] in [8]:
                    