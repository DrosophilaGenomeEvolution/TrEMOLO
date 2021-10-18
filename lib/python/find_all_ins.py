import pysam
import os
import sys
import argparse


parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion, at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")


#OPTION
parser.add_argument("-m", "--max_distance", dest='max_distance', type=int, default=30,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size", dest="min_size", type=int, default=1000,
                    help="minimum size sequence.")
args = parser.parse_args()

name_bamfile = args.bam_file

min_size     = args.min_size
max_distance = args.max_distance

bamfile = pysam.AlignmentFile(name_bamfile, "rb")


for e, read in enumerate(bamfile.fetch()):
    start_query     = read.query_alignment_start
    reference_start = read.reference_start
    seq             = read.seq
    read_name       = read.query_name
    REF             = read.reference_name #chrom
        
    count_ref  = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
    count_read = 0

    if read_name :

        for tupl in read.cigartuples:
            
            if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
                count_ref  += tupl[1]

            if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
                count_read += tupl[1]

            #if we have found INS to a good position
            if tupl[0] == 1 and tupl[1] >= min_size :
                
                if seq :
                    seq_vr = seq[count_read-tupl[1]:count_read]
                    print("\t".join([str(REF), str(reference_start + count_ref), str(reference_start + count_ref + 1), str(read_name), str(seq_vr)]))











