import pysam
import os
import sys
import re
import argparse
import numpy

#TODO CHECK ALL METHODS
parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion (soft reads), at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")

parser.add_argument("bed_file", metavar='<bed-file>', type=argparse.FileType('r'),
                    help="")


#OPTION
parser.add_argument("-m", "--window", dest='window', type=int, default=50,
                    help="Maximum distance to group SV together.")

args = parser.parse_args()

name_bamfile = args.bam_file
window       = args.window

bamfile      = pysam.AlignmentFile(name_bamfile, "rb")

#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False


for line in befile:

    #format bed requier chom, start, end, name, size, seq
    chrom   = line.split("\t")[0]
    start   = int(line.split("\t")[1])
    end     = int(line.split("\t")[2])
    name    = line.split("\t")[3]
    size    = int(line.split("\t")[4])
    bed_seq = line.split("\t")[5]

    #print("NAME: ", name)
    ref_name_pred = ""
    for read in bamfile.fetch(chrom, (start-1), (start+1)):
        start_query     = read.query_alignment_start
        reference_start = read.reference_start
        seq             = read.seq
        read_name       = read.query_name
        REF             = read.reference_name #chrom

        count_ref  = 0 #number of nucleotides on the ref before reaching the insertion site
        count_read = 0

        if read_name :
            found = False
            for tupl in read.cigartuples :
                
                if tupl[0] in [0, 2, 7] : #Check M,D,= CIGAR for position on ref
                    count_ref  += tupl[1]

                if tupl[0] in [0, 1, 7, 4] : #Check M,I,=,S CIGAR for postion on reads
                    count_read += tupl[1]

                position_read_start = count_read-tupl[1]
                position_read_end   = count_read
                position_ref  = reference_start + count_ref
                #if we have found SOFT to a good position
                if tupl[0] == 4 && start - window < position_ref && position_ref < start + window :
                    #wirte
                    
                    found = True

            if not found:
                #write





