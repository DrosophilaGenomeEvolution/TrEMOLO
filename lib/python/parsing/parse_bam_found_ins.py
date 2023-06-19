import pysam
import os
import sys
import argparse


parser = argparse.ArgumentParser(description="parse bam file to get sequence of insertion, at specific position in bed file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("bam_file", metavar='<bam-file>', option_strings=['bam-file'], type=argparse.FileType('r'),
                    help="Input bam file")

parser.add_argument("bed_file", metavar='<bed-file>', type=argparse.FileType('r'),
                    help="")


#OPTION
parser.add_argument("-m", "--max_distance", dest='max_distance', type=int, default=30,
                    help="Maximum distance to group SV together.")
parser.add_argument("-s", "--min-size-percent", dest="min_size_percent", type=int, default=0.9,
                    help="minimum percentage of TE size.")
args = parser.parse_args()

name_bamfile = args.bam_file
name_bedfile = args.bed_file.name

min_size_percent = args.min_size_percent
max_distance     = args.max_distance

bamfile = pysam.AlignmentFile(name_bamfile, "rb")
befile  = open(name_bedfile, "r")


for line in befile:

    #format bed requier chom, start, end, name, size, seq
    chrom   = line.split("\t")[0]
    start   = int(line.split("\t")[1])
    end     = int(line.split("\t")[2])
    name    = line.split("\t")[3]
    size    = int(line.split("\t")[4])
    bed_seq = line.split("\t")[5]

    #print("NAME: ", name)
    number_read_support = 1
    for read in bamfile.fetch(chrom, (start-1), (start+1)):
        start_query     = read.query_alignment_start
        reference_start = read.reference_start
        seq             = read.seq
        
        count_ref  = 0 #nombre de nucleotide sur la ref avant d'atteindre le site d'insertion
        count_read = 0
        for tupl in read.cigartuples:
            
            if tupl[0] in [0, 2, 7]: #Check M,D,= CIGAR for position on ref
                count_ref  += tupl[1]

            if tupl[0] in [0, 1, 7, 4]: #Check M,I,=,S CIGAR for postion on reads
                count_read += tupl[1]

            #if we have found INS to a good position
            if tupl[0] == 1 and tupl[1] > (size*min_size_percent) and abs((reference_start + count_ref)-start) < max_distance:
                
                if seq :
                    seq_vr = seq[count_read-tupl[1]:count_read] #get SEQ INS in reads

                    print(">" + name + ":" + str(number_read_support))
                    print(seq_vr)
                    number_read_support += 1

    #Put the sequence report by sniffles
    if number_read_support == 1:
        print(">" + name + ":" + str(number_read_support))
        print(bed_seq)
    else :
        print(">" + name + ":" + str(0))
        print(bed_seq)










