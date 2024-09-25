#!/home/mourdas/anaconda3/bin/python3
import os
import sys

##FASTA CLASSIC TO FASTA BASIC

fasta = open(sys.argv[1], "r")

line = fasta.readline()
while line :
    if len(line) and  line.startswith(">"):
        head = line.strip()
        line = fasta.readline()
        sequence = ""
        while len(line) and not line.startswith(">") :
            sequence += line.strip()
            line = fasta.readline()
        print(head)
        print(sequence.upper())
    else:
        line = fasta.readline()

fasta.close()
