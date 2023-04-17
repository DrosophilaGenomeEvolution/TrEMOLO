#!/home/mourdas/anaconda3/bin/python3
import os
import sys

##FASTA CLASSIC TO FASTA BASIC

fasta = open(sys.argv[1], "r")

line = fasta.readline()
while line :
    if len(line) and  line[0] == ">":
        head = line.strip()
        line = fasta.readline()
        sequence = ""
        while len(line) and line[0] != ">" :
            sequence += line.strip()
            line = fasta.readline()
        print(head)
        print(sequence.upper())
    if len(line) and line[0] != ">":
        line = fasta.readline()
    elif len(line) == 0:
        line = fasta.readline()

fasta.close()
