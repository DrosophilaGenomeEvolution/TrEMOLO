import sys
import os
import pandas
import re

name_fastq = sys.argv[1]
name_fasta = sys.argv[2]
file_out   = open(name_fasta, "w")
file       = open(name_fastq, "r")
line       = file.readline()

sequence          = ""
qual              = ""
while line:
    if line[0] == "@":
        sequence = ""
        qual     = ""
        file_out.write(">"+line[1:])
        line = file.readline()
        while line and line[0] != "+":
            sequence += line.replace("\n", "")
            line = file.readline()
    
    if line and line[0] == "+":
        file_out.write(sequence+"\n")
        line = file.readline()
        while line and len(qual) < len(sequence):
            qual += line.replace("\n", "")
            line = file.readline()

file_out.close()
file.close()