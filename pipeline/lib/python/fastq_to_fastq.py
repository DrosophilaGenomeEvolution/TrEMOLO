import sys
import os
import re

name_fastq = sys.argv[1]
output     = sys.argv[2]

print("[" + sys.argv[0] + " INFO]", name_fastq, output)

file       = open(name_fastq, "r")
file_out   = open(output, "w")
line       = file.readline()

sequence          = ""
qual              = ""
while line:
    if line[0] == "@":
        sequence = ""
        qual     = ""
        file_out.write("@"+line[1:])
        line = file.readline()
        while line and line[0] != "+":
            sequence += line.replace("\n", "")
            line = file.readline()
    
    if line and line[0] == "+":
        file_out.write(sequence+"\n")
        file_out.write("+\n")
        line = file.readline()
        while line and len(qual) < len(sequence):
            qual += line.replace("\n", "")
            line = file.readline()
        file_out.write(qual+"\n")

file_out.close()
file.close()






