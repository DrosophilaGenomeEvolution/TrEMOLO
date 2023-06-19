import sys

fai_index = sys.argv[1]

file = open(fai_index, "r")

lines = [i.strip() for i in file.readlines()]


tab_var = []
for l in lines:
    chrom = l.split("\t")[0]
    size  = l.split("\t")[1]
    tab_var.append([chrom, int(size)])

print("var BioCircosGenome = " + str(tab_var))


