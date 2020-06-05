import sys
import os
import re
import argparse

parser = argparse.ArgumentParser()

#MAIN ARGS
parser.add_argument("vcf_file", type=str,
                    help="VCF file to parse")

#OPTION
parser.add_argument("-d", "--directory_name", type=str,
                    help="name directory contains regions files", default="REGION")

args = parser.parse_args()


rep_region = args.directory_name
os.system("rm -fr " + rep_region)
os.system("mkdir -p " + rep_region)

file = open(args.vcf_file, "r")

tab_chrom = ["2L", "2R", "3L", "3R", "X", "4"]

line = file.readline()
while line:
    if re.search("^[^#]", line) :
        spl    = line.split("\t")
        chrom  = spl[0].split("_")[0]
        ID     = spl[2]
        type_v = spl[4]
        start  = spl[1]
        end    = spl[7].split(";")[3].split("=")[1]
        rname  = spl[7].split(";")[9].split("=")[1]
        lrname = rname.split(",")

        if type_v[0] == "<" and type_v not in ["<DEL>"] and chrom in tab_chrom:
            file_out  = open("./" + rep_region + "/reads_"+str(chrom)+":"+str(ID)+":"+str(start)+"-"+str(end)+".txt", "w")
            for index, value in enumerate(lrname):
                file_out.write(str(value)+"\n")

            file_out.close()
    line = file.readline()

file.close()