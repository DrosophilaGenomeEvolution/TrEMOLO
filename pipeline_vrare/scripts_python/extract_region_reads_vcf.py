"""
    GET REGION IN VCF
    ==========================
    :author:  Mourdas MOHAMED 
    :contact: mourdas.mohamed@igh.cnrs.fr
    :date: 01/06/2020
    :version: 0.1
    Script description
    ------------------
    extract_region_reads_vcf.py get the ID of reads support by SV in vcf file
    -------
    >>> extract_region_reads_vcf.py file.vcf -d REGION_DIRECTORY
    
    Help Programm
    -------------
    usage: extract_region_reads_vcf.py [-h] [-d DIRECTORY_NAME] [-i ID_SV]
                                   [-t TYPE] [-c CHROM]
                                   vcf_file

    positional arguments:
      vcf_file              VCF file to parse

    optional arguments:
      -h, --help            show this help message and exit
      -d DIRECTORY_NAME, --directory_name DIRECTORY_NAME
                            name directory contains regions files
      -i ID_SV, --id_sv ID_SV
                            file id strucural variant to keep
      -t TYPE, --type TYPE  not keep this type on vcf (give a list of arguments
                            separate the values ​​with commas "<DEL>,<INS>")
                            ["<DEL>"]
      -c CHROM, --chrom CHROM
                            chromosome (or part) to keep (give a list of arguments
                            separate the values ​​with commas "X,Y")
                            [2L,2R,3L,3R,4,X]
"""


#Import
import sys
import os
import re
import argparse

#parse args
parser = argparse.ArgumentParser(description="Get the ID of reads support by SV in vcf file")

#MAIN ARGS
parser.add_argument("vcf_file", type=str,
                    help="VCF file to parse")

#OPTION
parser.add_argument("-d", "--directory_name", type=str,
                    help="name directory contains regions files", default="REGION")
parser.add_argument("-i", "--id_sv", type=str,
                    help="file id strucural variant to keep")
parser.add_argument("-t", "--type", type=str, default='<DEL>',
                    help="not keep this type on vcf (give a list of arguments separate the values ​​with commas \"<DEL>,<INS>\") [\"<DEL>\"]")
parser.add_argument("-c", "--chrom", type=str, default='2L,2R,3L,3R,4,X',
                    help='chromosome (or part) to keep (give a list of arguments separate the values ​​with commas "X,Y") [2L,2R,3L,3R,4,X]')


args = parser.parse_args()


rep_region = args.directory_name
os.system("rm -fr " + rep_region)
os.system("mkdir -p " + rep_region)

#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False


list_id = None
if args.id_sv:
    file_id = open(args.id_sv, "r")
    list_id = [i.strip() for i in file_id.readlines()]

file = open(args.vcf_file, "r")

#regex_in_list(chrom, chrom_list)
type_list   = args.type.split(",")#NOT KEEP THIS
chrom_list  = args.chrom.split(",")#KEEP ON

print("exclude type : ", type_list)
print("chromosome liste : ", chrom_list)

dico_chrom = {}#just for resume

line = file.readline()
while line:
    if re.search("^[^#]", line) :
        spl    = line.split("\t")
        chrom  = spl[0]
        ID     = spl[2]
        type_v = spl[4]
        start  = spl[1]
        end    = spl[7].split(";")[3].split("=")[1]
        rname  = spl[7].split(";")[9].split("=")[1]
        lrname = rname.split(",")

        #print(ID, list_id[0])
        if list_id == None or (list_id and ID in list_id) :

            if type_v[0] == "<" and type_v not in type_list and regex_in_list(chrom, chrom_list) :
                #just for resume
                if chrom not in dico_chrom:
                    dico_chrom[chrom] = 1
                else:
                    dico_chrom[chrom] += 1

                file_out  = open("./" + rep_region + "/reads_"+str(chrom)+":"+str(ID)+":"+str(start)+"-"+str(end)+".txt", "w")
                for index, value in enumerate(lrname) :
                    file_out.write(str(value) + "\n")

                file_out.close()
    line = file.readline()

file.close()

print("resume number by chrom : ", dico_chrom)