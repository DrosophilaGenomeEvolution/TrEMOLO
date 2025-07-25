#!/usr/bin python3
# -*- coding: utf-8 -*-

###################################################################################################################################
#
# Copyright 2019-2020 IRD-CNRS-Lyon1 University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Mourdas MOHAMED
#                                                                                                                                   
####################################################################################################################################

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
                            [2L,2R,3L,3R,^4_,X_]
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
parser.add_argument("-c", "--chrom", type=str, default='2L,2R,3L,3R,^4_,X_',
                    help='chromosome (or part) to keep (give a list of arguments separate the values ​​with commas "X,Y") [2L,2R,3L,3R,^4_,X_]')
#Warning : not used now
parser.add_argument("-m", "--min_len_seq", type=int, default=1000,
                    help="minimum size of the sequence to keep (default: [1000])")
#Warning : not used now
parser.add_argument("--max-len-seq", dest='max_len_seq', type=int, default=-1,
                    help="maximum size of the sequence to keep (default: [-1])")

args = parser.parse_args()


repRegion = args.directory_name
os.system("rm -fr " + repRegion)
os.system("mkdir -p " + repRegion)

#check if an element of the array is at least part of the value
def regexInList(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False


listId = None
if args.id_sv:
    fileId = open(args.id_sv, "r")
    listId = [i.strip() for i in fileId.readlines()]

file = open(args.vcf_file, "r")

#regexInList(chrom, chromList)
typeList   = args.type.split(",")#NOT KEEP THIS
chromList  = args.chrom.split(",")#KEEP ON

print("[" + str(sys.argv[0]) + "]", "exclude type : ", typeList)
print("[" + str(sys.argv[0]) + "]", "chromosome liste : ", chromList)

dicoChrom = {}#just for resume

line = file.readline()

versionVcf =  line.split("=")[1].strip()
print("versionVcf=", versionVcf)
##fileformat=VCFv4.3
#Check format vcf file
if versionVcf != "VCFv4.3" and versionVcf != "VCFv4.2" and versionVcf != "VCFv4.1":
    print("[" + str(sys.argv[0]) + "] : ERROR format vcf must be VCFv4.3 not " + str(line.split("=")[1].strip()) )
    print("[" + str(sys.argv[0]) + "] : please change format of vcf file, or use Snifflesv1.0.10 for genrate the good vcf file")
    exit(1)



while line:
    if re.search("^[^#]", line) :
        spl    = line.split("\t")
        chrom  = spl[0]
        start  = spl[1]
        ID     = spl[2]
        typeV = spl[4]
        
        #end    = spl[7].split(";")[3].split("=")[1]
        expEnd = re.search(r'END=([^;]*);', line.strip())

        
        if versionVcf == "VCFv4.3" :
            ID = "sniffles." + typeV.replace("<", "").replace(">", "") + "." + ID#For False ID sniffles
            expRName = re.search(r'RNAMES=([^;]*);', line.strip())
        
        elif versionVcf == "VCFv4.2" :#svim
            expRName = re.search(r'READS=([^;\t]*)', line.strip())
        
        elif versionVcf == "VCFv4.1" :
            expRName = re.search(r'RNAMES=([^;\t]*)', line.strip())
            seq = spl[4]
            typeV  = "<" + re.search("SVTYPE=([A-Z]+)", spl[7]).group(1) + ">" if re.search("SVTYPE=([A-Z]+)", spl[7]) else "None"
            ID = "sniffles." + typeV.replace("<", "").replace(">", "") + "." + ID#For False ID sniffles


        if listId == None or (ID in listId) :

            if typeV[0] == "<" and typeV not in typeList and regexInList(chrom, chromList) and expEnd != None and expRName != None :
                end   = expEnd.group(1)
                rname = expRName.group(1)

                lrname = rname.split(",")
                #just for resume
                if chrom not in dicoChrom:
                    dicoChrom[chrom] = 1
                else:
                    dicoChrom[chrom] += 1

                fileOut  = open( repRegion + "/reads_"+str(chrom)+":"+str(ID)+":"+str(start)+"-"+str(end)+".txt", "w")
                for index, value in enumerate(lrname) :
                    fileOut.write(str(value) + "\n")

                fileOut.close()
    line = file.readline()

file.close()

print("[" + str(sys.argv[0]) + "] resume TE number by chromosome : ", dicoChrom)