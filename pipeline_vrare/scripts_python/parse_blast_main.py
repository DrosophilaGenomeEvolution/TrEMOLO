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
# Version 0.1 written by Mourdas Mohamed
#                                                                                                                                   
####################################################################################################################################


"""
    GET TE CANDIDATE
    ==========================
    :author:  Mourdas MOHAMED 
    :contact: mourdas.mohamed@igh.cnrs.fr
    :date: 01/06/2020
    :version: 0.1
    Script description
    ------------------
    parse_blast_main.py filters a blast file in output format 6 to keep the candidate TE active
    -------
    >>> parse_blast_main.py input.bln output_ALL_TE.csv
    
    Help Programm
    -------------
    usage: parse_blast_main.py [-h] [-t NAME_FILE_TE] [-p MIN_PIDENT]
                           [-s MIN_SIZE_PERCENT] [-r MIN_READ_SUPPORT]
                           blast_file output_file

    positional arguments:
      blast_file            blast file format outfmt 6
      output_file           name of tabular file out

    optional arguments:
      -h, --help            show this help message and exit
      -t NAME_FILE_TE, --name_file_te NAME_FILE_TE
                            file TE for get size
      -p MIN_PIDENT, --min_pident MIN_PIDENT
                            minimum percend of identity [94]
      -s MIN_SIZE_PERCENT, --min_size_percent MIN_SIZE_PERCENT
                            minimum precent size of TE [90]
      -r MIN_READ_SUPPORT, --min_read_support MIN_READ_SUPPORT
                        minimum read support number [1]
"""


import numpy
import pandas as pd
import re
import os
import sys
import argparse

# import tkinter
# import pylab as plt
# import warnings
# warnings.filterwarnings("ignore")


parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active")

#MAIN ARGS
parser.add_argument("blast_file", type=str,
                    help="input blast file format outfmt 6")
parser.add_argument("output_file", type=str,
                    help="name of tabular output file ")

#OPTION
parser.add_argument("-t", "--name_file_te", type=str, default=None,
                    help="file TE for get size")
parser.add_argument("-p", "--min_pident", type=int, default=94,
                    help="minimum percentage of identity [94]")
parser.add_argument("-s", "--min_size_percent", type=int, default=90,
                    help="minimum percentage of TE size [90]")
parser.add_argument("-r", "--min_read_support", type=int, default=1,
                    help="minimum read support number [1]")
parser.add_argument("-k", "--type_sv_keep", type=str, default="INS",
                    help="type_sv_keep")

args = parser.parse_args()


name_file        = args.blast_file
min_pident       = args.min_pident
min_size_percent = args.min_size_percent
min_read_support = args.min_read_support
type_sv_keep     = args.type_sv_keep

output           = args.output_file
name_out         = output.split("/")[-1].split(".")[0]
dir_out          = "/".join(output.split("/")[:-1])

print("------------------------")
print("name_out :", name_out)


if dir_out != "":
    dir_out += "/"
else :
    dir_out += "./"


size_et = {}

if args.name_file_te != None:
    #GET_SIZE TE Database format fasta
    file  = open(args.name_file_te, "r")
    lines = file.readlines()

    for i, l in enumerate(lines):
        if l[0] == ">":
            size_et[l[1:].strip()] = len(lines[i + 1].strip())

    file.close()

else :#default dm6 canonnical TE
    size_et = {'Idefix': 7411, '17.6': 7439, '1731': 4648, '297': 6995, '3S18': 6126, '412': 7567, 'aurora-element': 4263, 'Burdock': 6411, 'copia': 5143, 'gypsy': 7469, 'mdg1': 7480, 'mdg3': 5519, 'micropia': 5461, 'springer': 7546, 'Tirant': 8526, 'flea': 5034, 'opus': 7521, 'roo': 9092, 'blood': 7410, 'ZAM': 8435, 'GATE': 8507, 'Transpac': 5249, 'Circe': 7450, 'Quasimodo': 7387, 'HMS-Beagle': 7062, 'diver': 6112, 'Tabor': 7345, 'Stalker': 7256, 'gtwin': 7411, 'gypsy2': 6841, 'accord': 7404, 'gypsy3': 6973, 'invader1': 4032, 'invader2': 5124, 'invader3': 5484, 'gypsy4': 6852, 'invader4': 3105, 'gypsy5': 7369, 'gypsy6': 7826, 'invader5': 4038, 'diver2': 4917, 'Dm88': 4558, 'frogger': 2483, 'rover': 7318, 'Tom1': 410, 'rooA': 7621, 'accord2': 7650, 'McClintock': 6450, 'Stalker4': 7359, 'Stalker2': 7672, 'Max-element': 8556, 'gypsy7': 5486, 'gypsy8': 4955, 'gypsy9': 5349, 'gypsy10': 6006, 'gypsy11': 4428, 'gypsy12': 10218, 'invader6': 4885, 'Helena': 1317, 'HMS-Beagle2': 7220, 'Osvaldo': 1543}
    #print("SIZE TE DEFAULT", list(size_et.keys())[:5], "...")

#GET BLAST OUTFMT 

df = pd.read_csv(name_file, "\t", header=None)

#get prefix file
name_file_withou_ext = name_file.split("/")[-1].split(".")[0] + "_"
print("file name prefix :", name_file_withou_ext)

df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

#KEEP ONLY TE on the LIST (GET_SIZE CANNONICAL TE)
df = df[df["sseqid"].isin(size_et.keys())]
print("keep only TE on list, shape :", df.shape)

#CALCUL SIZE AND PERCENT SIZE OF TE
tab_percent = []
tab_size    = []
for index, row in enumerate(df[["sseqid", "qend", "qstart"]].values):
    size_element = abs(int(row[1])-int(row[2]))#qend - qstart
    tab_percent.append(round((size_element/size_et[row[0]]) * 100, 1))
    tab_size.append(size_element)
    
#keeps the TE according to the percentage of identity and the percentage of size
df["size_per"] = tab_percent
df["size_el"]  = tab_size
df = df[df["size_per"] >= min_size_percent]
df = df[df["pident"] >= min_pident]
df = df.sort_values(by=["sseqid"])
df = df[["sseqid", "qseqid", "pident", "size_per", "size_el", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]]


print("shape for min_size_percent>=" + str(min_size_percent) + ", min_pident>=" + str(min_pident), df.shape)

#keep the TE with the highest score
best_score_match = []
best_score_match_index = []
chaine = ""
for index, row in enumerate(df.values):
    
    chaine = df["qseqid"].values[index] + df["sseqid"].values[index]
    if chaine not in best_score_match:
        sstart       = df["sstart"].values[index]
        send         = df["send"].values[index]
        if send < sstart:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "-"
        else:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "+"

        qseqid       = df["qseqid"].values[index]
        info_qseqid  = qseqid.split(":")
        read_support = int(info_qseqid[5])

        df_tmp        = df[df["qseqid"] == qseqid]
        maxe_bitscore = max(df_tmp["bitscore"].values)
        df_best_score = df_tmp[df_tmp["bitscore"] == maxe_bitscore]
        
        chaine = df_best_score["qseqid"].values[0] + df_best_score["sseqid"].values[0]
        best_score_match.append(chaine)
        
        #only type SV choice 
        find_type_SV = re.search(type_sv_keep, df_best_score["qseqid"].values[0])
        if find_type_SV and read_support >= min_read_support:
            best_score_match_index.append(index)
    
df = df.iloc[best_score_match_index]


print("shape best score match insertion :", df.shape)

#(optional)
indexs = []
for index, value in enumerate(df.values):
    qseqid  = df["qseqid"].values[index]
    qsstart = int(qseqid.split(":")[2])#example (in brackets) = X:<INS>:(11703079):11710477:119253:1:PRECISE
    qsstop  = int(qseqid.split(":")[3])#example (in brackets) = X:<INS>:11703079:(11710477):119253:1:PRECISE
    qssize  = int(abs(qsstop - qsstart))
    sseqid  = df["sseqid"].values[index]
    if qssize <= size_et[sseqid] + 18:###
        indexs.append(index)

df = df.iloc[indexs]
print("shape by critere of size SV size TE :", df.shape)

df.to_csv(output, sep="\t", index=None)

#*********** HISTOGRAM *************

# print("BUILD HISTOGRAM IN :", dir_out + "HISTO")
# os.system("mkdir -p " + dir_out + "HISTO")

# reads_supports = []
# for i, qseqid in enumerate(df["qseqid"].values):
#     read_support = qseqid.split(":")[5]
#     reads_supports.append(int(read_support))
    
# list.sort(reads_supports)

# plt.grid(axis='y', alpha=0.95)

# bins = [x + 0.5 for x in range(0, max(reads_supports)+1)]
# maxe = max(reads_supports)
# if bins == 0:
#     bins = "auto"
#     maxe += 1
    
# plt.hist(reads_supports, bins=bins, color='#ff0000', alpha=0.7, rwidth=0.85) 
# plt.xlabel("Number of reads support")
# plt.ylabel("Number of TE")
# plt.xlim(min(reads_supports)-1, maxe+1)
# #plt.ylim(1, 80)
# plt.savefig(dir_out + "HISTO/"+ str(name_out) + "_TOTAL_histo.pdf")
# #plt.show() 


# for e, element in enumerate(df["sseqid"].unique()):
#     reads_supports = []
#     for i, qseqid in enumerate(df[df["sseqid"] == element]["qseqid"].values):
#         read_support = qseqid.split(":")[5]
#         reads_supports.append(int(read_support))

#     list.sort(reads_supports)
    
#     plt.grid(axis='y', alpha=0.95)
    
#     bins = [x + 0.5 for x in range(0, max(reads_supports)+1)]
#     maxe = max(reads_supports)
#     if bins == 0:
#         bins = "auto"
#         maxe += 1 
        
#     #print(element, bins, maxe)
#     plt.hist(reads_supports, bins=bins, color='#ff0000', alpha=0.7, rwidth=0.85) 
#     plt.xlabel("Number of reads support")
#     plt.ylabel("Number of TE")
#     plt.xlim(min(reads_supports)-1, maxe+1)
#     #plt.ylim(1, 80)
#     plt.savefig(dir_out + "HISTO/" + str(name_out) + "_" + str(element) + "_histo.pdf")
    #plt.show() 
    #
    
#*********** HISTOGRAM END *************


#*********** COUNT TE FOR GGPLOT *************

df_count = df.groupby(["sseqid"]).count()
d = {'x': [" "] * len(df_count.iloc[:, 0].values), 'y': df_count.index, 'z': df_count.iloc[:, 0].values}# x : Nothing, y : Name TE, z : Number of TE

df_o = pd.DataFrame(data=d)
df_o = df_o.sort_values(by="z")

print("build :", dir_out + name_out + "_hit_map.csv")
df_o.to_csv(dir_out + name_out + "_hit_map.csv", sep="\t", index=None)

print("Number total TE :", sum(df_o["z"].values))
print("list TE :")
df_o.columns = ["x", "TE", "NB_TE"]
print(df_o[["TE", "NB_TE"]])

