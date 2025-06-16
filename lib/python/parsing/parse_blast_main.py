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
    usage: parse_blast_main.py [-h] [-p MIN_PIDENT] [-s MIN_SIZE_PERCENT] [-r MIN_READ_SUPPORT] [-k TYPE_SV_KEEP] [-c] [--combine_name COMBINE_NAME] <blast-file> <db-file-TE> <output>

    filters a blast file in output format 6 to keep the candidate candidate TE active

    positional arguments:
      <blast-file>          Input blast file format 6 (-outfmt 6)
      <db-file-TE>          Multi fasta file TE for get size
      <output>              Name of tabular output file

    optional arguments:
      -h, --help            show this help message and exit
      -p MIN_PIDENT, --min-pident MIN_PIDENT
                            minimum percentage of identity (default: 94)
      -s MIN_SIZE_PERCENT, --min-size-percent MIN_SIZE_PERCENT
                            minimum percentage of TE size (default: 90)
      -r MIN_READ_SUPPORT, --min-read-support MIN_READ_SUPPORT
                            minimum read support number (default: 1)
      -k TYPE_SV_KEEP, --type-sv-keep TYPE_SV_KEEP
                            Type of SV (default: INS)
      -c, --combine         Combine parts blast TE (default: False)
      --combine_name COMBINE_NAME
                            Combine name output file (default: COMBINE_TE.csv)
"""


#import numpy
import pandas as pd
import re
import os
import sys
import argparse
import datetime



parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("blast_file", metavar='<blast-file>', option_strings=['blast-file'], type=argparse.FileType('r'),
                    help="Input blast file format 6 (-outfmt 6)")

parser.add_argument("db_file_te", metavar='<db-file-TE>', type=argparse.FileType('r'),
                    help="index (.fai) of multi fasta file TE for get size")

parser.add_argument("output_file", metavar='<output>', type=argparse.FileType('w'),
                    help="Name of tabular output file ")


#OPTION
parser.add_argument("-p", "--min-pident", dest='min_pident', type=int, default=94,
                    help="minimum percentage of identity")
parser.add_argument("-s", "--min-size-percent", dest="min_size_percent", type=int, default=90,
                    help="minimum percentage of TE size")
parser.add_argument("-r", "--min-read-support", dest="min_read_support", type=int, default=1,
                    help="minimum read support number")
parser.add_argument("-k", "--type-sv-keep", dest="type_sv_keep", type=str, default="INS",
                    help="Type of SV")
parser.add_argument("-c", "--combine", dest="combine", default=False, action='store_true',
                    help="Combine parts blast TE")
parser.add_argument("--combine_name", dest="combine_name", default="COMBINE_TE.csv",
                    help="Combine name output file")
args = parser.parse_args()


name_file        = args.blast_file.name
min_pident       = args.min_pident
min_size_percent = args.min_size_percent
min_read_support = args.min_read_support
type_sv_keep     = args.type_sv_keep
combine          = args.combine
combine_name     = args.combine_name

output           = args.output_file.name
name_out         = output.split("/")[-1].split(".")[0]
dir_out          = "/".join(output.split("/")[:-1])

print(f"[{str(sys.argv[0])}] : PREFIX OUTPUT : {name_out}")


if dir_out != "":
    dir_out += "/"
else :
    dir_out += "./"


size_et = {}

if args.db_file_te != None:
    #GET_SIZE TE Database format fasta
    file  = args.db_file_te
    lines = file.readlines()

    for i, line in enumerate(lines):
        TE, size = line.split("\t")[:2]
        size_et[TE] = int(size)

    file.close()

else :
    print(f"[{str(sys.argv[0])}] : ERROR fasta file Database TE Not Found")
    exit(1)

#GET BLAST OUTFMT 

df = pd.read_csv(filepath_or_buffer=name_file, sep="\t", header=None)


df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

#KEEP ONLY TE on the LIST (GET_SIZE CANNONICAL TE)
df = df[df["sseqid"].isin(size_et.keys())]
df.index = [int(i) for i in range(0, len(df.values))]
print("[" + str(sys.argv[0]) + "] : KEEP ONLY TE on DB, NUMBER :", len(df.values))

#IDENTITY
#df = df[df["pident"] >= min_pident]

#df = df[df["evalue"] == 0]

print("[" + str(sys.argv[0]) + "] : FILTER BEST bitscore...")
#keep the TE with the highest score
best_score_match = []
best_score_match_index = []
best_score_match_index_comb = []

chaine = ""

size_df    = len(df.values)
##TODO add sort by length maybe
dfs = df.sort_values(by=["bitscore"], ascending=False)
df["index"] = [0] * len(df.values)
for index, row in enumerate(dfs.values):
    
    qseqid = dfs["qseqid"].values[index]
    ID     = qseqid.split(":")[4]

    #ex 2R:<INS>:12536774:12536775:TrEMOLO.INS.2785:57:IMPRECISE:25:+ to 2R:<INS>:12536774:12536775:TrEMOLO.INS.2785:57:IMPRECISE
    chaine          = str(":".join(qseqid.split(":")[:7]))
    
    if chaine not in best_score_match :

        index_sort = dfs.index.values[index]
        best_score_match.append(chaine)

        sseqid       = df.iloc[index_sort]["sseqid"]
        qseqid       = df.iloc[index_sort]["qseqid"]

        sstart       = int(df.iloc[index_sort]["sstart"])
        send         = int(df.iloc[index_sort]["send"])

        if send < sstart :
            df.iloc[index_sort, df.columns.get_loc("qseqid")] = df.iloc[index_sort]["qseqid"] + ":" + "-"
        else :
            df.iloc[index_sort, df.columns.get_loc("qseqid")] = df.iloc[index_sort]["qseqid"] + ":" + "+"

        #GET RS
        info_qseqid  = qseqid.split(":")
        read_support = int(info_qseqid[5])

        #only type SV choice 
        find_type_SV = re.search(type_sv_keep, qseqid)
        if find_type_SV and read_support >= min_read_support:
            best_score_match_index.append(index_sort)
            df.loc[index_sort, "index"] = index_sort
            
            best_score_match_index_comb.append(index_sort)
            i = index_sort + 1
            while i < size_df and qseqid == df.iloc[i]["qseqid"] and sseqid == df.iloc[i]["sseqid"]:

                if send < sstart :
                    df.iloc[i, df.columns.get_loc("qseqid")] = df.iloc[i]["qseqid"] + ":" + "-"
                else :
                    df.iloc[i, df.columns.get_loc("qseqid")] = df.iloc[i]["qseqid"] + ":" + "+"

                df.loc[i, "index"] = i
                best_score_match_index_comb.append(i)
                i += 1


df_best = df.iloc[best_score_match_index]

df_b_to_combine = df.iloc[best_score_match_index_comb]
print("[" + str(sys.argv[0]) + "] : NUMBER BEST TE SCORE MATCH :", len(df_best.values))
#df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

if combine :
    
    print("[" + str(sys.argv[0]) + "] : COMBINE BLAST TE...")
    #COMBINE
    dic_comb = {"qseqid": [], "sseqid": [], "grain_pident":[], "qstart": [], "qend": [], "sstart": [], "send": []}
    qseqid_p = ""
    sseqid_p = ""
    pident_p = 0
    for i, v in enumerate(df_b_to_combine.values) :

        qseqid = df_b_to_combine["qseqid"].values[i]
        sseqid = df_b_to_combine["sseqid"].values[i]
        pident = df_b_to_combine["pident"].values[i]

        if qseqid_p != qseqid :
            if qseqid_p != "" :
                dic_comb["qseqid"].append(qseqid_p)
                dic_comb["sseqid"].append(sseqid_p)
                dic_comb["grain_pident"].append(pident_p)
                dic_comb["qstart"].append(qstart_global)
                dic_comb["qend"].append(qend_global)
                dic_comb["sstart"].append(sstart_global)
                dic_comb["send"].append(send_global)

            qseqid_p = qseqid
            sseqid_p = sseqid
            pident_p = pident

            sstart_best = df_b_to_combine["sstart"].values[i]
            send_best   = df_b_to_combine["send"].values[i]
            qstart_best = df_b_to_combine["qstart"].values[i]
            qend_best   = df_b_to_combine["qend"].values[i]

            sstart_global = sstart_best
            send_global   = send_best
            qstart_global = qstart_best
            qend_global   = qend_best
        else :
            if int(sstart_best) < int(send_best):#forward

                sstart = df_b_to_combine["sstart"].values[i]
                send   = df_b_to_combine["send"].values[i]
                qstart = df_b_to_combine["qstart"].values[i]
                qend   = df_b_to_combine["qend"].values[i]
                if sstart < sstart_global and send < send_global and qstart < qstart_global and qend < qend_global :
                    sstart_global = sstart
                    qstart_global = qstart
                elif sstart > sstart_global and send > send_global and qstart > qstart_global and qend > qend_global:
                    send_global   = send
                    qend_global   = qend
                
            else:#reverse

                sstart = df_b_to_combine["sstart"].values[i]
                send   = df_b_to_combine["send"].values[i]
                qstart = df_b_to_combine["qstart"].values[i]
                qend   = df_b_to_combine["qend"].values[i]
                if sstart > sstart_global and send > send_global and qstart < qstart_global and qend < qend_global :
                    sstart_global = sstart
                    qstart_global = qstart
                elif sstart < sstart_global and send < send_global and qstart > qstart_global and qend > qend_global:
                    send_global   = send
                    qend_global   = qend

            #remove :+/:- at no best value
            df.iloc[df_b_to_combine["index"].values[i], df.columns.get_loc("qseqid")] = df_b_to_combine["qseqid"].values[i][:-2]

    if qseqid_p != "" :   
        dic_comb["qseqid"].append(qseqid_p)
        dic_comb["sseqid"].append(sseqid_p)
        dic_comb["grain_pident"].append(pident_p)
        dic_comb["qstart"].append(qstart_global)
        dic_comb["qend"].append(qend_global)
        dic_comb["sstart"].append(sstart_global)
        dic_comb["send"].append(send_global)

    df_comb = pd.DataFrame(dic_comb)

    #CALCUL SIZE AND PERCENT SIZE
    tab_percent = []
    tab_size    = []
    for index, row in enumerate(df_comb[["sseqid", "qend", "qstart"]].values):
        size_element = abs(int(row[1])-int(row[2]))#qend - qstart
        tab_percent.append(round((size_element/size_et[row[0]]) * 100, 2))
        tab_size.append(size_element)

    df_comb["size_per"] = tab_percent
    df_comb["size_el"]  = tab_size
    df_comb = df_comb[df_comb["size_per"] >= min_size_percent]#min_size_v2 combine

    df_comb = df_comb[["sseqid", "qseqid", "grain_pident", "size_per", "size_el", "qstart", "qend", "sstart", "send"]]
    df_comb = df_comb.sort_values(by=["sseqid", "qseqid"])


    df_comb.to_csv(path_or_buf=combine_name, sep="\t", index=None)


    tab_percent = []
    tab_size    = []
    for index, row in enumerate(df_comb[["sseqid", "qend", "qstart"]].values):
        size_element = abs(int(row[1])-int(row[2]))#qend - qstart
        tab_percent.append(round((size_element/size_et[row[0]]) * 100, 2))
        tab_size.append(size_element)


    df_comb["size_per"] = tab_percent
    df_comb["size_el"]  = tab_size
    df_comb             = df_comb[df_comb["size_per"] >= min_size_percent]


##KEEP COMBINE OR NOT
if combine :
    df = df[df["qseqid"].isin(df_comb["qseqid"].values)]
else :
    df = df_best

#CALCUL SIZE AND PERCENT SIZE OF TE
tab_percent = []
tab_size    = []
for index, row in enumerate(df[["sseqid", "qend", "qstart"]].values):
    size_element = abs(int(row[1])-int(row[2]))#qend - qstart
    tab_percent.append(round((size_element/size_et[row[0]]) * 100, 2))
    tab_size.append(size_element)
    
#keeps the TE according to the percentage of identity and the percentage of size
df["size_per"] = tab_percent
df["size_el"]  = tab_size

if not combine :
    df = df[df["size_per"] >= min_size_percent]

df = df[df["pident"] >= min_pident]
df = df.sort_values(by=["sseqid"])
df = df[["sseqid", "qseqid", "pident", "size_per", "size_el", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]]
#df = df[["sseqid", "qseqid", "size_per", "size_el", "qstart", "qend", "sstart", "send"]]


print("[" + str(sys.argv[0]) + "] : NUMBER TE with minimum size percentage>=" + str(min_size_percent) + ", minimum percent identity>=" + str(min_pident)+" :", len(df.values))
df.to_csv(path_or_buf=output, sep="\t", index=None)

#*********** COUNT TE FOR GGPLOT *************

df_count = df.groupby(["sseqid"]).count()
d = {'x': [" "] * len(df_count.iloc[:, 0].values), 'y': df_count.index, 'z': df_count.iloc[:, 0].values}# x : Nothing, y : Name TE, z : Number of TE

df_o = pd.DataFrame(data=d)
df_o = df_o.sort_values(by="z")

print("[" + str(sys.argv[0]) + "] : BUILD :", dir_out + name_out + "_COUNT.csv")
df_o.to_csv(path_or_buf=dir_out + name_out + "_COUNT.csv", sep="\t", index=None)

print("[" + str(sys.argv[0]) + "] : NUMBER TE TOTAL :", sum(df_o["z"].values))
print("[" + str(sys.argv[0]) + "] : LIST of TE :")
df_o.columns = ["x", "TE", "NB_TE"]
print(df_o[["TE", "NB_TE"]].to_string(index=False))


