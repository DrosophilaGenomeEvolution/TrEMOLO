import pandas as pd
import os
import sys
import re
import argparse


parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("blast_file", metavar='<blast-file>', option_strings=['blast-file'], type=argparse.FileType('r'),
                    help="Input blast file format 6 (-outfmt 6)")

parser.add_argument("db_file_te", metavar='<db-file-TE>', type=argparse.FileType('r'),
                    help="Multi fasta file TE for get size")

parser.add_argument("output_file", metavar='<output>', type=argparse.FileType('w'),
                    help="Name of tabular output file ")


#OPTION
parser.add_argument("-p", "--min-pident", dest='min_pident', type=int, default=80,
                    help="minimum percentage of identity")
parser.add_argument("-s", "--min-size-percent", dest="min_size_percent", type=int, default=80,
                    help="minimum percentage of TE size")
parser.add_argument("-c", "--combine", dest="combine", default=False, action='store_true',
                    help="Combine parts blast TE")
parser.add_argument("--combine_name", dest="combine_name", default="COMBINE_TE.csv",
                    help="Combine name file out")

args = parser.parse_args()


input_file_bln   = args.blast_file.name
db_te            = args.db_file_te.name

min_pident       = args.min_pident
min_size_percent = args.min_size_percent
combine          = args.combine
combine_name     = args.combine_name

output_file_csv  = args.output_file.name
name_out         = args.output_file.name.split("/")[-1].split(".")[0]
dir_out          = "/".join(args.output_file.name.split("/")[:-1])

df = pd.read_csv(input_file_bln, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

size_et = {}
file    = open(db_te, "r")
lines   = file.readlines()

for i, l in enumerate(lines):
    if l[0] == ">":
        size_et[l[1:].strip()] = len(lines[i + 1].strip())

file.close()

df = df[df["sseqid"].isin(size_et.keys())]
print("[" + sys.argv[0] + "]", "keep only TE on list :", str(len(df.values)))


#keep the TE with the highest score
best_score_match       = []
best_score_match_index = []
chaine                 = ""
for index, row in enumerate(df.values):
    
    qseqid = df["qseqid"].values[index]

    #chaine = str(df["qseqid"].values[index])
    #TODO INDICE 7 a regler
    chaine          = str(":".join(qseqid.split(":")[:7]))
    df_tmp          = df[df["qseqid"] == qseqid]

    maxe_bitscore   = max(df_tmp["bitscore"].values) 
    df_best_score   = df_tmp[df_tmp["bitscore"] == maxe_bitscore]

    if chaine not in best_score_match :

        #sstart       = df["sstart"].values[index]
        #send         = df["send"].values[index]
        
        sstart       = df_best_score["sstart"].values[0]
        send         = df_best_score["send"].values[0]
        
        info_qseqid  = qseqid.split(":")

        index_tmp = df_best_score.index
        condition = df_best_score["bitscore"] == maxe_bitscore
        max_index = index_tmp[condition]

        indice    = max_index.tolist()[0]

        if send < sstart :
            df["qseqid"].values[indice] = df["qseqid"].values[indice] + ":" + "-"
        else :
            df["qseqid"].values[indice] = df["qseqid"].values[indice] + ":" + "+"
            
        #chaine = df_best_score["qseqid"].values[0] + df_best_score["sseqid"].values[0]
        best_score_match.append(chaine)

        #only type SV choice 
        best_score_match_index.append(indice)

    

df_best = df.iloc[best_score_match_index]


print("[" + sys.argv[0] + "]", " FILTER BEST MATCH TE : ", str(len(df_best.values)))
print("[" + sys.argv[0] + "]", "Combine TE...")
#COMBINE
dic_comb = {"qseqid": [], "sseqid": [], "grain_pident":[], "qstart": [], "qend": [], "sstart": [], "send": []}
for i, v in enumerate(df_best.values) :

    qseqid = df_best["qseqid"].values[i]
    sseqid = df_best["sseqid"].values[i]
    pident = df_best["pident"].values[i]
    df_tmp = df[df["qseqid"].str.contains(qseqid[:-2].replace(":", "\:").replace("+", "\+").replace("-", "\-"))]
    df_tmp = df_tmp[df_tmp["sseqid"] == sseqid]

    sstart_best = df_tmp["sstart"].values[0]
    ssend_best  = df_tmp["send"].values[0]
    qstart_best = df_tmp["qstart"].values[0]
    qend_best   = df_tmp["qend"].values[0]

    if int(sstart_best) < int(ssend_best):#forward
        df_tmp = df_tmp.sort_values(by=["sstart"])

        # if "2R_RaGOO_RaGOO_RaGOO:2371284-2378557" in qseqid :
        #     print(qseqid)
        #     print(df_tmp)

        sstart_global = sstart_best
        send_global   = ssend_best
        qstart_global = qstart_best
        qend_global   = qend_best

        for a, x in enumerate(df_tmp.values):
            sstart = df_tmp["sstart"].values[a]
            ssend  = df_tmp["send"].values[a]
            qstart = df_tmp["qstart"].values[a]
            qend   = df_tmp["qend"].values[a]
            if sstart < sstart_global and send < send_global and qstart < qstart_global and qend < qend_global :
                sstart_global = sstart
                qstart_global = qstart
            elif sstart > sstart_global and send > send_global and qstart > qstart_global and qend > qend_global:
                send_global   = send
                qend_global   = qend

        sstart     = sstart_global
        send       = send_global
        qstart_min = qstart_global
        qend_max   = qend_global
        
    else:#reverse
        
        sstart_global = sstart_best
        send_global   = ssend_best
        qstart_global = qstart_best
        qend_global   = qend_best

        for a, x in enumerate(df_tmp.values):
            sstart = df_tmp["sstart"].values[a]
            ssend  = df_tmp["send"].values[a]
            qstart = df_tmp["qstart"].values[a]
            qend   = df_tmp["qend"].values[a]
            if sstart > sstart_global and send > send_global and qstart < qstart_global and qend < qend_global :
                sstart_global = sstart
                qstart_global = qstart
            elif sstart < sstart_global and send < send_global and qstart > qstart_global and qend > qend_global:
                send_global   = send
                qend_global   = qend


        sstart     = sstart_global
        send       = send_global
        qstart_min = qstart_global
        qend_max   = qend_global

    dic_comb["qseqid"].append(qseqid)
    dic_comb["sseqid"].append(sseqid)
    dic_comb["grain_pident"].append(pident)
    dic_comb["qstart"].append(qstart_min)
    dic_comb["qend"].append(qend_max)
    dic_comb["sstart"].append(sstart)
    dic_comb["send"].append(send)
        
df_comb = pd.DataFrame(dic_comb)
print("[" + sys.argv[0] + "]", "COUNT TE COMBINE :", str(len(df_comb.values)))

#CALCUL SIZE AND PERCENT SIZE
tab_percent = []
tab_size    = []
for index, row in enumerate(df_comb[["sseqid", "qend", "qstart"]].values):
    size_element = abs(int(row[1])-int(row[2]))#qend - qstart
    tab_percent.append(round((size_element/size_et[row[0]]) * 100, 1))
    tab_size.append(size_element)


df_comb["size_per"] = tab_percent
df_comb["size_el"]  = tab_size

df_comb = df_comb[["sseqid", "qseqid", "grain_pident", "size_per", "size_el", "qstart", "qend", "sstart", "send"]]
df_comb = df_comb.sort_values(by="sseqid")

df_comb = df_comb[df_comb["size_per"] >= min_size_percent]#min_size_v2 combine

df_comb.to_csv(combine_name, sep="\t", index=None)

df = df[df["qseqid"].isin(df_comb["qseqid"].values)]

tab_percent = []
tab_size    = []
for index, row in enumerate(df[["sseqid", "qend", "qstart"]].values):
    size_element = abs(int(row[1])-int(row[2]))#qend - qstart
    tab_percent.append(round((size_element/size_et[row[0]]) * 100, 2))
    tab_size.append(size_element)
    
#keeps the TE according to the percentage of identity and the percentage of size
df["size_per"] = tab_percent
df["size_el"]  = tab_size
#df = df[df["size_per"] >= min_size_percent]
df = df[df["pident"] >= min_pident]
df = df.sort_values(by=["sseqid"])
df = df[["sseqid", "qseqid", "pident", "size_per", "size_el", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]]

print("[" + sys.argv[0] + "]", "TE with min_size_percent>=" + str(min_size_percent) + ", min_pident>=" + str(min_pident), " :" + str(len(df.values)))

df.to_csv(output_file_csv, sep="\t", index=None)


