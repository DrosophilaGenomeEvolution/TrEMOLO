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
best_score_match            = []
best_score_match_index      = []
best_score_match_index_comb = []
chaine                      = ""

size_df    = len(df.values)

dfs = df.sort_values(by=["bitscore"], ascending=False)
df["index"] = [0] * len(df.values)
for index, row in enumerate(dfs.values):
    
    qseqid = dfs["qseqid"].values[index]

    #TODO INDICE 5 a regler
    chaine          = str(":".join(qseqid.split(":")[:5]))

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


print("[" + sys.argv[0] + "]", "FILTER BEST MATCH TE : ", str(len(df_best.values)))
print("[" + sys.argv[0] + "]", "COMBINE TE...")
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


