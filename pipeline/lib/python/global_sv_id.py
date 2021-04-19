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

parser.add_argument("output_file", metavar='<output_file_csv>', type=argparse.FileType('w'),
                    help="Name of tabular output file ")


#OPTION
parser.add_argument("-p", "--min-pident", dest='min_pident', type=int, default=80,
                    help="minimum percentage of identity")
parser.add_argument("-s", "--min-size-percent", dest="min_size_percent", type=int, default=90,
                    help="minimum percentage of TE size")
parser.add_argument("-c", "--combine", dest="combine", default=False, action='store_true',
                    help="Combine parts blast TE")
parser.add_argument("--combine_name", dest="combine_name", default="COMBINE_TE.csv",
                    help="Combine name output file")

args = parser.parse_args()


input_file_bln   = args.blast_file.name
db_te            = args.db_file_te.name

min_pident       = args.min_pident
min_size_percent = args.min_size_percent
combine          = args.combine
combine_name     = args.combine_name

output_file_csv  = args.output_file.name
name_out         = output_file_csv.split("/")[-1].split(".")[0]
dir_out          = "/".join(output_file_csv.split("/")[:-1])

df = pd.read_csv(input_file_bln, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
#display(df.head())

size_et = {}
file  = open(db_te, "r")
lines = file.readlines()

for i, l in enumerate(lines):
    if l[0] == ">":
        size_et[l[1:].strip()] = len(lines[i + 1].strip())

file.close()

#size_et = {'Idefix': 7411, '17.6': 7439, '1731': 4648, '297': 6995, '3S18': 6126, '412': 7567, 'aurora-element': 4263, 'Burdock': 6411, 'copia': 5143, 'gypsy': 7469, 'mdg1': 7480, 'mdg3': 5519, 'micropia': 5461, 'springer': 7546, 'Tirant': 8526, 'flea': 5034, 'opus': 7521, 'roo': 9092, 'blood': 7410, 'ZAM': 8435, 'GATE': 8507, 'Transpac': 5249, 'Circe': 7450, 'Quasimodo': 7387, 'HMS-Beagle': 7062, 'diver': 6112, 'Tabor': 7345, 'Stalker': 7256, 'gtwin': 7411, 'gypsy2': 6841, 'accord': 7404, 'gypsy3': 6973, 'invader1': 4032, 'invader2': 5124, 'invader3': 5484, 'gypsy4': 6852, 'invader4': 3105, 'gypsy5': 7369, 'gypsy6': 7826, 'invader5': 4038, 'diver2': 4917, 'Dm88': 4558, 'frogger': 2483, 'rover': 7318, 'Tom1': 410, 'rooA': 7621, 'accord2': 7650, 'McClintock': 6450, 'Stalker4': 7359, 'Stalker2': 7672, 'Max-element': 8556, 'gypsy7': 5486, 'gypsy8': 4955, 'gypsy9': 5349, 'gypsy10': 6006, 'gypsy11': 4428, 'gypsy12': 10218, 'invader6': 4885, 'Helena': 1317, 'HMS-Beagle2': 7220, 'Osvaldo': 1543}

df = df[df["sseqid"].isin(size_et.keys())]
print("keep only TE on list, size :", len(df.values))


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
#df = df.sort_values(by=["sseqid"])
df = df[["sseqid", "qseqid", "pident", "size_per", "size_el", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]]


print("size for min_size_percent>=" + str(min_size_percent) + ", min_pident>=" + str(min_pident), len(df.values))

#keep the TE with the highest score
best_score_match = []
best_score_match_index = []
chaine = ""
for index, row in enumerate(df.values):

    qseqid = df["qseqid"].values[index]
    #ID     = qseqid.split(":")[4]
    #df_tmp          = df[df["qseqid"].str.contains(">:[0-9]*:[0-9]*:" + str(ID).replace(".", "\.") + ":")]

    #chaine = str(df["qseqid"].values[index])
    #TODO INDICE 7 a regler
    chaine          = str(":".join(qseqid.split(":")[:7]))

    df_tmp          = df[df["qseqid"] == qseqid]

    maxe_bitscore   = max(df_tmp["bitscore"].values) 
    df_best_score   = df_tmp[df_tmp["bitscore"] == maxe_bitscore]
    
    if "INSERTION" in qseqid or "Deletion" in qseqid:
        chaine = df["qseqid"].values[index]
    else:
        chaine = df["qseqid"].values[index] + df["sseqid"].values[index]

    if "Assemblytics_w_541" in chaine:
        print(chaine)
        print(df["sseqid"].values[index], maxe_bitscore, df["bitscore"].values[index])
        
    if chaine not in best_score_match and ("INSERTION" in chaine or "Deletion" in chaine):
        sstart       = df["sstart"].values[index]
        send         = df["send"].values[index]

        chaine  = df_best_score["qseqid"].values[0] + df_best_score["sseqid"].values[0]

        if "INSERTION" in qseqid or "Deletion" in qseqid:
            chaine = df_best_score["qseqid"].values[0]
        else:
            chaine = df_best_score["qseqid"].values[0] + df_best_score["sseqid"].values[0]

        if send < sstart:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "-"
        else:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "+"

        qseqid       = df["qseqid"].values[index]
        info_qseqid  = qseqid.split(":")

        df_tmp        = df[df["qseqid"] == qseqid]
        maxe_bitscore = max(df_tmp["bitscore"].values)
        df_best_score = df_tmp[df_tmp["bitscore"] == maxe_bitscore]
        
        
        best_score_match.append(chaine)
        
        #only type SV choice 
        best_score_match_index.append(index)
    elif df["size_per"].values[index] == 100 and df["pident"].values[index] == 100 and not ("INSERTION" in chaine or "Deletion" in chaine) :
        #print(chaine)
        sstart       = df["sstart"].values[index]
        send         = df["send"].values[index]
        qstart       = df["qstart"].values[index]
        qend         = df["qend"].values[index]
        size_el      = df["size_el"].values[index]
        chaine       = df_best_score["qseqid"].values[0] + df_best_score["sseqid"].values[0]

        if send < sstart:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "-"
        else:
            df["qseqid"].values[index] = df["qseqid"].values[index] + ":" + "+"

        qseqid       = df["qseqid"].values[index]
        info_qseqid  = qseqid.split(":")

        df_tmp        = df[df["qseqid"] == qseqid]
        maxe_bitscore = max(df_tmp["bitscore"].values)
        df_best_score = df_tmp[df_tmp["bitscore"] == maxe_bitscore]
        #Assemblytics_b_1178:+:INSERTION::3L_RaGOO_RaGOO_RaGOO:8598269-8605760:+
        start_qseqid   = info_qseqid[5].split("-")[0]
        end_qseqid     = info_qseqid[5].split("-")[1]

        start_qseqid   = str(int(start_qseqid) + int(qstart))
        end_qseqid     = str(int(start_qseqid) + int(size_el))

        info_qseqid[5] = str(start_qseqid) + "-" + str(end_qseqid)
        df["qseqid"].values[index] = ":".join(info_qseqid)
        
        
        if chaine not in best_score_match:
            best_score_match.append(chaine)
        
        #only type SV choice 
        best_score_match_index.append(index)


df = df.iloc[best_score_match_index]
df = df.sort_values(by=["sseqid"])
#display(df)

df.to_csv(output_file_csv, sep="\t", index=None)

