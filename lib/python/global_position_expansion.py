import pandas as pd
import sys
import os
import random

#name_INS = "position_Repeat_expension_best_score.csv"
name_INS   = sys.argv[1]
output_pos = sys.argv[2]

df = pd.read_csv(name_INS, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
#print(df.head())


tab_ID = []
tab_region = []
for i, v in enumerate(df.values):
    qseqid = df["qseqid"].values[i]
    sseqid = df["sseqid"].values[i]
    sstart = df["sstart"].values[i]
    qend   = df["qend"].values[i]
    qstart = df["qstart"].values[i]
    qend   = df["qend"].values[i]
    tab_ID.append(qseqid.split(":")[0])
    tab_region.append(qseqid.split(":")[1])
    
df["ID"]     = tab_ID
df["region"] = tab_region

#print(tab_ID[:5])
tab_ind = []
#print(df.head())
dico_ID = {}
dico_df_out = {"chrom":[], "position_TE":[], "ID":[]}

size_match_kepp = 3000
for i, v in enumerate(df.values):
    qseqid = df["qseqid"].values[i]
    sseqid = df["sseqid"].values[i]
    sstart = df["sstart"].values[i]
    send   = df["send"].values[i]
    qstart = df["qstart"].values[i]
    qend   = df["qend"].values[i]
    ID     = df["ID"].values[i]
    region = df["region"].values[i]
    
    df_tmp_ID = df[df["ID"] == ID]
    if region == "LEFT" and len(df_tmp_ID.values) >= 1 :
        df_tmp_region = df_tmp_ID[df_tmp_ID["region"] == "RIGTH"]
        if len(df_tmp_region.values) >= 1 :
            if ID not in dico_ID:
                dico_ID[ID] = 0
            else :
                dico_ID[ID] += 1

            #print(dico_ID)
            sstart_region = df_tmp_region["sstart"].values[0]
            qend_region   = df_tmp_region["qend"].values[0]
            qstart_region = df_tmp_region["qstart"].values[0]
            sseqid_region = df_tmp_region["sseqid"].values[0]
            #print(df_tmp_region)
            if abs(send - sstart_region) <= 100 or abs(qend - qstart) > size_match_kepp or abs(qend_region-qstart_region) > size_match_kepp :
                tab_ind.append(i)
                tab_ind.append(i+1)
                if abs(qend - qstart) > size_match_kepp :
                    pos = send
                else :
                    pos = sstart_region
                    sseqid = sseqid_region

                id_pos = int(qseqid.split(":")[4].split("-")[1])
                dico_df_out["ID"].append(qseqid.split(":")[0]+":"+qseqid.split(":")[3]+":"+str(id_pos))
                dico_df_out["position_TE"].append(pos)
                dico_df_out["chrom"].append(sseqid)
            
    #else:
        #df_tmp_region = df_tmp_ID[df_tmp_ID["region"] == "LEFT"]

        

df     = df.iloc[tab_ind]

if len(df.values):
    #print(df.head())
    df.to_csv("position_Repeat_expension_best_score_size.csv", sep="\t", index=None)
else:
    print("ERROR : position_Repeat_expension_best_score_size.csv")

df_out = pd.DataFrame(dico_df_out)

if len(df_out.values):
    df_out.to_csv(output_pos, sep="\t", index=None)
else:
    print("ERROR : " + str(output_pos))
