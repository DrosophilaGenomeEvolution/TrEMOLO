import pandas as pd
import sys
import os
import random
import warnings
warnings.filterwarnings("ignore")

#name_TE_FOUND       =  "TE_FOUND.csv"
name_TE_FOUND       =  sys.argv[1]
df_TE_FOUND         =  pd.read_csv(name_TE_FOUND, sep="\t", header=None)
df_TE_FOUND.columns =  ["chrom", "start", "end", "ID", "IDR", "type"]
#display(df_TE_FOUND.head())

#name_INS  = "INSERTION.csv"
name_INS  = sys.argv[2]

df_INS    = pd.read_csv(name_INS, sep="\t")
#display(df_INS.head())


#name_POS_expension       = "pos_expension.csv"
name_POS_expension       = sys.argv[3]
file                     = open(name_POS_expension, "r")

output = sys.argv[4]

expension_found = False
if len(file.readlines()):
    expension_found = True
    df_POS_expension         = pd.read_csv(name_POS_expension, sep="\t", header=None)
    df_POS_expension.columns = ["chrom", "POS", "ID", "ID_asm", "ID_POS"]
    #display(df_POS_expension.head())


tab_type = []
tab_id   = [] 
tab_id_query = []
for i, v in enumerate(df_INS.values):
    qseqid = df_INS["qseqid"].values[i]
    Type   = qseqid.split(":")[2]
    ID_Query   = qseqid.split(":")[6]
    tab_type.append(Type)
    tab_id.append(qseqid.split(":")[4]+":"+qseqid.split(":")[5].split("-")[0])
    tab_id_query.append(ID_Query)
    
df_INS["type"] = tab_type
df_INS["ID"]   = tab_id
df_INS["ID_Query"]   = tab_id_query

#display(df_TE_FOUND[df_TE_FOUND["type"] == "Insertion"].head())

####Insertion
df_INS_ins = df_INS[df_INS["type"] == "INSERTION"]
print(df_INS_ins.shape)

tab_position_ins = []
tab_chrom = []
tab_ID = []
tab_ID_Query = []
for i, v in enumerate(df_INS_ins.values) :
    qseqid   = df_INS_ins["qseqid"].values[i]
    ID_asm   = qseqid.split(":")[0]
    TE       = df_INS_ins["sseqid"].values[i]
    ID_Query = df_INS_ins["ID_Query"].values[i]

    df_tmp = df_TE_FOUND[df_TE_FOUND["IDR"] == ID_asm]
    
    tab_position_ins.append(df_tmp["start"].values[0])
    tab_chrom.append(df_tmp["chrom"].values[0])
    tab_ID.append(TE + "|" + ID_asm + "|" + str(ID_Query))
    tab_ID_Query.append(ID_Query)

df_INS_ins["chrom"]    = tab_chrom
df_INS_ins["pos_ref"]  = tab_position_ins
df_INS_ins["ID_IGV"]   = tab_ID
df_INS_ins["ID_Query"] = tab_ID_Query

#display(df_INS_ins.head())
df_INS_ins.to_csv("pos_TE_insertion.csv", sep="\t", index=None)

#### Expension

#TE_FOUND_MODIF
df_TE_FOUND = df_TE_FOUND[df_TE_FOUND["type"] != "Insertion"]
#display(df_TE_FOUND.head())

df_INS = df_INS[df_INS["type"] != "INSERTION"]
    
#display(df_INS.head())
if expension_found:
    
    tab_te = []
    tab_ID = []
    tab_ID_Query = []
    for i, v in enumerate(df_POS_expension.values) :
        ID_POS = df_POS_expension["ID_POS"].values[i]
        ID_asm = df_POS_expension["ID_asm"].values[i]
        df_tmp = df_INS[df_INS["ID"] == ID_POS]
        TE = df_tmp["sseqid"].values[0]
        ID_Query = df_tmp["ID_Query"].values[0]
        
        tab_te.append(TE)
        tab_ID.append(TE + "|" + ID_asm + "|" + str(ID_Query))
        tab_ID_Query.append(ID_Query)
        
    df_POS_expension["TE"] = tab_te
    df_POS_expension["ID_IGV"] = tab_ID
    df_POS_expension["ID_Query"] = tab_ID_Query
    #display(df_POS_expension.head())

    df_POS_expension.to_csv(output, sep="\t", index=None)
