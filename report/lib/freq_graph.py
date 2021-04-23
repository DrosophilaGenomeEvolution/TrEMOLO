import pandas as pd
import sys
import os


name_input  = sys.argv[1]
name_output = sys.argv[2]

df  = pd.read_csv(name_input, sep="\t")

tab_per = ["0-2", "2-9", "9-15", "15-50", "50-100"]
dico = {"x":[], "y":[], "condition":[]}

for i, v in enumerate(df["TE"].unique()):
    TE       = v

    COUNT_TE = 0
    df_tmp = df[df["TE"] == TE]
    #print(TE, "TE")
    for e in tab_per:
        mini = float(e.split("-")[0])
        maxi = float(e.split("-")[1])
        
        df_tmp_mini =  df_tmp[df_tmp["read_support_percent"] > mini]
        if len(df_tmp_mini.values) :
            df_tmp_maxi = df_tmp_mini[df_tmp_mini["read_support_percent"] <= maxi]
            if len(df_tmp_maxi.values) :
                group = str(int(mini))+"-"+str(int(maxi))+"(%)"
                COUNT_TE = len(df_tmp_maxi.values)
                #print("ok", len(df_tmp_maxi.values))
                #print(TE, group, COUNT_TE, i)
                dico["x"].append(TE)
                dico["condition"].append(group)
                dico["y"].append(COUNT_TE)

df_out = pd.DataFrame(dico)
df_out.to_csv(name_output, sep="\t", index=None)
    
