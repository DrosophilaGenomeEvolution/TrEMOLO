#!/home/mourdas/anaconda3/bin/python3
import pandas as pd
import sys
import re
import os
from collections import OrderedDict


name_input = sys.argv[1]
name_var   = sys.argv[2]

df = pd.read_csv(name_input, sep="\t")

#print(df.head())
already_work = []
list_count = []
for i, v in enumerate(df[["x", "y"]].values):
    TE     = v[0]
    value  = v[1]
    dico   = {"group": TE}
    df_tmp = df[df["x"]==TE]

    list_count.append(sum(df_tmp["y"].values))
    already_work.append(TE)

df["total"] = list_count

df = df.sort_values(by=["total"])

already_work = []
liste        = []

for i, v in enumerate(df[["x", "y"]].values):
    TE     = v[0]
    value  = v[1]
    dico   = {"group": TE}
    #dico   = OrderedDict()
    #dico["group"] = TE
    df_tmp = df[df["x"] == TE]
    
    if TE not in already_work:
        df_tmp = df_tmp.sort_values(by=["condition"], ascending=True)
        for e, w in enumerate(df_tmp.values):
            dico[df_tmp["condition"].values[e]] = str(df_tmp["y"].values[e])

        already_work.append(TE)
        liste.append(dico)

print("var " + name_var + " = " + str(liste).replace("'", "\""))
print(name_var + "[\"columns\"] = " + str(["group"] + list(df["condition"].unique())).replace("'", "\""))


