import os
import sys
import pandas as pd
import math


dfM = pd.read_csv(sys.argv[1], sep="\t")
dfS = pd.read_csv(sys.argv[2], sep="\t")
type_col1 = []
type_col2 = []

for i, name in enumerate(dfS["name"].values):
    if len(dfM[dfM["name"] == name].values):
        type_col1.append(dfM[dfM["name"] == name]["type1"].values[0])
        type_col2.append(dfM[dfM["name"] == name]["type2"].values[0])
    else:
        type_col1.append("NONE")
        type_col2.append("NONE")

dfS["type1"] = type_col1
dfS["type2"] = type_col2

dfS = dfS.loc[:, ["chrom", "x", "y1", "y2", "group", "name", "type1", "type2", "IN_OUT", "ID"]]

dfS.to_csv(sys.argv[3], sep="\t", index=False)