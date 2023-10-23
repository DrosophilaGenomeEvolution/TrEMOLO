import pandas as pd
import sys
import re
import os


name_input = sys.argv[1]
name_var   = sys.argv[2]

df = pd.read_csv(name_input, sep=",")

liste = []
for i, v in enumerate(df.values):
    dico = {}
    for e, w in enumerate(df.columns):
        dico[w] = str(v[e])
        
    liste.append(dico)

print("var " + name_var + " = " + str(liste).replace("'", "\""))
print(name_var + "[\"columns\"] = " + str(list(df.columns)).replace("'", "\""))


