#!/home/mourdas/anaconda3/bin/python3
import pandas as pd
import sys
import re
import os
from collections import OrderedDict
import json
import argparse

parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("input", metavar='<input-file>', option_strings=['input-file'], type=argparse.FileType('r'),
                    help="Input blast file format 6 (-outfmt 6)")

parser.add_argument("var", metavar='<var>', type=str,
                    help="Multi fasta file TE for get size")

#OPTION
parser.add_argument("-c", "--order-cond", dest='order_c', type=str, default="NONE",
                    help="minimum percentage of identity")
parser.add_argument("-x", "--order-x", dest="order_x", type=str, default="NONE",
                    help="minimum percentage of TE size")

args = parser.parse_args()

name_input = args.input.name
name_var   = args.var

order_x = args.order_x
order_c = args.order_c

df = pd.read_csv(name_input, sep="\t")

if order_c != "NONE":

    file  = open(order_c, "r")
    lines = [line.strip() for line in file.readlines()]

     # Create the dictionary that defines the order for sorting
    sorterIndex = dict(zip(lines, range(len(lines))))

    df['rank'] = df['x'].map(sorterIndex)

    df.sort_values(['rank'], ascending = [True], inplace = True)
    df.drop('rank', 1, inplace = True)

    file.close()

    order_c = []
    for line in lines:
        if line in list(df["condition"].unique()):
            order_c.append(line)

else :
    order_c = sorted(list(df["condition"].unique()))


if order_x == "NONE" :

    already_work = []
    list_count   = []
    for i, v in enumerate(df[["x", "y"]].values):
        TE     = v[0]
        value  = v[1]
        dico   = {"group": TE}
        df_tmp = df[df["x"]==TE]

        list_count.append(sum(df_tmp["y"].values))
        already_work.append(TE)

    df["total"] = list_count

    df = df.sort_values(by=["total"])

else :

    file  = open(order_x, "r")
    lines = [line.strip() for line in file.readlines()]

    # Create the dictionary that defines the order for sorting
    sorterIndex = dict(zip(lines, range(len(lines))))

    df['rank'] = df['x'].map(sorterIndex)

    df.sort_values(['rank'], ascending = [True], inplace = True)
    df.drop('rank', 1, inplace = True)

    file.close()


already_work = []
liste        = []

for i, v in enumerate(df[["x", "y"]].values):
    TE     = v[0]
    value  = v[1]
    #dico   = {"group": TE}
    dico   = OrderedDict()
    dico["group"] = TE
    df_tmp = df[df["x"] == TE]
    if TE not in already_work:
        df_tmp = df_tmp.sort_values(by=["condition"])
        for e, w in enumerate(df_tmp.values[::-1]):
            dico[df_tmp["condition"].values[e]] = str(df_tmp["y"].values[e])
        already_work.append(TE)
        liste.append(json.dumps(dico))

print("var " + name_var + " = " + str(liste).replace("'", ""))
print(name_var + "[\"columns\"] = " + str(["group"] + order_c).replace("'", "\""))
