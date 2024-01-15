import os
import sys
import pandas as pd
import math


df = pd.read_csv(sys.argv[1], sep="\t")

fileInit = open(sys.argv[2], "r")
generation = [ gen.strip().split(":")[1] for gen in fileInit.readlines() ]

gen_index = {g: i for i, g in enumerate(generation)}

def count_gen_between(gd, ge):
    return gen_index[ge] - gen_index[gd]


# Pré-traitement df
# df['gen_num'] = df['group'].apply(lambda x: int(x.strip("G")))
# df.sort_values(by=['name', 'gen_num'], inplace=True)


print("chrom\tx\ty1\ty2\tname\ttype1\ttype2\tgen_inter\tnb_gen")
for i, name in enumerate(df["name"].unique()):
    df_tmp = df[df["name"] == name]
    #print("name", name, df_tmp)
    tmp_generation = []
    for index, gen in enumerate(generation):
        if gen in df_tmp["group"].unique() :
            tmp_generation.append(gen)


    type_1 = "NONE"
    type_2 = "NONE"
    nb_gen = "NONE"
    gen_inter = "NONE"
    type_val1 = []
    type_val2 = []
    for index, gen in enumerate(tmp_generation):

        df_tmp_gen = df_tmp[df_tmp["group"] == gen]
        if index == 0:
            if int(tmp_generation[0].strip("G")) > 0 :
                num_gen = int(tmp_generation[0].strip("G"))
            else :
                num_gen = 1
            sum_1 = df_tmp_gen["y1"].values[0]/num_gen
            sum_2 = df_tmp_gen["y2"].values[0]/num_gen
        
        #print("init:", sum_1, sum_2)
        if index < len(tmp_generation) - 1 :
            df_tmp_gen_p = df_tmp[df_tmp["group"] == tmp_generation[index+1]]
            value1 = df_tmp_gen["y1"].values[0] if len(df_tmp_gen.values) > 0 else 0
            value_p1 = df_tmp_gen_p["y1"].values[0] if len(df_tmp_gen_p.values) > 0 else 0
            
            ponderation = int(tmp_generation[index+1].strip("G"))-int(tmp_generation[index].strip("G"))
            sum_1 += (((value_p1*100) - (value1*100))/ponderation)*count_gen_between(tmp_generation[index], tmp_generation[index+1])

            #print("p:", ponderation, sum_1, (value1*100), (value_p1*100), ((value_p1*100) - (value1*100)), count_gen_between(tmp_generation[index], tmp_generation[index+1]))
            
            if ((value_p1*100) - (value1*100))/ponderation > 0:
                type_val1.append(">")
            elif ((value_p1*100) - (value1*100))/ponderation == 0 :
                type_val1.append("=")
            else :
                type_val1.append("<")

            value2 = df_tmp_gen["y2"].values[0] if len(df_tmp_gen.values) > 0 else 0
            value_p2 = df_tmp_gen_p["y2"].values[0] if len(df_tmp_gen_p.values) > 0 else 0

            sum_2 += (((value_p2*100) - (value2*100))/ponderation)*count_gen_between(tmp_generation[index], tmp_generation[index+1])
            if ((value_p2*100) - (value2*100))/ponderation > 0:
                type_val2.append(">")
            elif ((value_p2*100) - (value2*100))/ponderation == 0 :
                type_val2.append("=")
            else :
                type_val2.append("<")

            type_1 = ""
            if "<" not in type_val1 and ">" in type_val1:
                type_1 = "+"
            elif ">" not in type_val1  and "<" in type_val1:
                type_1 = "-"
            elif "<" not in type_val1 and ">" not in type_val1 and "=" in type_val1:
                type_1 = "="
            else:
                type_1 = "~"


            type_2 = ""
            if "<" not in type_val2 and ">" in type_val2:
                type_2 = "+"
            elif ">" not in type_val2  and "<" in type_val2:
                type_2 = "-"
            elif "<" not in type_val2 and ">" not in type_val2 and "=" in type_val2:
                type_2 = "="
            else:
                type_2 = "~"
            
            nb_gen = str(len(tmp_generation))
            gen_inter = tmp_generation[0] + "-" + tmp_generation[-1]

    if len(type_val1):
        print("\t".join([name.split(":")[0], str(name.split(":")[1]), str(sum_1), str(sum_2), name, type_1, type_2, gen_inter, nb_gen]))
