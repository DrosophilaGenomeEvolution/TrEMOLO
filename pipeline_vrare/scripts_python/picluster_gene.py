#!/usr/bin python3
# -*- coding: utf-8 -*-

###################################################################################################################################
#
# Copyright 2019-2020 IRD-CNRS-Lyon1 University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to authors and IRD, CNRS, and Lyon 1 University  for all versions
# Version 0.1 written by Mourdas Mohamed
#                                                                                                                                   
####################################################################################################################################



import sys
import pandas as pd
import re
import os
import argpars


parser = argparse.ArgumentParser(description="Get picluster with flanking genes")

#MAIN ARGS
parser.add_argument("genome", type=str,
                    help="genome fasta file")
parser.add_argument("genes", type=str,
                    help="blastn file format 6 flanking genes against genome ")
parser.add_argument("elements", type=str,
                    help="blastn file format 6 TE against SV (Structural variant) ")
parser.add_argument("flank_genes", type=str,
                    help="name output fasta file")
parser.add_argument("all_te", type=str,
                    help="tabular of all potentiel candidate of TE generate by vrare.snk")


#OPTION
parser.add_argument("-c", "--chrom", type=str, default='2L,2R,3L,3R,4,X_',
                    help='chromosome (or part) to keep (give a list of arguments separate the values ​​with commas "X,Y") [2L,2R,3L,3R,4,X_]')

args = parser.parse_args()

##
genome_fasta          = args.genome
name_file_gene_bln    = args.genes
name_file_element_bln = args.elements
gene_flank_file       = args.flank_genes
name_file_all_TE      = args.all_te


chrom_tab = args.chrom.split(",")#KEEP ONLY

### GET SIZE GENOME

#check if an element of the array is at least part of the value
def regex_in_list(value, liste):
    for index, pattern in enumerate(liste):
        if re.search(pattern, value):
            return True

    return False


file = open(genome_fasta, "r")
line = file.readline()
dico_size_genome = {}

while line:
    if line[0] == ">":
        chrom = line[1:].strip()
        line  = file.readline()
        if regex_in_list(chrom, chrom_tab) :
            dico_size_genome[chrom] = len(line)

    line = file.readline()
    

print("size chrom :", dico_size_genome, "\n")


dico_start_chrom = {}
size_genome = 0
for i in dico_size_genome:
    dico_start_chrom[i] = size_genome
    size_genome += dico_size_genome[i]


print("size_genome : " + str(size_genome) + "\n")
print("chrom start :", dico_start_chrom, "\n")

####-----------------------------------------------------

### GET BEST MATCHS GENES


df = pd.read_csv(name_file_gene_bln, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

print("shape total :", df.shape)

name_genome = genome_fasta.split(".")[0]
print("name genome :", name_genome)
os.system("mkdir -p " + name_genome)

#premiere passe 
df_best_matchs   = df
df_for_pic_goupe = df

tab_take_top     = []
tab_take_top_ind = []
for index, value in enumerate(df_best_matchs.values):
    chaine = value[0] #sseqid + qseqid
    
    qseqid_best_top = df_best_matchs["qseqid"].values[index].split(":")
    gene            = qseqid_best_top[0]
    chrom_qseq      = qseqid_best_top[1][1:]
    chrom_sseq      = df_best_matchs["sseqid"].values[index].split("_")[0]

    #best top blast
    if chaine not in tab_take_top :#and chrom_qseq == chrom_sseq:
        tab_take_top.append(chaine)
        tab_take_top_ind.append(index)


df_best_matchs = df_best_matchs.iloc[tab_take_top_ind]
genes          = []
for index, value in enumerate(df_best_matchs.values):
    genes.append(df_best_matchs["qseqid"].values[index].split(":-")[0].lower())
    

df_best_matchs["gene"] = genes
print("shape of best matchs :", df_best_matchs.shape)
print("-------------END GET BEST MATCHS----------------")

####-----------------------------------------------------


##### COMBINE GENES


#Second passe
dico_out = {"qseqid":[], "sseqid":[], "qstart":[], "qend":[], "sstart":[], "send":[], "qsize":[], "ssize":[], "rsize":[]}
dico_out["persize"] = []
dico_out["gene"] = []
dico_out["samechrom"] = []

df = df_for_pic_goupe
df_min = df[df["evalue"] == 0]


def all_same_seq(df):
    all_same = True
    for e, w in enumerate(df.values):
        nb = len(df.values)
        qstart = df["qstart"].values[e]
        qend = df["qend"].values[e]
        if len(df[df["qstart"] == qstart ].values) != nb :
            all_same = False
            
    return all_same

for e, v in enumerate(df_best_matchs.values):
    qseqid = df_best_matchs["qseqid"].values[e]
    sseqid = df_best_matchs["sseqid"].values[e]
    test = df_min[df_min["qseqid"] == qseqid]
    test = test[test["sseqid"] == sseqid]
    #display(test)

    sstart_best = test["sstart"].values[0]
    ssend_best = test["send"].values[0]
    size        = test["qseqid"].values[0].split(":")[-1]

    if int(test["qstart"].values[0]) == 1 and int(test["qend"].values[0]) == int(size) :
        if int(test["sstart"].values[0]) < int(test["send"].values[0]):
            start_min = test["sstart"].values[0]
            send_max  = test["send"].values[0]
        else :
            start_min = test["send"].values[0]
            send_max  = test["sstart"].values[0]
    
    elif int(sstart_best) < int(ssend_best):#forward
        qstart_best = test["qstart"].values[0]
        marge_restant_start = 10000
        restant_start = abs(1 - int(qstart_best))
        test   = test[test["evalue"] <= 0]
        test   = test[test["sstart"] >= int(sstart_best) - int(restant_start) - marge_restant_start]
        test   = test.sort_values(by=["sstart"])
        #display(test)


        marge_send = 10000
        sstart = test["sstart"].values[0]
        test   = test[test["send"] <= int(sstart)+int(size)+marge_send]
        #print(int(sstart)+int(size)+marge_send)
        #display(test)
        all_same = all_same_seq(test)
        if all_same:
            max_score = max(test["bitscore"].values)
            test = test[test["bitscore"] == max_score]
            send_max  = test["send"].values[0]
            start_min = test["sstart"].values[0]
        else:
            send_max  = max(test["send"].values)
            start_min = min(test["sstart"].values)
        #print(start_min, send_max)

    else :#complement
        qstart_best = test["qstart"].values[0]
        marge_restant_start = 10000
        restant_start = abs(1 - int(qstart_best))
        test   = test[test["evalue"] <= 0]

        test   = test[test["send"] >= int(ssend_best) - int(restant_start) - marge_restant_start]
        test   = test.sort_values(by=["send"])
        #display(test)

        marge_send = 10000
        sstart = test["sstart"].values[0]
        test   = test[test["sstart"] <= int(sstart) + int(size) + marge_send]
        qend_max = max(test["qend"].values)
        limit_send = test[test["qend"] == qend_max]["send"].values[-1]
        test = test[test["send"] >= limit_send]

        min_qstart = min(test["qstart"].values)
        limi_max_sstart = test[test["qstart"] == min_qstart]["sstart"].values[-1]
        test = test[test["sstart"] <= limi_max_sstart]
        ##print(int(sstart)+int(size)+marge_send)
        #display(test)

        all_same = all_same = all_same_seq(test)
    
        if all_same:
            max_score = max(test["bitscore"].values)
            test = test[test["bitscore"] == max_score]
            send_max  = test["sstart"].values[0]
            start_min = test["send"].values[0]
        else:
            send_max  = max(test["sstart"].values)
            start_min = min(test["send"].values)
        #print(start_min, send_max)
    
    dico_out["qseqid"].append(test["qseqid"].values[0])
    dico_out["sseqid"].append(test["sseqid"].values[0])
    dico_out["qstart"].append(min(test["qstart"].values))
    dico_out["qend"].append(max(test["qend"].values))
    dico_out["sstart"].append(start_min)
    dico_out["send"].append(send_max)
    dico_out["qsize"].append(abs( max(test["qend"].values) - min(test["qstart"].values) ))
    dico_out["ssize"].append(abs(send_max-start_min))
    dico_out["rsize"].append(int(test["qseqid"].values[0].split(":")[-1]))
    dico_out["persize"].append((dico_out["ssize"][-1]/dico_out["rsize"][-1])*100)
    dico_out["gene"].append(qseqid.split(":-")[0].lower())
    dico_out["samechrom"].append(qseqid.split(":-")[1].split(":")[0] == sseqid.split("_")[0] )#WARING CHROM FOR DROSO : TODO
    
    
df_combine = pd.DataFrame(dico_out)
df_combine = df_combine[df_combine["samechrom"] == True]
print("shape of only samechrom=True :", df_combine.shape)
#df_combine = df_combine[df_combine["persize"] >= 0]

df_combine.to_csv(name_genome + "/" + "Gene_" + name_genome + "_combine.csv", sep=";", index=None)



####-----------------------------------------------------


### ALL TE DATAFRAME

df_for_cluster = pd.read_csv(name_file_all_TE, sep="\t")


####-----------------------------------------------------

### 


file_name_gene_flank = open(gene_flank_file, "r")
couples_genes        = []
for line in file_in:
    genes_fk = line.strip().split("\t")
    couples_genes.append([genes_fk[0].lower(), genes_fk[1].lower()])
    
print("couples genes :", couples_genes[:5],  len(couples_genes))
file_in.close()


####-----------------------------------------------------


df = df_combine
#df = df_best_matchs

dico_new = {'name_clust': [], 'chrom': [], 'start': [], 'stop': [], 'size': [], 'ref_pos_deb': [], 'ref_pos_fin': [], "ref_size": []}

for i, v in enumerate(couples_genes):
    first = v[0]
    secon = v[1]
    
    #start = df[df["gene"] == first]["sstart"].values[0]
    chif_first = re.search("^[0-9]+$", first)
    chif_secon = re.search("^[0-9]+$", secon)
    if len(df[df["gene"] == first].values):
        start = None
        stop  = None
        chrom = None
        qseqid_deb = "None"
        qseqid_fin = "None"
        ref_size   = 0
        

        if secon != "" and secon != first  :
            if not chif_first and not chif_secon and len(df[df["gene"] == secon].values) and len(df[df["gene"] == first].values):
                
                chrom  = df[df["gene"] == first]["sseqid"].values[0]
                chrom2 = df[df["gene"] == secon]["sseqid"].values[0]
                #print(chrom, chrom2)
                qseqid_deb = df[df["gene"] == first]["qseqid"].values[0]
                qseqid_fin = df[df["gene"] == secon]["qseqid"].values[0]
                #print(len(df[df["gene"] == secon].values), qseqid_deb)
                ref_size   = abs(int(qseqid_deb.split(":-")[1].split(":")[1].split("..")[1].replace(")",""))-int(qseqid_fin.split(":-")[1].split(":")[1].split("..")[0].replace("complement(","")))
                chrom_ref  = qseqid_deb.split(":-")[1].split(":")[0]
                chrom_ref2 = qseqid_fin.split(":-")[1].split(":")[0]
                if chrom != chrom2 or chrom_ref != chrom_ref2:
                    print("chrom : ", chrom, chrom2, chrom_ref, chrom_ref2)
                    chrom = None
                start = max(df[df["gene"] == first]["send"].values[0], df[df["gene"] == first]["sstart"].values[0])
                stop  = min(df[df["gene"] == secon]["send"].values[0], df[df["gene"] == secon]["sstart"].values[0])
            elif chif_secon:
                
                qseqid_deb = df[df["gene"] == first]["qseqid"].values[0]
                
                chrom    = df[df["gene"] == first]["sseqid"].values[0]
                start    = max(df[df["gene"] == first]["send"].values[0], df[df["gene"] == first]["sstart"].values[0])
                stop     = int(start) + int(secon)
                stop     = min(int(stop), int(dico_size_genome[chrom]) - 1)
                ref_size = abs(int(stop)-int(start))
                #secon = ref_size
            elif chif_first and len(df[df["gene"] == secon].values):
                qseqid_fin = df[df["gene"] == secon]["qseqid"].values[0]
                chrom = df[df["gene"] == secon]["sseqid"].values[0]
                stop = min(df[df["gene"] == secon]["send"].values[0], df[df["gene"] == secon]["sstart"].values[0])
                start = int(stop)-int(first)
                ref_size = abs(int(stop)-int(start))
                
            if not chrom:
                print(first, secon)
            else:
                dico_new["start"].append(min(int(start), int(stop)))
                dico_new["stop"].append(max(int(start), int(stop)))
                #dico_new["chrom"].append(chrom.split("_")[0])
                dico_new["chrom"].append(chrom)
                dico_new["ref_pos_deb"].append(qseqid_deb)
                dico_new["ref_pos_fin"].append(qseqid_fin)
                dico_new["size"].append(abs(max(int(start), int(stop))-min(int(start), int(stop))))
                dico_new["name_clust"].append("-:-".join([str(first), str(secon)]))
                dico_new["ref_size"].append(ref_size)
                
        elif secon == first:
            chrom = df[df["gene"] == first]["sseqid"].values[0]
            start = df[df["gene"] == first]["sstart"].values[0]
            stop  = df[df["gene"] == secon]["send"].values[0]
            dico_new["start"].append(min(int(start), int(stop)))
            dico_new["stop"].append(max(int(start), int(stop)))
            #dico_new["chrom"].append(chrom.split("_")[0])
            dico_new["chrom"].append(chrom)
            dico_new["size"].append(abs(max(int(start), int(stop))-min(int(start), int(stop))))
            dico_new["name_clust"].append("-:-".join([str(first), str(secon)]))
            dico_new["ref_pos_deb"].append(qseqid_deb)
            dico_new["ref_pos_fin"].append(qseqid_fin)
            dico_new["ref_size"].append(ref_size)


df_pic = pd.DataFrame(data=dico_new)
df_pic = df_pic.drop_duplicates()
print("shape cluster with flanking genes :", df_pic.shape)

df_pic.to_csv(name_genome + "/" + "PiCluster_Gene_" + name_genome + ".csv", sep="\t", index=None)



####-----------------------------------------------------


###


df = df_for_cluster #df TE

index_te_in_cluster     = []
index_te_not_in_cluster = []
tab_cluster_in           = []
for i, v in enumerate(df.values):
    qseqid    = df["qseqid"].values[i]
    chrom     = qseqid.split(":")[0]#.split("_")[0]
    start     = min(int(qseqid.split(":")[2]), int(qseqid.split(":")[3]))
    end       = max(int(qseqid.split(":")[2]), int(qseqid.split(":")[3]))
    in_cluter = False

    df_tmp    = df_pic[df_pic["chrom"] == chrom]

    for e, vc in enumerate(df_tmp.values):
        start_cluster = min(int(df_tmp["start"].values[e]), int(df_tmp["stop"].values[e]))
        end_cluster   = max(int(df_tmp["start"].values[e]), int(df_tmp["stop"].values[e]))

        #TODO
        if (start <= end_cluster and start >= start_cluster) or (end <= end_cluster and end >= start_cluster):
            index_te_in_cluster.append(i)
            tab_cluster_in.append(df_tmp["name_clust"].values[e])
            in_cluter = True
            break# :(

    if in_cluter == False :
        index_te_not_in_cluster.append(i)


df = df.iloc[index_te_in_cluster]
df["name_cluster"] = tab_cluster_in
print("shape of cluster contains TE :", df.shape)
#df = df.iloc[index_te_not_in_cluster]
df_for_count = df


####-----------------------------------------------------

### COUNT FOR GGPLOT

df_cout = df_for_count.groupby(["name_cluster", "sseqid"]).count()


tab_x = []
tab_y = []
tab_z = []
for i, value in enumerate(df_cout_out.values):
    tab_y.append(df_cout.index[i][0])
    tab_x.append(df_cout.index[i][1])
    tab_z.append(value[0])
    
prefixe_count = name_genome + "_IN_CLUSTER_COV"

print("prefix file : ", prefixe_count)
d = {'x': tab_x, 'y': tab_y, 'z': tab_z}

df_out = pd.DataFrame(data=d)
df_out = df_out.sort_values(by="z")
print(df_out.shape)

#df_out.to_csv(name_genome + "/" + str(prefixe_count) + "_hit_map.csv", sep="\t", index=None)
#print("file out >> :", name_genome + "/" + str(prefixe_count) + "_hit_map.csv")
print("Number of TE in cluster :", sum(df_out["z"].values))


## PRINT 
for i, v in enumerate(df_o["x"].unique()):
    print(v, sum(df_o[df_o["x"] == v]["z"]))
    
print(list(df_o["x"].unique()))


tab_add = []
nb_add  = 0
size    = len(df_o.values)
df_out  = pd.DataFrame({"x":df_o["x"].values, "y":df_o["y"].values, "z":df_o["z"].values})

for i, clust in enumerate(df_pic["name_clust"].unique()):
    df_tmp = df_o[df_o["y"] == clust]
    elements = ['gtwin', 'blood', 'ZAM', 'roo', '412']
    for e, elem in enumerate(elements):
        if elem not in list(df_tmp["x"].unique()):
            df_out.loc[size + nb_add] = list([elem, clust, 0])
            nb_add += 1
            
print(len(df_out.values), nb_add)
print(df_pic["name_clust"].unique())

df_out.to_csv(name_genome + "/" + str(prefixe_count) + "_hit_map.csv", sep="\t", index=None)
print("file out >> :", name_genome + "/" + str(prefixe_count) + "_hit_map.csv")


####-----------------------------------------------------
