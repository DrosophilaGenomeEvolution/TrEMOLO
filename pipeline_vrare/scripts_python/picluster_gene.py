import sys
import re
import pandas as pd


file_in = open("pair_genes.csv", "r")
lines = file_in.readlines()
couples_genes = []
for i, v in enumerate(lines):
    couples_genes.append([v.split("\t")[0].strip().lower(), v.split("\t")[1].strip().lower()])
    
#print(couples_genes[:5],  len(couples_genes))
file_in.close()


import pandas as pd
import re

#name_file = "./BLAST/Gene_Ds_SJ_27.bln"
name_file = sys.argv[1]
df = pd.read_csv(name_file, "\t", header=None)
name_keep = name_file.split("/")[-1].split(".")[0]
print(name_keep)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
#display(df.head())
#df = df[["sseqid", "qseqid", "pident", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]]

tab_take_top = []
tab_take_top_ind = []
tab_gene = []
for i, v in enumerate(df.values):
    chaine = v[0]
    sp = df["qseqid"].values[i].split(":")
    gene = sp[0]
    chrom = sp[1][1:]
    chrom2 = df["sseqid"].values[i].split("_")[0]
    #chrom2 = df["sseqid"].values[i].split("_")[1]
    tab_gene.append(gene.strip().lower())
    
    #best top blast
    #if chaine not in tab_take_top and chrom == chrom2:
    if chaine not in tab_take_top :
        tab_take_top.append(chaine)
        tab_take_top_ind.append(i)

df["gene"] = tab_gene
df = df.iloc[tab_take_top_ind]
#display(df.head())
#print(df.shape)
#print(df[df["gene"] == "pld"])
#print(df[df["gene"] == "jing"])

import pandas as pd
import sys
import re

sys.stdout = open(name_keep+"_info.txt", "w")

tab_name_cluster = []
tab_start = []
tab_stop = []
tab_chrom = []
tab_size = []
tab_ref_pos_deb = []
tab_ref_pos_fin = []
tab_ref_size    = []
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
        ref_size = 0
        if secon != "" and secon != first  :
            if not chif_first and not chif_secon and len(df[df["gene"] == secon].values):
                chrom = df[df["gene"] == first]["sseqid"].values[0]
                chrom2 = df[df["gene"] == secon]["sseqid"].values[0]
                qseqid_deb = df[df["gene"] == first]["qseqid"].values[0]
                qseqid_fin = df[df["gene"] == secon]["qseqid"].values[0]
                ref_size = abs(int(qseqid_deb.split(":")[2].split("..")[1].replace(")",""))-int(qseqid_fin.split(":")[2].split("..")[0].replace("complement(","")))
                chrom_ref = qseqid_deb.split(":")[1].replace("-", "")
                chrom_ref2 = qseqid_fin.split(":")[1].replace("-", "")
                if chrom != chrom2 or chrom_ref != chrom_ref2:
                    if chrom_ref == chrom_ref2 :
	                    print("chrom_gene1 :", chrom+"\n", "chrom_gene2 :", chrom2+"\n", "chrom_ref_gene1 :", chrom_ref+"\n", "chrom_ref_gene2 :", chrom_ref2+"\n")
	                    print("gene1 : "+first+"\n", "gene2 : "+secon+"\n---------------------")
                    chrom = None
                start = max(df[df["gene"] == first]["send"].values[0], df[df["gene"] == first]["sstart"].values[0])
                stop  = min(df[df["gene"] == secon]["send"].values[0], df[df["gene"] == secon]["sstart"].values[0])
            elif chif_secon:
                
                qseqid_deb = df[df["gene"] == first]["qseqid"].values[0]
                
                chrom = df[df["gene"] == first]["sseqid"].values[0]
                start = max(df[df["gene"] == first]["send"].values[0], df[df["gene"] == first]["sstart"].values[0])
                stop  = int(start) + int(secon)
                ref_size = abs(int(secon))
            elif chif_first and len(df[df["gene"] == secon].values):
                qseqid_fin = df[df["gene"] == secon]["qseqid"].values[0]
                chrom = df[df["gene"] == secon]["sseqid"].values[0]
                stop = min(df[df["gene"] == secon]["send"].values[0], df[df["gene"] == secon]["sstart"].values[0])
                start = int(stop)-int(first)
                ref_size = abs(int(first))
            if not chrom:
                #print(first, secon)
                print("")
            else:
                tab_start.append(min(int(start), int(stop)))
                tab_stop.append(max(int(start), int(stop)))
                tab_chrom.append(chrom.split("_")[0])
                #tab_chrom.append(chrom)
                tab_ref_pos_deb.append(qseqid_deb)
                tab_ref_pos_fin.append(qseqid_fin)
                tab_size.append(abs(max(int(start), int(stop))-min(int(start), int(stop))))
                tab_name_cluster.append("-:-".join([first, secon]))
                tab_ref_size.append(ref_size)
                
        elif secon == first:
            chrom = df[df["gene"] == first]["sseqid"].values[0]
            start = df[df["gene"] == first]["sstart"].values[0]
            stop  = df[df["gene"] == secon]["send"].values[0]
            tab_start.append(min(int(start), int(stop)))
            tab_stop.append(max(int(start), int(stop)))
            tab_chrom.append(chrom.split("_")[0])
            #tab_chrom.append(chrom)
            tab_size.append(abs(max(int(start), int(stop))-min(int(start), int(stop))))
            tab_name_cluster.append("-:-".join([first, secon]))
            tab_ref_pos_deb.append(qseqid_deb)
            tab_ref_pos_fin.append(qseqid_fin)
            tab_ref_size.append(ref_size)
            #print(first, secon, "hereeeeee")
            

df_pic = pd.DataFrame(data={'name_clust': tab_name_cluster, 'chrom': tab_chrom, 'start':tab_start, 'stop': tab_stop, 'size':tab_size, 'ref_pos_deb':tab_ref_pos_deb, 'ref_pos_fin':tab_ref_pos_fin, "ref_size":tab_ref_size})
#display(df_pic)

df_pic.to_csv("PiCluster_" + name_keep + ".csv", sep="\t", index=None)
#df_pic[df_pic["name_clust"] == "papss-:-su(z)12"]
print("nombre picluster : ", df_pic.shape[0])
sys.stdout.close()