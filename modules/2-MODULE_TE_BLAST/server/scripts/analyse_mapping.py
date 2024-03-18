import sys
import os
import pandas as pd
import json
import numpy as np
import argparse


parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("blast_file", metavar='<blast-file>', option_strings=['blast-file'], type=argparse.FileType('r'),
                    help="Input blast file format 6 (-outfmt 6)")

parser.add_argument("db_file_te", metavar='<db-file-TE>', type=argparse.FileType('r'),
                    help="Multi fasta file TE for get size")

parser.add_argument("query_file_te", metavar='<query-file-TE>', type=argparse.FileType('r'),
                    help="Multi fasta file TE for get size")

parser.add_argument("output_json", metavar='<output-json>', type=argparse.FileType('w'),
                    help="output-json file")

#OPTION
parser.add_argument("-p", "--min-pident", dest='min_pident', type=int, default=94,
                    help="minimum percentage of identity")
parser.add_argument("-s", "--min-size-percent", dest="min_size_percent", type=int, default=90,
                    help="minimum percentage of TE size")
parser.add_argument("-r", "--min-read-support", dest="min_read_support", type=int, default=1,
                    help="minimum read support number")
parser.add_argument("-k", "--type-sv-keep", dest="type_sv_keep", type=str, default="INS",
                    help="Type of SV")
parser.add_argument("-c", "--combine", dest="combine", default=False, action='store_true',
                    help="Combine parts blast TE")
parser.add_argument("--combine_name", dest="combine_name", default="COMBINE_TE.csv",
                    help="Combine name output file")
parser.add_argument("--nb_elements", dest="nb_elements", type=int, default=None,
                    help="nombre d'element dans la data list")
args = parser.parse_args()


name_file        = args.blast_file.name
name_querty_file = args.query_file_te.name
name_output_file = args.output_json.name

nb_elements      = args.nb_elements
min_pident       = args.min_pident
min_size_percent = args.min_size_percent
min_read_support = args.min_read_support
type_sv_keep     = args.type_sv_keep
combine          = args.combine
combine_name     = args.combine_name

#print("[" + str(sys.argv[0]) + "] : PREFIX OUTPUT :", name_out)

 
# 1.   qseqid      query or source (e.g., gene) sequence id
# 2.   sseqid      subject  or target (e.g., reference genome) sequence id
# 3.   pident      percentage of identical matches
# 4.   length      alignment length (sequence overlap)
# 5.   mismatch    number of mismatches
# 6.   gapopen     number of gap openings
# 7.   qstart      start of alignment in query
# 8.   qend        end of alignment in query
# 9.   sstart      start of alignment in subject
# 10.  send        end of alignment in subject
# 11.  evalue      expect value
# 12.  bitscore    bit score

#awk 'OFS="\t"{if(substr($0, 1, 1) == ">"){id=$0}else if(length($0) != 0){print id, length($0)}}' SEQUENCE_INDEL.fasta | sed 's/^>//g' > SIZE_QUERY.txt

size_et = {}

if args.db_file_te != None:
    #GET_SIZE TE Database format fasta
    file  = args.db_file_te
    lines = file.readlines()

    for i, l in enumerate(lines):
        if l[0] == ">":
            size_et[l[1:].strip()] = len(lines[i + 1].strip())

    file.close()

else :
    print("[" + str(sys.argv[0]) + "] : ERROR fasta file Database TE Not Found")
    exit(1)




size_query = {}

if args.query_file_te != None:
    #GET_SIZE TE Database format fasta
    file  = args.query_file_te
    lines = file.readlines()

    for i, l in enumerate(lines):
        if l[0] == ">":
            size_query[l[1:].strip()] = len(lines[i + 1].strip())

    file.close()

else :
    print("[" + str(sys.argv[0]) + "] : ERROR fasta file Database TE Not Found")
    exit(1)


df = pd.read_csv(name_file, sep="\t", header=None)
df.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

#KEEP ONLY TE on the LIST (GET_SIZE CANNONICAL TE)
df = df[df["sseqid"].isin(size_et.keys())]

df.index = [int(i) for i in range(0, len(df.values))]

#print(df["qseqid"].unique())

dico = {"data": []}
count_1 = 0
qseqid_pred = ""
sseqid_pred = ""
SIZE_HTML   = 500
size_html   = SIZE_HTML
array_TE    = []
array_Chrom = []
for index, value in enumerate(df.values):
	qseqid = df["qseqid"].values[index]
	sseqid = df["sseqid"].values[index]
	chrom  = qseqid.split(":")[0]
	start  = int(qseqid.split(":")[2])
	end    = int(qseqid.split(":")[3])

	if chrom not in array_Chrom :
		array_Chrom.append(chrom)

	if sseqid not in array_TE :
		array_TE.append(sseqid)

	if sseqid != sseqid_pred or qseqid != qseqid_pred :
		dico["data"].append({"qseqid":qseqid, "sseqid":sseqid, "chrom":chrom, "start":start})
		qseqid_pred = qseqid
		sseqid_pred = sseqid

	size_cumul = [0, 0]
	strand     = 1
	if df["sstart"].values[index] > df["send"].values[index] :
		strand = -1

	sstart = min(df["sstart"].values[index], df["send"].values[index])
	send   = max(df["sstart"].values[index], df["send"].values[index])
	
	qstart = min(df["qstart"].values[index], df["qend"].values[index])
	qend   = max(df["qstart"].values[index], df["qend"].values[index])

	bitscore = df["bitscore"].values[index]
	evalue   = df["evalue"].values[index]
	gapopen  = df["gapopen"].values[index]
	mismatch = df["mismatch"].values[index]
	pident   = df["pident"].values[index]

	subject_seq_size = size_et[sseqid]
	query_seq_size   = size_query[qseqid]

	# 
	size_html        = (subject_seq_size/query_seq_size)*SIZE_HTML

	if "positions" not in dico["data"][len(dico["data"])-1] :
		dico["data"][len(dico["data"])-1]["positions"] = [{"sstart":sstart, "send":send, "subject_seq_size":subject_seq_size, "query_seq_size":query_seq_size, "strand":strand, "qstart": qstart, "qend": qend, "mismatch": mismatch, "gapopen":gapopen, "evalue": evalue, "bitscore":bitscore, "pident":pident}]
	else :
		dico["data"][len(dico["data"])-1]["positions"].append({"sstart":sstart, "send":send, "subject_seq_size":subject_seq_size, "query_seq_size":query_seq_size, "strand":strand, "qstart": qstart, "qend": qend, "mismatch": mismatch, "gapopen":gapopen, "evalue": evalue, "bitscore":bitscore, "pident":pident})

	index_tab = len(dico["data"][len(dico["data"])-1]["positions"])-1
		
	#empty deb
	if strand == 1:
		percent = round((sstart-1)/size_et[sseqid], 2)	#percent empty deb
	else:
		percent = round((size_et[sseqid]-send)/size_et[sseqid], 2)	#percent empty deb

	empty_size = [int(SIZE_HTML * percent), int(size_html * percent)]

	#
	if empty_size[0] > 0:
		dico["data"][len(dico["data"])-1]["positions"][index_tab]["no_match_start"] = 0
		dico["data"][len(dico["data"])-1]["positions"][index_tab]["no_match_send"] = sstart

		dico["data"][len(dico["data"])-1]["positions"][index_tab]["no_match_size"] = [empty_size[0], empty_size[1]]
		size_cumul[0] += empty_size[0]
		size_cumul[1] += empty_size[1]


	#qpercent empty
	qpercent_empty = round((qstart-1)/size_query[qseqid], 2)
	posq           = [-empty_size[0]+int(SIZE_HTML * qpercent_empty), -empty_size[1]+int(SIZE_HTML * qpercent_empty)]
	q_empty_size   = int(SIZE_HTML * qpercent_empty)

	#qpercent full
	qsize      = abs(qstart - qend)
	qpercent   = round(qsize/size_query[qseqid], 2)
	qsize_html = int(SIZE_HTML * qpercent)

	#TE
	size    = abs(send - sstart)
	percent = round(size/size_et[sseqid], 2)

	#print("q-:", qsize, qpercent, posq, int(size_html * percent), posq-int(size_html * percent)+size_html, empty_size, size_et[sseqid], size_query[qseqid])

	dico["data"][len(dico["data"])-1]["positions"][index_tab]["posq"] = [posq[0], posq[1]] 


	dico["data"][len(dico["data"])-1]["positions"][index_tab]["size"] = [int(SIZE_HTML * percent), int(size_html * percent)]
	size_cumul[0] += int(SIZE_HTML * percent)
	size_cumul[1] += int(size_html * percent)

	if abs(size_cumul[0] - SIZE_HTML) >= 4:
		dico["data"][len(dico["data"])-1]["positions"][index_tab]["rest_size"] = [abs(size_cumul[0] - SIZE_HTML), abs(size_cumul[1] - size_html)]


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)




if nb_elements == None :
	output_file = open(str(name_output_file) + ".json", "w")

	print(len(dico["data"]))
	##to json
	json_dico = json.dumps(dico, cls=NpEncoder)

	output_file.write(json_dico)
	output_file.close()
else :
	size_data = len(dico["data"])
	print(size_data, nb_elements, type(size_data), type(nb_elements), int((size_data/nb_elements)-1))
	limit = int((size_data/nb_elements))
	for i in range(limit-1) :
		output_file = open(str(name_output_file) + "_part" + str(i) + ".json", "w")

		json_dico = json.dumps({"data": dico["data"][i*nb_elements:(i+1)*nb_elements]}, cls=NpEncoder)

		output_file.write(json_dico)
		output_file.close()

#print("var data = " + json_dico)
#print("var array_TE = " + str(array_TE) )
#print("var array_Chrom = " + str(array_Chrom) )

##to html
# for dic in dico["data"] :
# 	#print(dic, "\n")
# 	for pos in dic["positions"] :
# 		#print("pos", pos)
# 		print('<div class="leftGr grpos">')
# 		if "no_match_send" in pos :
# 			print('\t<div style="width: ' + str(pos["no_match_size"]-4) + 'px; padding-left: ' + str(pos["no_match_size"]-4) + 'px; margin-left: ' + str(pos["posq"]) + 'px;" class="h4 white"></div>')
# 			print('\t<div class="red"> <a href="#" qseqid="' + dic["qseqid"] + '" sseqid="' + dic["sseqid"] + '" strand="' + str(pos["strand"]) + '" style="width: ' + str(pos["size"]) + 'px; "></a> </div>')
# 		else :
# 			if "posq" in pos :
# 				print('\t<div class="red" style="margin-left: ' + str(pos["posq"]) + 'px;"> <a href="#" qseqid="' + dic["qseqid"] + '" sseqid="' + dic["sseqid"] + '" strand="' + str(pos["strand"]) + '" style="width: ' + str(pos["size"]) + 'px; "></a> </div>')
# 			else :
# 				print('\t<div class="red"> <a href="#" qseqid="' + dic["qseqid"] + '" sseqid="' + dic["sseqid"] + '" strand="' + str(pos["strand"]) + '" style="width: ' + str(pos["size"]) + 'px;"></a> </div>')
		
# 		if "rest_size" in pos :
# 			print('\t<div style="width: ' + str(pos["rest_size"]-4) + 'px; padding-left: ' + str(pos["rest_size"]-4) + 'px;" class="h4 white"></div>')
# 		print('</div>')