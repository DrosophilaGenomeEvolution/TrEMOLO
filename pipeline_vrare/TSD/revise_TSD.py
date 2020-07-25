import sys
import os
import subprocess



file = open(sys.argv[1], "r")
name_file_tsd = sys.argv[2]

def rev_comp(seq):
	seq_out = ""
	comp = {"A":"T", "T":"A", "G":"C", "C":"G"}
	for i, v in enumerate(seq[::-1]) :
		seq_out += comp[v]

	return seq_out

def grep(motif, file, options=""):
	proc = subprocess.Popen(["grep "+ options + " " + motif + " " + file], stdout=subprocess.PIPE, shell=True)
	(out, err) = proc.communicate()
	return out.split("\n")

#print(grep("KO", "total_results_tsd_ZAM_KO*"))

line = file.readline()

kmer_size = int(sys.argv[3])
total = 0
nb_OK = 0
nb_KO = 0

def search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision):

	global nb_OK
	global nb_KO
	global kmer_size

	#empty_site_seq_l_revcomp = rev_comp(empty_site_seq_r)
	#empty_site_seq_r_revcomp = rev_comp(empty_site_seq_l)

	kmer_l = empty_site_seq_l.strip()[-kmer_size:]
	kmer_r = empty_site_seq_r.strip()[:kmer_size]

	#kmer_l_revcomp = empty_site_seq_l_revcomp.strip()[-kmer_size:]
	#kmer_r_revcomp = empty_site_seq_r_revcomp.strip()[:kmer_size]

	# print("")
	# print("ID: ", ID)
	# print(">" + chrom + ":" + str(position))
	# print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
	# print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
	
	# print(empty_site_seq_l, empty_site_seq_r)
	# print(kmer_l, kmer_r)

	#print(empty_site_seq_l_revcomp, empty_site_seq_r_revcomp)
	#print(kmer_l_revcomp, kmer_r_revcomp)

	#GET SEQ TE candidate and FLANK TE in TSD results
	#tab_old_line   = grep("KO:" + str(ID), name_file_tsd, "-A 2 -B 2")
	header_1  = tab_old_line[0].strip()
	header_2  = tab_old_line[1].strip()
	header_3  = tab_old_line[2].strip()
	flank     = tab_old_line[-3]
	sequence  = tab_old_line[-2]

	ETC  = sequence.split("|")[1]
	FK_L = flank.split(":")[0][1:]
	FK_R = flank.split(":")[1]

	new_pos_fkr  = 0
	rest_pos_fkr = 0
	new_pos_fkl  = 0
	rest_pos_fkl = 0

	dico = {"kl": {}, "kr": {}}
	for i in range(len(FK_L)-kmer_size):
		if kmer_r == FK_L[i:i+kmer_size]:
			#print("match FKL kr :", i)
			dico["kr"]["FKL"] = [FK_L[:i+kmer_size], FK_L[i+kmer_size:]]

		if kmer_l == FK_L[i:i+kmer_size]:
			#print("match FKL kl :", i)
			dico["kl"]["FKL"] = [FK_L[:i+kmer_size], FK_L[i+kmer_size:]]

		if kmer_r == FK_R[i:i+kmer_size]:
			#print("match FKR kr :", i)
			dico["kr"]["FKR"] = [FK_R[:i], FK_R[i:]]

		if kmer_l == FK_R[i:i+kmer_size]:
			#print("match FKR kl :", i)
			dico["kl"]["FKR"] = [FK_R[:i], FK_R[i:]]
			new_pos_fkr = i

	# FASTA/reads_2L:18734:21516216-21524592.fasta
	# ----EMPTY SITE----
	# >2L_RaGOO_RaGOO:21516216
	# CAGAGCAAAGGAGGTTGGTAGGCAG+CGCG+C-:-GAGCCATTTTTAACAGAAAAAAGTGTTCTC
	# GAGAACACTTTTTTCTGTTAAAAATGGCTC-:-GCGCGCTGCCTACCAACCTCCTTTGCTCTG
	# ----END EMPTY SITE----
	# >2L:<INS>:21516216:21524592:18734:1:PRECISE:2-8381
	# (GCGCGA, GCGCGA, [KO:18734], 23, 1, 6)

	if len(dico["kl"]) == 2:
		#print("---------", dico["kl"])
		ETC  = dico["kl"]["FKL"][1] + ETC
		ETC += dico["kl"]["FKR"][0]

		FK_L = dico["kl"]["FKL"][0]
		FK_R = dico["kl"]["FKR"][1]

		seq_and_tsd = FK_L[:len(FK_L) - kmer_size] + "++:" + FK_L[len(FK_L) - kmer_size:] + ":++" + "--|" + ETC.strip() + "|--" +  "++:" + FK_R[:kmer_size] + ":++" + FK_R[kmer_size:] 
		#print("FK_L", FK_L, FK_R)
		print(header_1)
		print(header_2)
		print(header_3)
		print("TSM=" + FK_L[len(FK_L) - kmer_size:], "[KO->OK:" + str(ID)+ "]", kmer_size)
		print(seq_and_tsd)
		nb_OK += 1

	elif len(dico["kr"]) == 2:
		#print("---------", dico["kr"])
		ETC  = dico["kr"]["FKL"][1] + ETC
		ETC += dico["kr"]["FKR"][0]

		FK_L = dico["kr"]["FKL"][0]
		FK_R = dico["kr"]["FKR"][1]

		seq_and_tsd = FK_L[:len(FK_L) - kmer_size] + "++:" + FK_L[len(FK_L) - kmer_size:] + ":++" + "--|" + ETC.strip() + "|--" +  "++:" + FK_R[:kmer_size] + ":++" + FK_R[kmer_size:] 
		#print("FK_L", FK_L, FK_R)
		print(header_1)
		print(header_2)
		print(header_3)
		print("TSM=" + FK_L[len(FK_L) - kmer_size:], "[KO->OK:" + str(ID)+ "]", kmer_size)
		print(seq_and_tsd)

		nb_OK += 1

	else:
		if precision == "PRECISE":#PRECISE
			#print(precision)
			print(tab_old_line[0])
			print(tab_old_line[1])
			print(tab_old_line[2])
			print(tab_old_line[3])
			print(tab_old_line[4])
			nb_KO += 1
		else :
			return False
			#REcursion
			# empty_site_seq_l
			# boole = search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision)
			# if len(empty_site_seq_r) == kmer_size or boole :
			# 	return True
			
	return len(dico["kr"]) == 2 or len(dico["kl"]) == 2




ID = ""
while line :

	if line[0] == ">" :
		ID       = line[1:].strip().split(":")[0]
		chrom    = line.strip().split(":")[3]
		position = line.strip().split(":")[4].split("-")[1] 

		empty_site_seq_l = file.readline().strip()
		file.readline()
		empty_site_seq_r = file.readline().strip()

		tab_old_line   = grep("KO:" + str(ID), name_file_tsd, "-A 2 -B 2")
		#header_1  = tab_old_line[0].strip()
		header_2  = tab_old_line[1].strip()
		#header_3  = tab_old_line[2].strip()

		precision = header_2.split(":")[6]

		kmer_l = empty_site_seq_l.strip()[-kmer_size:]
		kmer_r = empty_site_seq_r.strip()[:kmer_size]

		



		#search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, "PRECISE")
		if precision == "PRECISE":#PRECISE
			print("")
			print("ID: ", ID)
			print(">" + chrom + ":" + str(position))
			print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
			print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
			
			print(empty_site_seq_l, empty_site_seq_r)
			print(kmer_l, kmer_r)
			# 	print(precision)
			search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision)
		else:
			if not search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision):
				empty_site = empty_site_seq_l + empty_site_seq_r
				for i in range(kmer_size, len(empty_site)-kmer_size ):
					empty_site_seq_l_slide = empty_site[:i]
					empty_site_seq_r_slide = empty_site[i:]
					boole = search(ID, empty_site_seq_l_slide, empty_site_seq_r_slide, tab_old_line, name_file_tsd, precision)
					if boole :
						print("")
						print("ID: ", ID)
						print(">" + chrom + ":" + str(position))
						print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
						print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
						
						print(empty_site_seq_l, empty_site_seq_r)
						print(kmer_l, kmer_r)
						print("##determining new position of empty site")
						print("POS_EMPTY="+str(i), empty_site_seq_l_slide, empty_site_seq_r_slide)
						break
					elif len(empty_site_seq_r_slide) == kmer_size+1:
						#print(precision)
						print(tab_old_line[0])
						print(tab_old_line[1])
						print(tab_old_line[2])
						print(tab_old_line[3])
						print(tab_old_line[4])

		total += 1
		# else :#IMPRECIS determined position of empty site
		# 	print(precision)
		# 	#for i in range():
		# 		#search(ID, empty_site_seq_l, empty_site_seq_r, name_file_tsd)

		#print(FK_G, FK_R)

		#print(grep("KO:" + str(ID), name_file_tsd, "-A 2")[-2])

	line = file.readline()


#OK/total : 25/51
# KO/total : 23/51
# OK+KO/total : 48/51
# K-O/total : 3/51
# OK+K-O/total : 28/51
# OK% : 49%

print("\n\n")
print("OK/total :" + str(nb_OK) + "/" + str(total))
print("KO/total :" + str(nb_KO) + "/" + str(total))




#proc = subprocess.Popen(["grep KO total_results_tsd_ZAM_KO.txt"], stdout=subprocess.PIPE, shell=True)
#(out, err) = proc.communicate()
#print("program output:", str(out.split("\n")), err)


