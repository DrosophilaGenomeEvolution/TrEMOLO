import sys
import os
import subprocess



file = open(sys.argv[1], "r")#Empty sites
name_file_tsd = sys.argv[2]#TSD

def calcule_suff(P):
    m = len(P) - 1

    suff = [0] * (m + 1)
    suff[m] = m
    g = m
    f = m

    for i in range(m-1, 0, -1):
        #print("here", i, m)
        if i > g and suff[i+m-f] != i - g:
            suff[i] = min(suff[i + m - f], i-g)
        else :
            f = i
            g = min(g, i)

            while g > 0 and P[g] == P[g + m - f]:
                g = g - 1
            
            suff[i] = f - g

    return suff


def calcule_D(P):

    m = len(P) - 1
    D = [0] * (m + 1)
    suff = calcule_suff(P)

    #print("suff: ", suff[1:])
    #python m+1 = m
    for i in range(1, m + 1):
        D[i] = m

    i = 1
    for j in range(m-1, -1, -1):
        if j == 0 or suff[j] == j:
            while i <= m-j:
                D[i] = m - j
                i = i + 1


    for j in range(1, m):
        D[m-suff[j]] = m - j


    return D


def calcule_R(P):
    m = len(P) - 1

    R = {"A":0, "T":1, "G":2, "C":3}

    for i in range(1, m + 1):
        R[P[i]] = i

    return R

#Boyer_moore
def BM(P, T):
    P = "#" + P
    T = "#" + T

    all_position = []
    
    D = calcule_D(P)
    R = calcule_R(P)

    #print("D:", D[1:])
    #print("R:", R)
    m = len(P) - 1
    n = len(T) - 1

    #print(P, T, m, n)

    i = 0
    pos = 1
    nb_compare = 0

    while  pos <= n - m + 1:
        i = m
        while i>0 and P[i] == T[pos+i-1]:

            nb_compare += 1
            i -= 1

        if i > 0:
            nb_compare += 1

        if i == 0:
            #print(P, "apparait ", pos)
            all_position.append(pos - 1)
            pos = pos + D[1]
            
        else:
            #print("Echec ", i, " motif et ", (pos + i -1), " du texte ", i, P[i])
            pos = pos + max(D[i], i-R[T[pos+i-1]])

    #print("comparaison", nb_compare)
    return all_position


def rev_comp(seq):
    seq_out = ""
    comp = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for i, v in enumerate(seq[::-1]) :
        seq_out += comp[v]

    return seq_out

def grep(motif, file, options=""):
    proc = subprocess.Popen(["grep "+ options + " " + motif + " " + file], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    #print("youhou", str(out).split("\n"))
    return str(out).split("\\n")

#print(grep("KO", "total_results_tsd_ZAM_KO*"))

line = file.readline()

kmer_size = int(sys.argv[3])
total = 0
nb_OK = 0
nb_KO = 0

def search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision, printing=True):

    global nb_OK
    global nb_KO
    global kmer_size
    
    chain_to_show = ""

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

    positions_kl_fl = BM(kmer_l, FK_L)
    positions_kl_fr = BM(kmer_l, FK_R)
    if len(positions_kl_fl) and len(positions_kl_fr) :
        pos_kl_fl = positions_kl_fl[-1]#best positino left
        pos_kl_fr = positions_kl_fr[0]#...right

        dico["kl"]["FKL"] = [FK_L[:pos_kl_fl + kmer_size], FK_L[pos_kl_fl + kmer_size:]]
        dico["kl"]["FKR"] = [FK_R[:pos_kl_fr], FK_R[pos_kl_fr:]]


    positions_kr_fl = BM(kmer_r, FK_L)
    positions_kr_fr = BM(kmer_r, FK_R)
    if len(positions_kr_fl) and len(positions_kr_fr) :
        pos_kr_fl = positions_kr_fl[-1]
        pos_kr_fr = positions_kr_fr[0]

        dico["kr"]["FKL"] = [FK_L[:pos_kr_fl + kmer_size], FK_L[pos_kr_fl + kmer_size:]]
        dico["kr"]["FKR"] = [FK_R[:pos_kr_fr], FK_R[pos_kr_fr:]]

    # FASTA/reads_2L:18734:21516216-21524592.fasta
    # ----EMPTY SITE----
    # >2L_RaGOO_RaGOO:21516216
    # CAGAGCAAAGGAGGTTGGTAGGCAG+CGCG+C-:-GAGCCATTTTTAACAGAAAAAAGTGTTCTC
    # GAGAACACTTTTTTCTGTTAAAAATGGCTC-:-GCGCGCTGCCTACCAACCTCCTTTGCTCTG
    # ----END EMPTY SITE----
    # >2L:<INS>:21516216:21524592:18734:1:PRECISE:2-8381
    # (GCGCGA, GCGCGA, [KO:18734], 23, 1, 6)

    left_or_right = "k"
    if len(dico["kl"]) == 2:
        left_or_right += "l"

    elif len(dico["kr"]) == 2:
        left_or_right += "r"

    else:
        if precision == "PRECISE":#PRECISE

            chain_to_show += tab_old_line[0] + "\n"
            chain_to_show += tab_old_line[1] + "\n"
            chain_to_show += tab_old_line[2] + "\n"
            chain_to_show += tab_old_line[3] + "\n"
            chain_to_show += tab_old_line[4] + "\n"

            nb_KO += 1
        else :
            if printing:
                return False
            else:
                return False, ""
            #REcursion
            # empty_site_seq_l
            # boole = search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision)
            # if len(empty_site_seq_r) == kmer_size or boole :
            #   return True

    if left_or_right != "k":
        # print(dico)
        # print("kmer:", kmer_l, kmer_r)
        # print("fank:", FK_L, FK_R)

        # print("BM KL")
        # print(BM(kmer_l, FK_L))
        # print(BM(kmer_l, FK_R))

        # print("BM FR")
        # print(BM(kmer_r, FK_L))
        # print(BM(kmer_r, FK_R))
        # print("------------")

        #print("---------", dico["kl"])
        #redefine TE
        ETC  = dico[left_or_right]["FKL"][1] + ETC
        ETC += dico[left_or_right]["FKR"][0]

        #redefine FLANK
        FK_L = dico[left_or_right]["FKL"][0]
        FK_R = dico[left_or_right]["FKR"][1]

        

        seq_and_tsd = FK_L[:len(FK_L) - kmer_size] + "++:" + FK_L[len(FK_L) - kmer_size:] + ":++" + "--|" + ETC.strip() + "|--" +  "++:" + FK_R[:kmer_size] + ":++" + FK_R[kmer_size:] 
        #print("FK_L", FK_L, FK_R)
        # print(header_1)
        # print(header_2)
        # print(header_3)
        # print("TSM=" + FK_L[len(FK_L) - kmer_size:], "[KO->OK:" + str(ID)+ "]", kmer_size)
        # print(seq_and_tsd)

        chain_to_show += header_1 + "\n"
        chain_to_show += header_2 + "\n"
        chain_to_show += header_3 + "\n"
        if len(FK_L[len(FK_L) - (kmer_size * 2) + 2:len(FK_L) - (kmer_size)] + FK_R[:kmer_size + 2]) == kmer_size + 4:
            chain_to_show += " ".join(["SVI=" + FK_L[len(FK_L) - (kmer_size * 2) + 2:len(FK_L) - (kmer_size)] + FK_R[:kmer_size + 2], "[KO->OK:" + str(ID)+ "]", str(kmer_size) + "\n"])
        chain_to_show += " ".join(["TSM=" + FK_L[len(FK_L) - kmer_size:], "[KO->OK:" + str(ID)+ "]", str(kmer_size) + "\n"])
        chain_to_show += seq_and_tsd + "\n"

        if printing:
            nb_OK += 1
        
    if printing:
        print(chain_to_show)
        return len(dico["kr"]) == 2 or len(dico["kl"]) == 2
    else:
        return [len(dico["kr"]) == 2 or len(dico["kl"]) == 2, chain_to_show]


def find_decale_svi(empty_site, start, stop):
    itera = 1
    if start > stop:
        itera = -1

    pos = {"text":"None", "find":False}
    for i in range(start, stop, itera):
        empty_site_seq_l_slide = empty_site[:i]
        empty_site_seq_r_slide = empty_site[i:]
        find, chaine = search(ID, empty_site_seq_l_slide, empty_site_seq_r_slide, tab_old_line, name_file_tsd, "IMPRECIS", False)
        if find :
            cr = ""
            cr += " ".join(["ID: ", str(ID) + "\n"])
            cr += ">" + chrom + ":" + str(position) + "\n"
            cr += ", ".join(["EMPTY_LEFT=" + empty_site_seq_l, "EMPTY_RIGTH=" + empty_site_seq_r + "\n"])
            cr += ", ".join(["kmer_LEFT="+kmer_l, "kmer_RIGTH=" + kmer_r + "\n"])
            cr += ", ".join([empty_site_seq_l, empty_site_seq_r + "\n"])
            cr += kmer_l + ", " + kmer_r + "\n"
            cr += "##determining new position of empty site" + "\n"
            cr += ", ".join(["POS_EMPTY="+str(i), empty_site_seq_l_slide, empty_site_seq_r_slide + "\n"])
            cr += ">" + chrom + ":" + str(position) + "\n"

            chaine = cr + chaine
            pos = {"pos":i, "text":chaine, "find":find}
            return pos

        elif i + itera == stop:

            # print("")
            # print("ID: ", ID)
            # print(">" + chrom + ":" + str(position))
            # print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
            # print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
            
            # print(empty_site_seq_l, empty_site_seq_r)
            # print(kmer_l, kmer_r)

            # #print(precision)
            # print(tab_old_line[0])
            # print(tab_old_line[1])
            # print(tab_old_line[2])
            # print(tab_old_line[3])
            # print(tab_old_line[4])

            cr = ""
            cr += " ".join(["ID: ", str(ID) + "\n"])
            cr += ">" + chrom + ":" + str(position) + "\n"
            cr += ", ".join(["EMPTY_LEFT=" + empty_site_seq_l, "EMPTY_RIGTH=" + empty_site_seq_r + "\n"])
            cr += ", ".join(["kmer_LEFT="+kmer_l, "kmer_RIGTH=" + kmer_r + "\n"])
            cr += ", ".join([empty_site_seq_l, empty_site_seq_r + "\n"])
            cr += kmer_l + ", " + kmer_r + "\n"

            cr += tab_old_line[0] + "\n"
            cr += tab_old_line[1] + "\n"
            cr += tab_old_line[2] + "\n"
            cr += tab_old_line[3] + "\n"
            cr += tab_old_line[4] + "\n"

            pos = {"pos":i, "text":cr, "find":False}

    return pos


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
            #   print(precision)
            #search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision)
            find = search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, "IMPRECIS")
            if not find:
                print("precise not find")
                origine_position_sniffle = len(empty_site_seq_l)
                #empty_site = empty_site_seq_l[origine_position_sniffle - (kmer_size*2):] + empty_site_seq_r[:kmer_size*2] ##all region
               # empty_site = empty_site_seq_l + empty_site_seq_r ##all region
                #origine_position_sniffle = len(empty_site_seq_l[origine_position_sniffle - (kmer_size*2):])
                #print("empty --", empty_site)
                #pos_l = find_decale_svi(empty_site, origine_position_sniffle, origine_position_sniffle - (kmer_size * 2))

                #pos_r = find_decale_svi(empty_site, origine_position_sniffle, origine_position_sniffle + kmer_size * 2)

                empty_site = empty_site_seq_l + empty_site_seq_r ##all region
                origine_position_sniffle = len(empty_site_seq_l)
                pos_l = find_decale_svi(empty_site, origine_position_sniffle, kmer_size - 1)

                pos_r = find_decale_svi(empty_site, origine_position_sniffle, len(empty_site) - kmer_size)

                if pos_r["find"] and pos_l["find"] and abs(pos_r["pos"] - origine_position_sniffle) < abs(pos_l["pos"] - origine_position_sniffle) :
                    print(pos_r["text"])
                    nb_OK += 1

                elif pos_l["find"]:
                    print(pos_l["text"])
                    nb_OK += 1

                elif pos_r["find"]:
                    print(pos_r["text"])
                    nb_OK += 1
                else:
                    print(False)
                    print(pos_r["text"])
                    nb_KO += 1

        else:
            print("")
            print("ID: ", ID)
            print(">" + chrom + ":" + str(position))
            print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
            print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
            print(empty_site_seq_l, empty_site_seq_r)
            print(kmer_l, kmer_r)

            pos_l, pos_r = [{"find":False}, {"find":False}]
            if not search(ID, empty_site_seq_l, empty_site_seq_r, tab_old_line, name_file_tsd, precision):

                empty_site = empty_site_seq_l + empty_site_seq_r ##all region
                origine_position_sniffle = len(empty_site_seq_l)
                pos_l = find_decale_svi(empty_site, origine_position_sniffle, kmer_size - 1)

                pos_r = find_decale_svi(empty_site, origine_position_sniffle, len(empty_site) - kmer_size)

                #position more wear
                if pos_r["find"] and pos_l["find"] and abs(pos_r["pos"] - origine_position_sniffle) < abs(pos_l["pos"] - origine_position_sniffle) :
                    print(pos_r["text"])
                    nb_OK += 1

                elif pos_l["find"]:
                    print(pos_l["text"])
                    nb_OK += 1

                elif pos_r["find"]:
                    print(pos_r["text"])
                    nb_OK += 1
                else:
                    print(False)
                    print(pos_r["text"])
                    nb_KO += 1

                # for i in range(kmer_size, len(empty_site) - kmer_size ):
                #   empty_site_seq_l_slide = empty_site[:i]
                #   empty_site_seq_r_slide = empty_site[i:]
                #   boole = search(ID, empty_site_seq_l_slide, empty_site_seq_r_slide, tab_old_line, name_file_tsd, precision)
                #   if boole :
                #       print("")
                #       print("ID: ", ID)
                #       print(">" + chrom + ":" + str(position))
                #       print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
                #       print("kmer_LEFT="+kmer_l, "kmer_RIGTH=" + kmer_r)
                        
                #       print(empty_site_seq_l, empty_site_seq_r)
                #       print(kmer_l, kmer_r)
                #       print("##determining new position of empty site")
                #       print("POS_EMPTY="+str(i), empty_site_seq_l_slide, empty_site_seq_r_slide)
                #       break
                #   elif len(empty_site_seq_r_slide) == kmer_size+1:

                #       print("")
                #       print("ID: ", ID)
                #       print(">" + chrom + ":" + str(position))
                #       print("EMPTY_LEFT="+empty_site_seq_l, "EMPTY_RIGTH="+empty_site_seq_r)
                #       print("kmer_LEFT="+kmer_l, "kmer_RIGTH="+kmer_r)
                        
                #       print(empty_site_seq_l, empty_site_seq_r)
                #       print(kmer_l, kmer_r)

                #       #print(precision)
                #       print(tab_old_line[0])
                #       print(tab_old_line[1])
                #       print(tab_old_line[2])
                #       print(tab_old_line[3])
                #       print(tab_old_line[4])

        total += 1
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


