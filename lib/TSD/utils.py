import subprocess

#Boyer moore
def calcule_suff(P):
    m = len(P) - 1

    suff = [0] * (m + 1)
    suff[m] = m
    g = m
    f = m

    for i in range(m-1, 0, -1):
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

    m = len(P) - 1
    n = len(T) - 1

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
            all_position.append(pos - 1)
            pos = pos + D[1]
            
        else:
            pos = pos + max(D[i], i-R[T[pos+i-1]])

    return all_position

#END Boyer moore

def rev_comp(seq):
    seq_out = ""
    comp = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
    for i, v in enumerate(seq[::-1]) :
        if v in comp:
            seq_out += comp[v]
        else:
            seq_out += comp["N"]

    return seq_out

def grep(motif, file, options=""):
    proc = subprocess.Popen(["grep "+ options + " " + motif + " " + file], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    return str(out).split("\\n")




#find the real position of SVI
def find_svi_position(TSD, empty_site, pos_svi_init, position_on_genome):
    positions = BM(TSD, empty_site)
    
    dist_mini = len(empty_site)
    decalage  = 0
    new_pos   = position_on_genome

    for index, pos in enumerate(positions) :
        if  dist_mini > abs(pos - pos_svi_init + len(TSD)) :
            dist_mini = abs(pos - pos_svi_init + len(TSD))
            decalage  = pos - pos_svi_init

    if len(positions) :
        new_pos = position_on_genome + decalage + len(TSD)

    return new_pos

