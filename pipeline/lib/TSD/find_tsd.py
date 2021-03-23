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

"""
    GET REGIONy
    ==========================
    :author:  Mourdas MOHAMED 
    :contact: mourdas.mohamed@igh.cnrs.fr
    :date: 01/06/2020
    :version: 0.1
    Script description
    ------------------
    find_tsd.py find TSD of TE
    -------
    >>> find_tsd.py ...
    
    Help Programm
    -------------
    usage: find_tsd.py [-h] flank_fasta size_flank id_elem strands size_tsd

    positional arguments:
      flank_fasta  fasta file of flank
      size_flank   size of flankan
      id_elem      id of variant rare
      strands      strands of sequence
      size_tsd     name fasta file out

    optional arguments:
      -h, --help   show this help message and exit
"""



import os
import sys
import re
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import argparse

parser = argparse.ArgumentParser()

#MAIN ARGS
parser.add_argument("flank_fasta", type=str,
                    help="fasta file of flank")
parser.add_argument("file_seq_TE", type=str,
                    help="file of sequence TE")
parser.add_argument("size_flank", type=int,
                    help="size of flankan")


parser.add_argument("id_elem", type=str,
                    help="id of variant rare")
parser.add_argument("strands", type=str,
                    help="strands of sequence")
parser.add_argument("size_tsd", type=int,
                    help="name fasta file out")

args  = parser.parse_args()

file  = open(args.flank_fasta, "r")
lines = file.readlines()

id_elem          = args.id_elem
file_sequence_TE = open(args.file_seq_TE, "r")
seq_lines        = file_sequence_TE.readlines()
head_seq         = seq_lines[0].strip()
TE               = seq_lines[1].strip()

def rev_comp(seq):
    seq_out = ""
    comp = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for i, v in enumerate(seq[::-1]) :
        seq_out += comp[v]

    return seq_out


#gere le forward et reverse
strands = False #False + et True -
if args.strands == "-":
    strands = True
    sequence_1 = lines[3].strip()
    sequence_2 = lines[1].strip()
    sequence_1 = rev_comp(sequence_1)
    sequence_2 = rev_comp(sequence_2)
    TE = rev_comp(TE)
else :
    sequence_1 = lines[1].strip()
    sequence_2 = lines[3].strip()


def align(alignments, size):
    liste = None
    i = 0
    while i < len(alignments) and liste == None:
        if alignments[i][-1] == size :
            liste = list(alignments[i])
        i += 1  

    return liste, (i-1)


size_k       = int(args.size_flank)
size_tsd     = int(args.size_tsd)
size_flank   = size_k
find         = False
report_seq   = True
tab_ko       = []
tab_gap      = []


while not find and size_k > 1 :

    for i in reversed(range(len(sequence_1)-size_k+1)):
        chaine = ""
        for e in range(len(sequence_2)-size_k+1):

            seq_and_tsd = sequence_1[0:i] + "++:" + sequence_1[i:i+size_k] + ":++" + sequence_1[i+size_k:].strip() + "--|" + TE.strip() + "|--" + sequence_2[0:e] +  "++:" + sequence_2[e:e+size_k] + ":++" + sequence_2[e+size_k:].strip() 
            head = ""
            if strands :
                head = ":".join([sequence_1.strip(), sequence_2.strip(), str(i), str(e), str(size_k), str(size_flank), "(-)"])
            else:
                head = ":".join([sequence_1.strip(), sequence_2.strip(), str(i), str(e), str(size_k), str(size_flank), "(+)"])
            
            #print("---", sequence_1[i:i+size_k], sequence_2[e:e+size_k])
            alignments = pairwise2.align.globalxx(sequence_1[i:i+size_k], sequence_2[e:e+size_k])
            aln, ind = align(alignments, size_k)
            if aln and aln[2] == size_k :
                if i == len(sequence_1.strip()) - size_k and e == 0 and ((size_tsd != -1 and size_k == size_tsd) or size_tsd == -1) :
                    #print(sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[OK]", i, e, size_k)
                    chaine += "("+", ".join([sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[OK:"+id_elem+"]", str(i), str(e), str(size_k)])+")\n"
                    find = True
                else :
                    #print(sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[KO]", i, e, size_k)
                    chaine += "("+", ".join([sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[KO:"+id_elem+"]", str(i), str(e), str(size_k)])+")\n"

                if report_seq :
                    #print(">"+head+"\n"+seq_and_tsd+"\n")
                    chaine += ">"+head+"\n"+seq_and_tsd+"\n"

                tab_ko.append(chaine)
                if find:
                    print(chaine)
                    tmp_file = open("all_flank.fasta", "a")
                    tmp_file.write(sequence_1[i:i+size_k]+"\n"+sequence_2[e:e+size_k]+"\n")
                    tmp_file.close()
                    break
            elif aln and aln[2] == size_k - 1 and size_k >= 4 :
                if i == len(sequence_1.strip()) - size_k and e == 0 and ((size_tsd != -1 and size_k == size_tsd) or size_tsd == -1) :
                    #print(sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[OK]", i, e, size_k)
                    chaine += "("+", ".join([sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[K-O:"+id_elem+"]", str(i), str(e), str(size_k)])+")\n"
                    
                    if report_seq :
                        #print(">"+head+"\n"+seq_and_tsd+"\n")
                        chaine += ">"+head+"\n"+seq_and_tsd+"\n"
                    chaine += format_alignment(*alignments[ind])
                    tab_gap.append(chaine)
            # else :
            #   chaine += "("+", ".join([sequence_1[i:i+size_k], sequence_2[e:e+size_k], "[KO:"+id_elem+"]", str(i), str(e), str(size_k)])+")\n"
            #   seq_and_tsd = sequence_1[0:i] + "++:" + sequence_1[i:i+size_k] + ":++" + sequence_1[i+size_k:].strip() + "--|" + TE.strip() + "|--" + sequence_2[0:e] +  "++:" + sequence_2[e:e+size_k] + ":++" + sequence_2[e+size_k:].strip() 
            #   if report_seq :
            #           #print(">"+head+"\n"+seq_and_tsd+"\n")
            #           chaine += ">"+head+"\n"+seq_and_tsd+"\n"
            #   tab_ko.append(chaine)
            if find :
                break
        if find :
            break

    size_k -= 1


if not find :

    if len(tab_gap) :
        print(tab_gap[0])

    elif len(tab_ko) :
        print(tab_ko[0])


file.close()
file_sequence_TE.close()