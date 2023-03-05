# replace.py
# USAGE: 
#
import sys
import re
import argparse

parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("term_file", metavar='<term-file>', option_strings=['term-file'], type=argparse.FileType('r'),
                    help="")

parser.add_argument("target_file", metavar='<target-file>', type=argparse.FileType('r'),
                    help="")

parser.add_argument("column", metavar='<coulmn>', type=int,
                    help="")


args = parser.parse_args()


column_term = args.column

file_term = args.term_file
file_term.readline() ##pass header

##get pseudo = original
dico = {}
for line in file_term :
    dico[line.strip().split("\t")[1]] = line.strip().split("\t")[0]
      
file_term.close()

file_target = args.target_file
for line in file_target :
    pseudo   = line.strip().split("\t")[column_term - 1]
    if pseudo in dico:
        original = dico[pseudo]
        new_line = re.sub(r"\b" + str(pseudo) + r"\b", str(original).replace("\\", "_"), line)
        print(new_line.strip())
    else:
        print(line.strip())


file_target.close()