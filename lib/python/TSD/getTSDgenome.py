import os
import sys
import argparse
from multiprocessing import Pool


#our
from utils import *

parser = argparse.ArgumentParser()

parser.add_argument("bed_seq", metavar='<bed-seq>', type=argparse.FileType('r'),
                    help="bed with sequence")

parser.add_argument("output_file", metavar='<output>', type=argparse.FileType('w'),
                    help="Name of tabular output file ")

#OPTION
parser.add_argument("-t", "--threads", dest='threads', type=int, default=4,
                    help="Number of threads")

args  = parser.parse_args()


def process_line(line):
    if len(line.strip().split("\t")) == 7 :
        infos, TE, ID, FK_L, FK_R, FKGL, FKGR = line.strip().split("\t")

        decalageR = 0
        decalageL = 0
        FK_M      = FK_R
        FK_S      = FK_L
        mode      = "R"
        if len(FK_R) > len(FK_L):
            mode = "L"
            FK_M = FK_L
            FK_S = FK_R

        sizeK = len(FK_M)
        find  = False
        while not find and sizeK > 3 :
            for i in range(len(FK_M)-sizeK+1):
                TSD = FK_M[i:i+sizeK]
                positionsMotif = BM(TSD, FK_S)
                if len(positionsMotif) > 0:
                    find = True
                    if mode=="R":
                        decalageL = len(FK_L)-sizeK-max(positionsMotif)
                        decalageR = i
                    else:
                        decalageL = len(FK_L)-sizeK-i
                        decalageR = min(positionsMotif)

                    position_on_genome = int(infos.split(":")[2].split("-")[0])
                    sizeTE = int(infos.split(":")[0])
                    sizeTE += int(decalageL) + int(decalageR)
                    print("NEXT:", f'{TE}|{ID}', positionsMotif, i, TSD, "OK", FK_L, FK_R, FK_M, FK_S, mode, decalageL, decalageR, infos)
                    return [f'{TE}|{ID}', TSD, str(position_on_genome-decalageL), str(sizeTE), infos]
                    break

            sizeK-=1
    else:
        print("Error number of column : ", str(line.strip()), file=sys.stderr)


def main():
    with Pool(processes=int(args.threads))  as pool:
        results = pool.map(process_line, args.bed_seq)

    for result in results:
        if result:
            args.output_file.write("\t".join(result)+"\n")


if __name__ == "__main__":
    main()
