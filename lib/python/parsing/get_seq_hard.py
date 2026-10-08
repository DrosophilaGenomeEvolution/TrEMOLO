import argparse


parser = argparse.ArgumentParser(description="Get HARD sequence ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#MAIN ARGS
parser.add_argument("fil_index", metavar='<idx-file>', type=argparse.FileType('r'),
                    help="Input index TrEMOLO fasta")

parser.add_argument("name_fil_fasta", metavar='<name-file-fasta>', type=str,
                    help="Multi fasta file HARD")

parser.add_argument("fil_hard", metavar='<file-hard>', type=argparse.FileType('r'),
                    help="HARD TrEMOLO file")

#OPTION
parser.add_argument("-s", "--max-size", dest='max_size', type=int, default=None,
                    help="maximum sequence size")
args = parser.parse_args()


fil_index      = args.fil_index
name_fil_fasta = args.name_fil_fasta
fil_hard       = args.fil_hard
max_size       = args.max_size

#index to dic
dico_index_fasta = {}
line = fil_index.readline()
while line :
    dico_index_fasta[line.split(":")[1].strip()] = line.split(":")[0].strip()
    line = fil_index.readline()

fil_index.close()

# Build byte offsets, retaining only read names rather than the entire FASTA.
sequence_offsets = {}
with open(name_fil_fasta, "rb") as source:
    while True:
        header = source.readline()
        if not header:
            break
        if header.startswith(b">"):
            sequence_offsets[header[1:].strip().decode()] = source.tell()
fasta_source = open(name_fil_fasta, "rb")

#GET SEQ
line = fil_hard.readline()
while line :
    if line[0] != "#":
        spl = line.split("\t")
        chrom     = spl[0]
        start     = spl[1]
        end       = spl[2]
        id_reads  = spl[4]
        size_hard = int(spl[5])
        side      = spl[7]
        #print(id_reads, str(size_hard), str(side))
        num_line = int(dico_index_fasta[id_reads])
        #print("line:", num_line)
        if id_reads in sequence_offsets:
            fasta_source.seek(sequence_offsets[id_reads])
            seq = fasta_source.readline().strip().decode()
            if not max_size :
                if side == "L" :
                    #print("side:", str(side), "seq:", seq[:size_hard])
                    print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[:size_hard], str(side), str(size_hard)]))
                else :
                    #print("side:", str(side), "seq:", seq[-size_hard:])
                    print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[-size_hard:], str(side), str(size_hard)]))
            else :
                if size_hard > max_size :
                    if side == "L" :
                        #print("side:", str(side), "seq:", seq[:size_hard])
                        print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[:size_hard][-max_size:], str(side), str(max_size)]))
                    else :
                        #print("side:", str(side), "seq:", seq[-size_hard:])
                        print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[-size_hard:][:max_size], str(side), str(max_size)]))
                else :
                    if side == "L" :
                        #print("side:", str(side), "seq:", seq[:size_hard])
                        print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[:size_hard], str(side), str(size_hard)]))
                    else :
                        #print("side:", str(side), "seq:", seq[-size_hard:])
                        print("\t".join([str(chrom), str(start), str(end), str(id_reads), seq[-size_hard:], str(side), str(size_hard)]))

    line = fil_hard.readline()


fil_hard.close()



fasta_source.close()
