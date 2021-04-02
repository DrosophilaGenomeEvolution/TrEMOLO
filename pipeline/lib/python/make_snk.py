import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser(description="build snakemake file")

#MAIN ARGS
parser.add_argument("name_instruction", type=str,
                    help="file instruction_file snakemake file")

parser.add_argument("name_file_rules", type=str,
                    help="file contains a list of rules (of snakemake file)")

parser.add_argument("name_out", type=str,
                    help="name of output snakemake file ")

parser.add_argument("-t", "--name_template", type=str, help="template file contain additionnel instruction file python for snakemake file")

args = parser.parse_args()

#ARGS INIT
name_instruction = args.name_instruction
name_file_rules  = args.name_file_rules
name_template    = args.name_template
name_out         = args.name_out

instruction_file = open(args.name_instruction, "r")
instructions     = instruction_file.readline().replace(" ", "").strip().split(">")

list_name_rule   = [i.split(":")[0] for i in instructions]

print("[" + sys.argv[0] + "] LIST RULES : ", " -> ".join(instructions))
os.system("rm -f " + name_out)


if name_template:
    os.system("cat " + name_template + " > " + name_out)


def cmd(cmd):
    proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    return str(out.decode("utf-8")).split("\\n")


nb_line_file_rules = cmd("wc -l " + name_file_rules + " | cut -d \" \" -f 1")[0]


for i, instruct in enumerate(instructions[::-1]):
    #os.system("rm -f tmpo"+str(i)+".snk")
    name_rule = instruct.split(":")[0]

    os.system("grep \"rule " + str(name_rule) + "\" " + name_file_rules + " -A " + str(int(nb_line_file_rules)) + " | grep \"^#END " + name_rule + "$\" -B " + str(int(nb_line_file_rules)) + " | grep -v \"^#END\" > tmpo"+str(i)+".snk")
    
    if len(instruct.split(":")) == 1  :
        if i == len(instructions) - 1 :
            os.system("sed -i 's/input_link=\[\],//g' tmpo" + str(i) + ".snk")
            os.system("sed -i 's/output_link=\[\],/temp(touch(\"" + name_rule + "\")),/g' tmpo" + str(i) + ".snk")
        else  :
            os.system("sed -i 's/output_link=\[\],/temp(touch(\"" + name_rule + "\")),/g' tmpo" + str(i) + ".snk")
            os.system("sed -i 's/input_link=\[\],/cout=[\"" + list_name_rule[::-1][i+1] + "\"],/g' tmpo"+ str(i) + ".snk")
    #suppression des input et output improviser
    elif len(instruct.split(":")) > 1 and instruct.split(":")[1] == "N":
        os.system("sed -i 's/input_link=\[.*\],//g' tmpo"+str(i)+".snk")
        os.system("sed -i 's/output_link=\[.*\],//g' tmpo" + str(i) + ".snk")
    elif len(instruct.split(":")) > 1 and instruct.split(":")[1] == "NI":
        os.system("sed -i 's/input_link=\[.*\],//g' tmpo"+str(i)+".snk")
        os.system("sed -i 's/output_link=\[\],/temp(touch(\"" + name_rule + "\")),/g' tmpo" + str(i) + ".snk")
    elif len(instruct.split(":")) > 1 and instruct.split(":")[1] == "NO":
        os.system("sed -i 's/output_link=\[.*\],//g' tmpo" + str(i) + ".snk")
        os.system("sed -i 's/input_link=\[\],/cout=[\"" + list_name_rule[::-1][i+1] + "\"],/g' tmpo"+ str(i) + ".snk")

    #print("grep \"rule " + str(name_rule) +"\" " + name_file_rules + " -A " + str(nb_line_file_rules) + " | grep \"^#END " + name_rule + "$\" -B 100000000 | grep -v \"^#END\" >> " + name_out)
    #os.system("grep \"rule "+ str(name_rule) +"\" list_rules.txt -A " + str(nb_line_file_rules) + " | grep \"^#END " + name_rule + "$\" -B 10000000 | grep -v \"^#END\" >> tmp.snk")
    os.system("cat tmpo"+str(i)+".snk >> " + name_out)
    os.system("rm -f tmpo*.snk")





