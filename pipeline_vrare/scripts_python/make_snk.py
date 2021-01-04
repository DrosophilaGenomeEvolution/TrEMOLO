import os
import sys
import argparse

parser = argparse.ArgumentParser(description="filters a blast file in output format 6 to keep the candidate candidate TE active")

#MAIN ARGS
parser.add_argument("name_instruction", type=str,
                    help="file instruction snake file")

parser.add_argument("name_file_rules", type=str,
                    help="file contains a list of rules (of snakefile)")



parser.add_argument("name_out", type=str,
                    help="name of tabular output file ")

parser.add_argument("-t", "--name_template", type=str, help="template file conatain additionnel instruction python for snakefile")

args = parser.parse_args()


name_instruction = args.name_instruction
name_file_rules  = args.name_file_rules
name_template    = args.name_template
name_out         = args.name_out

instruction      = open(args.name_instruction, "r")

list_instruction = instruction.readline().replace(" ", "").strip().split(">")

print("list instruction: ", list_instruction)

list_name_rule = [i.split(":")[0] for i in list_instruction]

os.system("rm -f " + name_out)

if name_template:
    os.system("cat " + name_template + " > " + name_out)

for i, instruct in enumerate(list_instruction[::-1]):
    #os.system("rm -f tmpo"+str(i)+".snk")
    name_rule = instruct.split(":")[0]

    #print(i, len(list_instruction)-i-1,  name_rule, i == len(list_instruction) - 1, len(list_instruction) )

    os.system("grep \"rule "+ str(name_rule) +"\" " + name_file_rules + " -A 100000000 | grep \"^#END " + name_rule + "$\" -B 100000000 | grep -v \"^#END\" > tmpo"+str(i)+".snk")
    if i == len(list_instruction) - 1 :
        os.system("sed -i 's/input_link=\[\],//g' tmpo"+str(i)+".snk")
        #os.system("sed -i 's/output_link=\[\],//g' tmpo"+str(i)+".snk")
        os.system("sed -i 's/output_link=\[\],/temp(touch(\"" + name_rule + "\")),/g' tmpo" + str(i) + ".snk")
    else :
        if i != len(list_instruction) - 1:
            os.system("sed -i 's/input_link=\[\],/cout=[\"" + list_name_rule[::-1][i+1] + "\"],/g' tmpo"+str(i)+".snk")


    #suppression des input et output improviser
    if len(instruct.split(":")) > 1 and instruct.split(":")[1] == "N": 
        os.system("sed -i 's/input_link=\[.*\],//g' tmpo"+str(i)+".snk")
        os.system("sed -i 's/output_link=\[.*\],//g' tmpo" + str(i) + ".snk")

    #print("grep \"rule " + str(name_rule) +"\" " + name_file_rules + " -A 100000000 | grep \"^#END " + name_rule + "$\" -B 100000000 | grep -v \"^#END\" >> " + name_out)
    #os.system("grep \"rule "+ str(name_rule) +"\" list_rules.txt -A 100000000 | grep \"^#END " + name_rule + "$\" -B 10000000 | grep -v \"^#END\" >> tmp.snk")
    os.system("cat tmpo"+str(i)+".snk >> " + name_out)





