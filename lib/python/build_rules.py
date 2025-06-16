import os
import sys
import argparse
import subprocess
from random import randrange

parser = argparse.ArgumentParser(description="build snakemake file")

#MAIN ARGS
parser.add_argument("name_instruction", type=str,
                    help="file instruction_file snakemake file")

parser.add_argument("name_file_rules", type=str,
                    help="file contains a list of rules (of snakemake file)")

parser.add_argument("name_out", type=str,
                    help="name of output snakemake file ")

parser.add_argument("-t", "--name_template", type=str, help="template file contain additionnel instruction file python for snakemake file")
parser.add_argument("-n", "--name_ID", type=str, help="put ID depending to the name")

args = parser.parse_args()

#ARGS INIT
name_instruction = args.name_instruction
name_file_rules  = args.name_file_rules
name_template    = args.name_template
name_out         = args.name_out
work             = ""
if args.name_ID == None:
    ID = randrange(200)
else:
    ID = os.path.basename(args.name_ID.rstrip("/"))
    work = args.name_ID.rstrip("/") + "/"

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


import re
import os

final_content = ""

for i, instruct in enumerate(instructions[::-1]):
    name_rule = instruct.split(":")[0]

    try:
        # Load rules from the source file
        with open(name_file_rules) as file:
            content = file.readlines()

        # Extract the corresponding rule
        start_pattern = rf"rule {re.escape(name_rule)}"
        end_pattern = rf"^#END {re.escape(name_rule)}$"
        inside_rule = False
        rule_content = []

        # Get rule
        for line in content:
            if re.match(start_pattern, line):
                inside_rule = True
            if inside_rule:
                rule_content.append(line)
            if re.match(end_pattern, line):
                inside_rule = False
                break

        if not rule_content:
            print(f"No rule found for {name_rule}")
            continue

        # Process the extracted content
        updated_content = []
        for line in rule_content:
            # Handle instructions with no additional flag
            if len(instruct.split(":")) == 1:
                if i == len(instructions) - 1:
                    # Remove lines containing `input_link=[]`
                    if re.search(r"input_link=\[\],", line):
                        continue
                    # Replace lines containing `output_link=[]`
                    line = re.sub(r'output_link=\[\],', f'temp(touch("{work}tmp_TrEMOLO_output_rule/rule_tmp_{name_rule}_{ID}")),', line)
                else:
                    # Replace lines containing `output_link=[]`
                    line = re.sub(r'output_link=\[\],', f'temp(touch("{work}tmp_TrEMOLO_output_rule/rule_tmp_{name_rule}_{ID}")),', line)
                    # Replace lines containing `input_link=[]`
                    line = re.sub(r'input_link=\[\],', f'cout=["{work}tmp_TrEMOLO_output_rule/rule_tmp_{list_name_rule[::-1][i+1]}_{ID}"],', line)

            elif len(instruct.split(":")) > 1:
                flag = instruct.split(":")[1]
                if flag == "N":
                    # Remove lines containing `input_link=[]` or `output_link=[]`
                    if re.search(r"input_link=\[.*\],", line) or re.search(r"output_link=\[.*\],", line):
                        continue
                elif flag == "NI":
                    # Remove lines containing `input_link=[]`
                    if re.search(r"input_link=\[.*\],", line):
                        continue
                    # Replace lines containing `output_link=[]`
                    line = re.sub(r'output_link=\[\],', f'temp(touch("{work}tmp_TrEMOLO_output_rule/rule_tmp_{name_rule}_{ID}")),', line)
                elif flag == "NO":
                    # Remove lines containing `output_link=[]`
                    if re.search(r"output_link=\[.*\],", line):
                        continue
                    # Replace lines containing `input_link=[]`
                    line = re.sub(r'input_link=\[\],', f'cout=["{work}tmp_TrEMOLO_output_rule/rule_tmp_{list_name_rule[::-1][i+1]}_{ID}"],', line)

            # Add the processed line
            updated_content.append(line)

        # Update the `step` field in the lines
        updated_content = [re.sub(r'step=0,', f'step={len(instructions) - i},', line) for line in updated_content]

        final_content += ''.join(updated_content) + "\n"

    except Exception as e:
        print(f"Error processing instruction {instruct}: {e}")

with open(name_out, "a") as output_file:
    output_file.write(final_content)
