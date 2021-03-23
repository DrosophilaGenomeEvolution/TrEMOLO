
NAME_FILE_RULES=$1
INSTRUCTION_FILE=$2
NAME_OUTPUT=$3

nb_line_file=`wc -l ${NAME_FILE_RULES} | cut -d \" \" -f 1`
i=0
for rule in `cat INSTRUCTION_FILE | tr -d " " | tr ">" "\n" `;
do
    echo $i;
    nb_line_file
    grep "rule $rule" ${NAME_FILE_RULES} -A $nb_line_file | grep "^#END $rule$" -B $nb_line_file | grep -v "^#END" > tmpo${i}.snk
    if [[ `cat INSTRUCTION_FILE | wc -l` -eq $(($i-1)) ]]; then
        sed -i 's/input_link=\[\],//g' tmpo${i}.snk
        sed -i 's/output_link=\[\],/temp(touch(${rule})),/g' tmpo${i}.snk
    else
        sed -i 's/input_link=\[\],/cout=[\"" + list_name_rule[::-1][i+1] + "\"],/g' tmpo${i}.snk
    fi

    i=$(($i+1))
done;

cat tmpo*.snk >> ${NAME_OUTPUT}

os.system("grep \"rule " + str(name_rule) + "\" " + name_file_rules + " -A " + str(int(nb_line_file_rules)) + " | grep \"^#END " + name_rule + "$\" -B " + str(int(nb_line_file_rules)) + " | grep -v \"^#END\" > tmpo"+str(i)+".snk")
if i == len(instructions) - 1 :
    os.system("sed -i 's/input_link=\[\],//g' tmpo" + str(i) + ".snk")
    os.system("sed -i 's/output_link=\[\],/temp(touch(\"" + name_rule + "\")),/g' tmpo" + str(i) + ".snk")
elif i != len(instructions) - 1 :
    os.system("sed -i 's/input_link=\[\],/cout=[\"" + list_name_rule[::-1][i+1] + "\"],/g' tmpo"+str(i)+".snk")

#suppression des input et output improviser
if len(instruct.split(":")) > 1 and instruct.split(":")[1] == "N": 
    os.system("sed -i 's/input_link=\[.*\],//g' tmpo"+str(i)+".snk")
    os.system("sed -i 's/output_link=\[.*\],//g' tmpo" + str(i) + ".snk")

#print("grep \"rule " + str(name_rule) +"\" " + name_file_rules + " -A " + str(nb_line_file_rules) + " | grep \"^#END " + name_rule + "$\" -B 100000000 | grep -v \"^#END\" >> " + name_out)
#os.system("grep \"rule "+ str(name_rule) +"\" list_rules.txt -A " + str(nb_line_file_rules) + " | grep \"^#END " + name_rule + "$\" -B 10000000 | grep -v \"^#END\" >> tmp.snk")
os.system("cat tmpo"+str(i)+".snk >> " + name_out)
