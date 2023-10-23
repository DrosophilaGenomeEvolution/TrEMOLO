# replace.py
# USAGE: python replace.py bad-word good-word target-file.txt
#
import sys
import re

search_term  = sys.argv[1]
replace_term = sys.argv[2]
target_file  = sys.argv[3]

with open(target_file, 'r') as file:
        content = file.read()

###content = content.replace(sys.argv[1], sys.argv[2])

content = re.sub(r"\b" + sys.argv[1] + r"\b", str(sys.argv[2]).replace("\\", "_"), content)

with open(target_file, 'w') as file:
        file.write(content)
        