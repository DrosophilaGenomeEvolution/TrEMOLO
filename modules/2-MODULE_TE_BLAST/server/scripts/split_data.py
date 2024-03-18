import os
import sys
import json
import numpy as np

name_input  = sys.argv[1]
div         = int(sys.argv[2])
name_output = sys.argv[3]

json_file = open(name_input, "r")

json_object = json.loads(json_file.readline())


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)



size = len(json_object["data"])

for i in range(div):
	output_file = open(str(name_output) + "_split" + str(i) + ".json", "w")

	print(len(json_object["data"][i*int(size/div):(i+1)*int(size/div)]))
	
	##to json
	json_dico = json.dumps({"data": json_object["data"][i*int(size/div):(i+1)*int(size/div)]}, cls=NpEncoder)

	output_file.write(json_dico)
	output_file.close()



