import sys
import json

data_file = sys.argv[1]
print(data_file)

if data_file.endswith(".dat"):
    ext = ".dat"
    sep = None
elif data_file.endswith(".csv"):
    ext = ".csv"
    sep = ","
    
outfile = data_file.replace(ext,"_NR.json")

data = []

with open(data_file) as PAfile:
    lines=PAfile.readlines()
    for line in lines:
        if "#" in line or line == "\n":
            continue
        values = [ float(i) for i in line.split(sep) ] 

        data.append(values)
        
        with open(outfile, 'w') as f:
            json.dump(data,f)
    
print(outfile)
