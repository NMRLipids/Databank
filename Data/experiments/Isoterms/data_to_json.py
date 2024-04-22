import sys
import json
import pandas as pd

data_file = sys.argv[1]
print(data_file)

# DAT file
if '.dat' in data_file:
    outfile = data_file.replace(".dat","_Isoterms.json")
    
    data = []

    with open(data_file) as PAfile:
        lines=PAfile.readlines()
        for line in lines:
            if "#" in line or line == "\n":
                continue
            values = [ float(i) for i in line.split() ] 

            data.append(values)
            
            with open(outfile, 'w') as f:
                json.dump(data,f)
    
# CSV file
if '.csv' in data_file:
    outfile = data_file.replace(".csv",".json")
    
    data = pd.read_csv( data_file, header=None )
    with open(outfile, 'w') as f:
        json.dump(  data.to_numpy().tolist(),f )
    
print(outfile)
