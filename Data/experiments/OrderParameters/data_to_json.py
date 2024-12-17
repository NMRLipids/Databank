#!/usr/bin/env python3
""""
Converting OP dat file into properly organized JSON format.
- Amount of each molecule in the membrane is given as a ratio.
- For a membrane of single molecule it's 1.
"""

import json
import re
import math
import sys

from DatabankLib.jsonEncoders import CompactJSONEncoder

data_file = sys.argv[1]
print("Converting from: ", data_file)
outfile = data_file.replace(".dat", ".json")
print("Converting to: ", outfile)

data = {}

DEFAULT_OP_ERROR = 0.02

with open(data_file) as OPfile:
    lines = OPfile.readlines()
    for line in lines:
        line = re.sub(r'#.*$', '', line).strip()
        if line == "":
            continue
        lSplit = line.split()
        print(lSplit)
        assert len(lSplit) in [3, 4]
        if math.isnan(float(lSplit[2])):
            continue
        OPname = lSplit[0] + " " + lSplit[1]
        err = DEFAULT_OP_ERROR if len(lSplit) == 3 else float(lSplit[3])
        OPvalues = [float(lSplit[2]),  err]
        data[str(OPname)] = [OPvalues]
        # TODO: remove this double list thing!

with open(outfile, 'w') as f:
    json.dump(data, f, cls=CompactJSONEncoder)
