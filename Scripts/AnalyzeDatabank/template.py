import os, sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from DatabankLib.databankLibrary import initialize_databank, lipids_dict

systems = initialize_databank()

for system in systems:
    print(system)
        

