import os
import sys
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import MDAnalysis as mda
import urllib.request
import seaborn as sns


from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral

import copy
sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, read_mapping_file, read_mapping_filePAIR, make_positive_angles, DihedralFromAtoms, calcDihedrals


dihedrals=[['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O1_M'],
					 ['M_G3_M', 'M_G3O1_M', 'M_G3P2_M', 'M_G3P2O2_M'],
					 ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
					 ['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
					 ['M_G3O1_M', 'M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M'],
					 ['M_G3P2_M', 'M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M'],
					 ['M_G3O3_M', 'M_G3C4_M', 'M_G3C5_M', 'M_G3N6_M'],
					 ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C6_M'],
					 ['M_G3O3_M','M_G3C4_M','M_G3C5_M','M_G3C5O1_M'],
					 ['M_G1C17_M', 'M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M'],
					 ['M_G1C16_M', 'M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M'],
					 ['M_G1C15_M', 'M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M'],
					 ['M_G1C14_M', 'M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M'],
					 ['M_G1C13_M', 'M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M'],
					 ['M_G1C12_M', 'M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M'],
					 ['M_G1C11_M', 'M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M'],
					 ['M_G1C10_M', 'M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M'],
					 ['M_G1C9_M', 'M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M'],
					 ['M_G1C8_M', 'M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M'],
					 ['M_G1C7_M', 'M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M'],
					 ['M_G1C6_M', 'M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M'],
					 ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2O1_M'],
					 ['M_G1C5_M', 'M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M'],
					 ['M_G1C4_M', 'M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M'],
					 ['M_G1C3_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],
					 ['M_G2C19_M', 'M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M'],
					 ['M_G2C18_M', 'M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M'],
					 ['M_G2C17_M', 'M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M'],
					 ['M_G2C16_M', 'M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M'],
					 ['M_G2C15_M', 'M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M'],
					 ['M_G2C14_M', 'M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M'],
					 ['M_G2C13_M', 'M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M'],
					 ['M_G2C12_M', 'M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M'],
					 ['M_G2C11_M', 'M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M'],
					 ['M_G2C10_M', 'M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M'],
					 ['M_G2C9_M', 'M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M'],
					 ['M_G2C8_M', 'M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M'],
					 ['M_G2C7_M', 'M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M'],
					 ['M_G2C6_M', 'M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M'],
					 ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2O1_M'],
					 ['M_G2C5_M', 'M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M'],
					 ['M_G2C4_M', 'M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M'],
					 ['M_G2C3_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M'],
					 ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G3_M'],
					 ['M_G1O1_M', 'M_G1_M', 'M_G2_M', 'M_G2O1_M'],
					 ['M_G2C2_M', 'M_G2O1_M', 'M_G2_M', 'M_G1_M'],
					 ['M_G1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
					 ['M_G2O1_M', 'M_G2_M', 'M_G3_M', 'M_G3O1_M'],
					 ['M_G2_M', 'M_G3_M', 'M_G3O1_M', 'M_G3P3_M'],
					 ['M_G1C2O1_M', 'M_G1C2_M', 'M_G1O1_M', 'M_G1_M'],           
					 ['M_G2C2O1_M', 'M_G2C2_M', 'M_G2O1_M', 'M_G2_M']
]

lipids = {'POPC','POPS','POPE','POPG'}
#lipids = {'POPG'}
#lipids = {'POPG','POPE','POPC'}
for i in dihedrals:
	calcDihedrals(lipids,i)

#mapping_file = "./mapping_files/mappingPOPEcharmm.txt"
#dihedrals = parseDihedralInput(mapping_file)
#parseDihedralInput(mapping_file)






