import os
import numpy as np

import yaml
import json
import matplotlib.pyplot as plt
import MDAnalysis as mda
import urllib.request
import seaborn as sns

#from OrderParameter import *

from MDAnalysis.analysis.base import AnalysisBase
from MDAnalysis.lib.distances import calc_dihedrals
from MDAnalysis.analysis.dihedrals import Dihedral

import copy

class DihedralFromAtoms(AnalysisBase):
		"""Calculate dihedral angles for specified atomgroups.

		Dihedral angles will be calculated for each atomgroup that is given for
		each step in the trajectory. Each :class:`~MDAnalysis.core.groups.AtomGroup`
		must contain 4 atoms.

		Note
		----
		This class takes a list as an input and is most useful for a large
		selection of atomgroups. If there is only one atomgroup of interest, then
		it must be given as a list of one atomgroup.

		"""

		def __init__(self, atomgroups, orders, **kwargs):
				"""Parameters
				----------
				atomgroups : list
						a list of atomgroups for which the dihedral angles are calculated

				Raises
				------
				ValueError
						If any atomgroups do not contain 4 atoms

				"""
				super(DihedralFromAtoms, self).__init__(
						atomgroups[0].universe.trajectory, **kwargs)
				self.atomgroups = atomgroups

				if any([len(ag) != 4 for ag in atomgroups]):
						raise ValueError("All AtomGroups must contain 4 atoms")

				if len(atomgroups) != len(orders):
					raise ValueError("Order data should be provided for every atom group")

				self.ag1 = mda.AtomGroup([atomgroups[i][orders[i][0]] for i in range(len(atomgroups))])
				self.ag2 = mda.AtomGroup([atomgroups[i][orders[i][1]] for i in range(len(atomgroups))])
				self.ag3 = mda.AtomGroup([atomgroups[i][orders[i][2]] for i in range(len(atomgroups))])
				self.ag4 = mda.AtomGroup([atomgroups[i][orders[i][3]] for i in range(len(atomgroups))])


		def _prepare(self):
				self.angles = []

		def _single_frame(self):
				angle = calc_dihedrals(self.ag1.positions, self.ag2.positions,
															 self.ag3.positions, self.ag4.positions,
															 box=self.ag1.dimensions)
				self.angles.append(angle)

		def _conclude(self):
				self.angles = np.rad2deg(np.array(self.angles))

# Download link
def download_link(doi, file):
	if "zenodo" in doi.lower():
		zenodo_entry_number = doi.split(".")[2]
		return 'https://zenodo.org/record/' + zenodo_entry_number + '/files/' + file
	else:
		print ("DOI provided: {0}".format(doi))
		print ("Repository not validated. Please upload the data for example to zenodo.org")
		return ""
		
# read mapping file
def read_mapping_file(mapping_file, atom1):
	with open(mapping_file, 'rt') as mapping_file:
		for line in mapping_file:
			if atom1 in line:
				m_atom1 = line.split()[1]
	return m_atom1

def read_mapping_filePAIR(mapping_file, atom1, atom2):
	with open(mapping_file, 'rt') as mapping_file:
		print(mapping_file)
		for line in mapping_file:
			if atom1 in line:
				m_atom1 = line.split()[1]
#                    print(m_atom1)
			if atom2 in line: 
				m_atom2 = line.split()[1]
#                    print(m_atom2)
	return m_atom1, m_atom2

def make_positive_angles(x):
	for i in range(len(x)):
		if x[i] < 0:
			x[i] = x[i] + 360
		else:
			x[i] = x[i]
	return x


def calcDihedrals(lipids,DIHatoms):
	colors = {'POPC' :'black','POPS':'red','POPE':'blue','POPG':'green'}

	found = 0
	
	for subdir, dirs, files in os.walk(r'../../Data/Simulations/'):
		for filename in files:
			if filename.endswith("README.yaml"):
				READMEfilepath = subdir + '/README.yaml'
				found +=1
				print(READMEfilepath)
				with open(READMEfilepath) as yaml_file:
					readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
					for molname in lipids:
						doi = readme.get('DOI')
						trj = readme.get('TRJ')
						tpr = readme.get('TPR')
						trj_name = subdir + '/' + readme.get('TRJ')[0][0]
						tpr_name = subdir + '/' + readme.get('TPR')[0][0]
						gro_name = subdir + '/conf.gro'
						trj_url = download_link(doi, trj[0][0])
						tpr_url = download_link(doi, tpr[0][0])

						#Download tpr and xtc files to same directory where dictionary and data are located
						if (not os.path.isfile(tpr_name)):
							response = urllib.request.urlretrieve(tpr_url, tpr_name)
												
						if (not os.path.isfile(trj_name)):
							response = urllib.request.urlretrieve(trj_url, trj_name)
						
						if sum(readme['N' + molname]) > 0:
							print('Analyzing '+molname+' in '+READMEfilepath)
							#fig= plt.figure(figsize=(12,9))
							if (not os.path.isfile(gro_name)):
								os.system('echo System | gmx trjconv -f {} -s {}  -dump 0 -o {}'.format(trj_name,tpr_name,gro_name))
						
							xtcwhole=subdir + '/whole.xtc'
							if (not os.path.isfile(xtcwhole)):
								os.system('echo System | gmx trjconv -f {} -s {} -o {} -pbc mol '.format(trj_name,tpr_name,xtcwhole))
					
							try:
								traj = mda.Universe(tpr_name,xtcwhole)
							except FileNotFoundError or OSError:
								continue
												
								#             #try:
								#             #    traj = mdtraj.load(xtcwhole, top = gro_name)
								#                 #print(lipid)
								#             #except FileNotFoundError or OSError:
								#             #    continue
							
							mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT'][molname] # readme.get('MAPPING')[0][0]
							try:
								atom1 = read_mapping_file(mapping_file, DIHatoms[0])
								atom2 = read_mapping_file(mapping_file, DIHatoms[1])
								atom3 = read_mapping_file(mapping_file, DIHatoms[2])
								atom4 = read_mapping_file(mapping_file, DIHatoms[3])
								print(atom1,atom2,atom3,atom4)
							except:
							#print(atom1 + " and " + atom2 + " not found in the mapping file.")
								print("Some atom not found in the mapping file.")
								continue

							#print(traj.residues.n_residues)

							ags = []
							orders = []
							for residue in traj.select_atoms("name {} {} {} {}".format(atom1,atom2,atom3,atom4)).residues:
								#print(residue.resid)
								#print(residue.atoms.names)
								atoms=traj.select_atoms("name {} {} {} {} and resid {}".format(atom1,atom2,atom3,atom4,str(residue.resid)))
								if len(atoms) == 4:
									#print(atoms.names)
									ags.append([])
									ags[-1] = copy.deepcopy(atoms)
									orders.append([
									np.where(atoms.names==atom1)[0][0],
									np.where(atoms.names==atom2)[0][0],
									np.where(atoms.names==atom3)[0][0],
									np.where(atoms.names==atom4)[0][0]])
							#print(len(ags))
							R = DihedralFromAtoms(ags,orders).run()
							dihRESULT = R.angles.T
							
							dihRESULT = [make_positive_angles(x) for x in dihRESULT ]
							distSUM = np.zeros(360)
							for i in dihRESULT:
								distSUM += np.histogram(i,np.arange(0,361,1),density=True)[0]
							#dihRESULT = [make_positive_angles(x) for x in dihRESULT ]
							#dist = [ 0 for i in range(len(dihRESULT))]
							#distSUM = [ 0 for i in range(360)]
							#for i in range(len(dihRESULT)):
							#	dist[i] =  plt.hist(dihRESULT[i], range(360),density=True);
							#	distSUM = np.add(distSUM,dist[i][0])
												
										
							distSUM = [x / len(dihRESULT) for x in distSUM]
							xaxis = (np.histogram(i,np.arange(0,361,1),density=True)[1][1:]+
								np.histogram(i,np.arange(0,361,1),density=True)[1][:-1])/2.
							#xaxis = [ 0 for i in range(len(dist[0][1])-1)]
							#for i in range(len(dist[0][1])-1):
							#	xaxis[i]=(dist[0][1][i])
											
							#plt.plot(xaxis,distSUM,color=colors[molname], label = readme.get('SYSTEM'))[0] 
							#plt.legend()
							#plt.xlabel("Angle (°)")
							#plt.show()
									
							dihedralFOLDERS = subdir.replace("Simulations","dihedral")
							os.system('mkdir -p {}'.format(dihedralFOLDERS))
							os.system('cp {} {}'.format(READMEfilepath,dihedralFOLDERS))
							outfile=open(str(dihedralFOLDERS) + '/' + molname + "_" + DIHatoms[0] + "_" + DIHatoms[1] + "_" + DIHatoms[2] + "_" + DIHatoms[3] +'.dat','w')
							for i in range(len(xaxis)):
									outfile.write(str(xaxis[i]) + " " + str(distSUM[i])+'\n')
							outfile.close()
							#plt.close()

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






