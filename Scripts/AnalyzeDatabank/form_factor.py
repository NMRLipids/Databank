#!/usr/bin/env python
# coding: utf-8

# In[2]:

import re
import MDAnalysis as mda
import numpy as np
import sys
import matplotlib.pyplot as plt
nav = 6.02214129e+23
import json
from json import JSONEncoder
import time
import gc
gc.collect()
import os
import yaml

from databankLibrary import lipids_dict


# In[2]:

# to write data in numpy arrays into json file
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)
        
def loadMappingFile(path_to_mapping_file):
    # load mapping file into a dictionary
    mapping_dict = {}
    with open('../BuildDatabank/mapping_files/'+path_to_mapping_file, "r") as yaml_file:
        mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
    yaml_file.close()
    
    return mapping_dict
        
def getLipids(readme, molecules=lipids_dict.keys()):
    lipids = 'resname '

    for key1 in readme['COMPOSITION'].keys():
        #print('TST',key1)
        if key1 in molecules:
            mapping_file = readme['COMPOSITION'][key1]['MAPPING']
            
            mapping_dict = loadMappingFile(mapping_file)
            
            switch = 0
            for key2 in mapping_dict:
                #print(mapping_dict[key2]['RESIDUE'])
                #if 'POPC' in key1:
                #    res = mapping_dict[key2]['RESIDUE']
                #    print('TST',res)
                try:
                    res = mapping_dict[key2]['RESIDUE']
                    #print('TST',res)
                except KeyError:
                    switch = 1
                    continue
                else:
                    if res not in lipids:
                        lipids = lipids + res + ' or resname '
            
            if switch == 1:
                lipids = lipids + readme['COMPOSITION'][key1]['NAME'] + ' or resname '
                #break

                     
    lipids = lipids[:-12]
  #  print(lipids)
    return lipids

def getWater(readme, molecules=lipids_dict.keys()):
    waters = 'resname ' + readme['COMPOSITION']['SOL']['NAME']
    return waters
    
 


#for testing purposes to safe time for loading trajectory and creating dictonary
#the electron.dat will be decripted in the future as the values will be used from the universal mapping file


# modify ion electron numbers Na, Cl, Ca, Mg, K

electron_dictionary={"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
                   "Na":10,"Mg":10,"Al":13,"Si":14,"P":15,"S":16,"Cl":18,"Ar":18,
                   "K":18,"Ca":18,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,
                     "Cu":29,"Zn":30,"Ga":31,"Ge":32,"Vi":0,"Cs":54,"D":0}
	
def filterHbonds(mapping_names):
    
    re_H = re.compile(r'M_([A-Z]{1,2}[0-9]{1,4})*H[0-4]{1,2}_M')
    
    filtered = list(filter(re_H.match, mapping_names))
    
    return filtered

#def listNamePairs(mapping_dict):
#    pairs = []
#    residues = []
        
#    for mapping_name in mapping_dict.keys():
#        pairs.append([mapping_name, mapping_dict[mapping_name]['ATOMNAME']])
        
#        try:
#            residues.append(mapping_dict[mapping_name]['RESIDUE'])
#        except KeyError:
#            continue
    
#    if residues:
#        residues = list(dict.fromkeys(residues))
    
    
    
#    return pairs, residues #returning a dictionary with residues as keys and lists as values maybe better?


def listNamePairs(mapping_dict,molecule):
    pairs_dictionary = {}
    residues = []
        
        
    for mapping_name in mapping_dict.keys():
        
        try:
            residues.append(mapping_dict[mapping_name]['RESIDUE'])
        except KeyError:
            continue
            
   # print(residues)
    
    if residues:
        pairs_dictionary = dict.fromkeys(residues)
        print(pairs_dictionary.keys())
        for res in pairs_dictionary.keys():
            pairs = []
            for mapping_name in mapping_dict.keys():
                if res == mapping_dict[mapping_name]['RESIDUE']:
                    pairs.append([mapping_name, mapping_dict[mapping_name]['ATOMNAME']]) 
            pairs_dictionary[res] = pairs

    else:
        pairs_dictionary = dict.fromkeys([molecule])  
        pairs = []  
        for mapping_name in mapping_dict.keys():
            pairs.append([mapping_name, mapping_dict[mapping_name]['ATOMNAME']]) 
        pairs_dictionary[molecule] = pairs
    
    
    return pairs_dictionary


class FormFactor:
    """
    Calculates form factors from density profiles.
    
    Further development will include thickness of the membrane 
    from the intersection of lipid density and water density.
    Already enables to calculate electrom/mass/number densities
    
    Density could be calculated only from the final form factor profile 
    - testing is needed to see stability for rare species.
    
    Examination of FF error spikes is needed!
    
    """
    def __init__(self, path, conf,traj,nbin,output,readme,density_type="electron"):
        self.path = path
        self.conf = conf
        self.traj = traj
        self.readme = readme
        start_time=time.time()
        try:
            #print(self.conf,self.traj)
            self.u = mda.Universe(self.conf,self.traj)
            print()
        except:
            gro = self.path + '/conf.gro'
            print("Generating conf.gro because MDAnalysis cannot read tpr version")
            os.system('echo System | gmx trjconv -s '+ self.conf + ' -f '+ self.traj + ' -dump 0 -o ' + gro)
            self.conf = gro
            self.u = mda.Universe(self.conf,self.traj)
            
        print("Loading the trajectory takes {:10.6f} s".format(time.time()-start_time))
        
        
        
        
        #number of bins
        self.nbin = nbin
        #the totatl box size in [nm] - will be probably removed and tpr box size or 
        #transform of final FF will be used instead

        self.output = path + output
        self.lipids = getLipids(readme)
        self.waters = getWater(readme)
        
        self.density_type = density_type
        
        self.system_mapping = self.temporary_mapping_dictionary()
        
        self.calculate_weight()

        self.calculate_density()
        
        
        
        
    def temporary_mapping_dictionary(self):
        mapping_dictionary={}
        for key1 in self.readme['COMPOSITION'].keys(): #
        
            key2 = self.readme['COMPOSITION'][key1]['NAME']
            mapping_file = self.readme['COMPOSITION'][key1]['MAPPING']
            mapping_dict = loadMappingFile(mapping_file)
        
            # put mapping files into a larger dictionary where molecule name is key and the 
            # dictionary in the correponding mapping file is the value
            mapping_dictionary[key2]= mapping_dict #try if this works
        
        return mapping_dictionary

    def getElectrons(self,mapping_name):
        name1 = re.sub(r'[0-9]*_M$', '', mapping_name) # removes numbers and '_M' from the end of string
        name2 = re.sub(r'^M_([A-Z]{1,2}[0-9]{1,4})*','',name1) # removes M_X12Y23... sequence from the beginning
        # name2 is the atom and electrons are assigned to this atom

        if name2 == 'G': # G is a glycerol carbon so change G to C
            name2 = 'C'

        try:
            el = electron_dictionary[name2]
        except KeyError:
            print (f'ERROR: This mapping name cannot be read by our rules: {mapping_name}', file=sys.stderr)
            print ('Consider changing naming in your mapping file.')
            sys.exit(1)

        return el

    def residueElectronsAll(self,molecule):
        print(molecule)
        "return list of electrons"
        electrons = []
        #get the name of molecule used in simulation files
        molname = self.readme['COMPOSITION'][molecule]['NAME']
      
        pairs_residue = listNamePairs(self.system_mapping[molname], molname)
        
        #print(pairs_residue)
        
        # if lipid is split to multiple residues
        selection_txt = ""
        
        try:
            selection_txt = "moltype " + self.readme['COMPOSITION'][molecule]['MOLTYPE'] + " and "
        except KeyError:
            pass
            
        
        for res in pairs_residue.keys(): #KORJAA TOI res tossa explicit_atoms
          #  print(res)                           #   when lipids are split into more than one residue in topology with same residue names e.g. PGR,PA,OL and PE, PA,OL
                                                  #    u.select_atoms("resname " + res) does not distinguish which PA or PL belongs to POPG or POPE
            selection_txt = selection_txt + "resname " + res
            explicit_atoms = list(self.u.select_atoms(selection_txt).atoms.names) #explicit atoms of the residue
            for atom in explicit_atoms:
                #find generic mapping name matching to forcefield atom name
                index_list = [atom in pairs for pairs in pairs_residue[res]] #find mapping name
                #print(atom, index_list)
                #print(atom_i1)
                atom_i1 = index_list.index(True)
                #print("atom_i1")
  
                mapping_name = pairs_residue[res][atom_i1][0] 
                e_atom_i = self.getElectrons(mapping_name) #get number of electrons in an atom i of residue
                #print(mapping_name + " " + str(e_atom_i))
                electrons.append(e_atom_i)

    #    print(electrons)
    #    print(len(electrons))
        return electrons
    
        
    def assignElectrons(self):
        # check if simulation is united atom or all atom
        try:
            UA = self.readme['UNITEDATOM_DICT']
        except KeyError:
            UA = False
        else:
            UA = True
        
        weights=[]
        l_weights = [] #separate list for lipid electrons
        w_weights = [] #separate list for water electrons
        
        if UA:  #UNITED ATOM SIMULATIONS
            for molecule1 in self.readme['COMPOSITION'].keys():
                print(molecule1)
                
                if molecule1 in lipids_dict:
                    molecule2 = self.readme['COMPOSITION'][molecule1]['NAME']
                    electrons=[]

                    #print(molecule2, self.system_mapping[molecule2])
                    pairs_residue = listNamePairs(self.system_mapping[molecule2], molecule2)
                    atomsH = filterHbonds(self.system_mapping[molecule2].keys()) #list of hydrogen atoms
                    #print(pairs_residue)
                    
                    selection_txt = ""
        
                    try:
                        selection_txt = "moltype " + self.readme['COMPOSITION'][molecule2]['MOLTYPE'] + " and "
                    except KeyError:
                        pass
                    
                    # extract explicit atoms and get the mapping names
                    for res in pairs_residue.keys():
                        explicit_atoms = list(self.u.select_atoms("resname " + res).atoms.names) #explicit atoms of the residue
                        for atom in explicit_atoms:
                            #print(atom)
                            #find generic mapping name matching to forcefield atom name
                            index_list = [atom in pairs for pairs in pairs_residue[res]] 
                            atom_i1 = index_list.index(True)
                            mapping_name = pairs_residue[res][atom_i1][0]
                            name1 = re.sub(r'_M','',mapping_name[::1]) #remove _M from the end of mapping name
                
                            #count how many times name1 occurs in atomsH i.e. how many H are bound to it
                            numberH = sum(1 for atomH in atomsH if name1 == re.sub(r'H[1-4]_M','',atomH[::1]))
                            number_e = self.getElectrons(mapping_name)
                            e_atom_i = number_e + numberH 
                            #print(name1 + " " + str(numberH))
                            electrons.append(e_atom_i)
                            
                    
                    weights.extend(electrons)
                    l_weights.extend(electrons)
                else:
                    w = self.residueElectronsAll(molecule1)
                    weights.extend(w)    
                    if molecule1 == 'SOL':
                        w_weights.extend(w)
                    
        else:   #ALL ATOM SIMULATIONS
            for molecule in self.readme['COMPOSITION'].keys():
                w = self.residueElectronsAll(molecule)
                weights.extend(w)
                if molecule in lipids_dict:
                    l_weights.extend(w)
                elif molecule == 'SOL':
                    w_weights.extend(w)
                
                    
           # how to do lipids split into three residues??    MAYBE WORKS 
              
              
      #  print("number of SOL: " + str(cSOL))
      #  print("number of NA: " + str(cNA))
             
        return weights, l_weights, w_weights
        
              
    def calculate_weight(self):
        """
        Creates an array of weights for all atoms in the simulation.
        
        For electron densities:
         - creates dictonary of atom types and number of electrons
           loaded form an external file --> in the future it backmaps 
                                            the atom names to mapping files
                                            and automaticaly assigns # of electrons
                                            
        Number densities:
         - all weights are 1
         
        Mass densities:
         - reads masses from u.atoms.masses
        """
        start_time=time.time()
        if self.density_type=="electron":
            #weights to calculate the electron density
            self.wght, self.lipid_wght, self.water_wght = self.assignElectrons() # function that maps the electrons to atoms
            print("lenght of wght: " + str(len(self.wght)))
            print("lenght of lipid_wght: " + str(len(self.lipid_wght)))
            #print(self.wght)
            print("atoms in system: " + str(len(self.u.atoms.names)))
            
            if len(self.wght) != len(self.u.atoms.names):
                print("ERROR: Number of atoms mismatch")
                print("ERROR: If lipids are split to many residues add moltype of lipids from tpr as 'MOLTYPE' to README.yaml under 'COMPOSITION'.")
            
           # print(self.lipid_wght)
            #print(self.water_wght)
            
           
        if self.density_type=="number":
            self.wght=np.ones(self.u.atoms.names.shape[0])
            clipids = self.u.select_atoms(self.lipids)
            self.lipid_wght=np.ones(clipids.atoms.names.shape[0])
            cwaters = self.u.select_atoms(self.waters)
            self.lipid_wght=np.ones(cwaters.atoms.names.shape[0])
            
        if self.density_type=="mass":
            self.wght=self.u.atoms.masses
            clipids = self.u.select_atoms(self.lipids)
            self.lipid_wght=np.ones(clipids.atoms.names.shape[0])
            cwaters = self.u.select_atoms(self.waters)
            self.lipid_wght=np.ones(cwaters.atoms.names.shape[0])           
        
        #print(self.lipids)
        #print(self.u.atoms.resnames)
            
        print("Creating the electron mapping dictonary takes {:10.6f} s".format(time.time()-start_time))
        
        

    def calculate_density(self):

        
        c = self.u.select_atoms(self.lipids)

        print(c)
        
        box_z = self.u.dimensions[2] # + 10 if fails
        print(box_z)
        d =  box_z/ 10/self.nbin #     # bin width
        boxH = box_z/10
        print(boxH)
        x = np.linspace(-boxH/2,boxH/2,self.nbin+1)[:-1] + d/2
        density_z_centered = np.zeros(self.nbin)
        density_z_no_center = np.zeros(self.nbin)
        density_lipids_center = np.zeros(self.nbin)
        density_waters_center = np.zeros(self.nbin)
        
        #for running FF - not needed now
        #fa=[]
        #fb=[]

        """Calculte density profiles and FF from individual frames"""
        start_time=time.time()
        min_z=10000000
        frame=0

#        for mol in self.readme['COMPOSITION']:
 #           print(self.readme['COMPOSITION'][mol]['MAPPING'])

        ElectronNumbers = {}

        try:
            UA = self.readme['UNITEDATOM_DICT']
        except KeyError:
            UA = False
        else:
            if self.readme['UNITEDATOM_DICT'] == None:
                #print(self.readme['UNITEDATOM_DICT'])
                UA = False 
            else:
                UA = True
        
        #print(self.system_mapping)

            
        #for mol in self.system_mapping:
        for key1 in self.readme['COMPOSITION'].keys(): #
            mol = self.readme['COMPOSITION'][key1]['NAME']
            print(mol)
            #ElectronNumbers[mol] = {}
            for universalAN in self.system_mapping[mol]:
                AtomName = self.system_mapping[mol][universalAN]['ATOMNAME']
                try:
                    ResName = self.system_mapping[mol][universalAN]['RESIDUE']
                except:
                    ResName = mol

                #print(ElectronNumbers.keys())
                if ResName not in ElectronNumbers.keys():
                    ElectronNumbers[ResName] = {}
                
                #print(mol)

                if UA and key1 in lipids_dict:
                    #print(self.readme['UNITEDATOM_DICT'])
                    UAlipidjsonNAME = './lipid_json_buildH/' + self.readme['UNITEDATOM_DICT'][key1] + '.json'
                    with open(UAlipidjsonNAME) as json_file:
                        UAlipidjson = json.load(json_file)
                    #print(UAlipidjson)

                    numberH = 0
                    try:
                        #print(UAlipidjson[AtomName][0])
                        if UAlipidjson[AtomName][0] == "CH":
                            numberH = 1
                        if UAlipidjson[AtomName][0] == "CH2":
                            numberH = 2
                        if UAlipidjson[AtomName][0] == "CH3":
                            numberH = 3
                    except:
                        pass
                        #print(AtomName)

                    ElectronNumbers[ResName][AtomName] = self.getElectrons(universalAN) + numberH

                else:
                    #print(universalAN, self.system_mapping[mol][universalAN]['ATOMNAME'])
                    ElectronNumbers[ResName][AtomName] = self.getElectrons(universalAN)
                    #print(universalAN, self.getElectrons(universalAN))


                    
        print(ElectronNumbers)
                
        for ts in self.u.trajectory:
            #count the index of the frame, numbered from 0, used to be used for the density profile averaging
            #posible not needed now
            frame += 1 #ts.frame
            #print(frame)
      
            

            
            #reads the dimension in z-direction
            box_z = ts.dimensions[2]
            if box_z/10<min_z:
                min_z=box_z/10
            
            #print(min_z)

            
            #reads the coordinates of all of the atoms
            crds = self.u.atoms.positions

            
            crds_no_center = c.atoms.positions[:,2]/10
          #  lipid_wght=np.ones(c.atoms.names.shape[0])
            
            #calculates the center of mass of the selected atoms that the density should be centered around and 
            #takes the z-coordinate value
            ctom = c.atoms.center_of_mass()[2]
            
            #moves the center of mass of the selected centering group into box/2
            crds[:,2] += box_z/2 - ctom

            #clipids_center += box_z/2 - ctom
            #cwaters_center += box_z/2 - ctom            

            #clipids_center /= 10
            #cwaters_center /= 10            
           
            
            """shifts the coordinates in the universe by the value of the center of mass"""
            self.u.atoms.positions = crds
            
            """puts the atoms back to the original box dimension; it possibly does not take PBC into account
            #therefore it may brake some of the water molecules; try it, come to the issue later"""
            self.u.atoms.pack_into_box()

            clipids = self.u.select_atoms(self.lipids)
            cwaters = self.u.select_atoms(self.waters)

            #print(self.lipids)
            
            """shif the coordinates so that the center in z-dimention is in 0; 
            #divide by 10 to get the coordinates in nm, since now the crds are only the z coordinates"""
            crds = (self.u.atoms.positions[:,2] - box_z/2)/10
            clipids_center = (clipids.atoms.positions[:,2] - box_z/2)/10
            cwaters_center = (cwaters.atoms.positions[:,2] - box_z/2)/10

            if frame == 1:
                weightsALL = []
                for i in self.u.atoms:
                    #print(i.resname)
                    weightsALL.append(ElectronNumbers[i.resname][i.name])

                weightsLIPIDS = []
                for i in clipids.atoms:
                    weightsLIPIDS.append(ElectronNumbers[i.resname][i.name])

                weightsWATERS = []
                for i in cwaters.atoms:
                    weightsWATERS.append(ElectronNumbers[i.resname][i.name])

                
                #print(self.u.atoms)
            #print(weights)
            
         #   water_wght=np.ones(cwaters.atoms.names.shape[0])
            
            """calculates the volume of the bin; d- the "height" of a bin; assumes in [nm] """
            # ts.dimension[0], ts.dimension[1] - the x and y dimension; in [A] --> devides by 100
            vbin = d*np.prod(ts.dimensions[:2])/100 #needed for density, correct!
            
            #start_time2=time.time()

            
                
            
            
            
            """Do running ff - seems to produce consistent results""" #we will not use now           
            #box_z /= 10 #gets the box-z size in nm
            #d_ff =  box_z/ self.nbin
            #vbin_FF =d_ff*np.prod(ts.dimensions[:2])/100 # volume of the current bin in running FF
            #ff_density=np.histogram(crds,bins=self.nbin,range=(-box_z/2,box_z/2),weights=self.wght/vbin_FF)[0]
            #FF_range = np.linspace(0,999,1000)  
            #fa_run,fb_run=self.fourier(ff_density,ts.dimensions[2]/10,FF_range,d_ff)
            #fourrier_result2= np.sqrt(np.multiply(fa_run,fa_run)+np.multiply(fb_run,fb_run))
            #fourrier_data2 = np.vstack((FF_range*0.1*0.01,fourrier_result2)).transpose()
            #self.plot_fourier_final(fourrier_data2)

            
            #fa.append(fa_run)
            #fb.append(fb_run)
            #print("One furier loop takes {:10.6f} s".format(time.time()-start_time2))
          
 
            #print(len(clipids_center))

            """calculates the total density profile; keep for now"""


            #density_z_centered += np.histogram(crds,bins=self.nbin,range=(-boxH/2,boxH/2),weights=self.wght/vbin)[0]
            #density_z_no_center += np.histogram(crds_no_center,bins=self.nbin,range=(0,boxH),weights=self.lipid_wght/vbin)[0]
            #density_lipids_center += np.histogram(clipids_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=self.lipid_wght/vbin)[0] # ELECTRON WEIGHTS 
            #density_waters_center += np.histogram(cwaters_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=self.water_wght/vbin)[0] # ELECTRON WEIGHTS

            density_z_centered += np.histogram(crds,bins=self.nbin,range=(-boxH/2,boxH/2),weights=weightsALL/vbin)[0]
            density_z_no_center += np.histogram(crds_no_center,bins=self.nbin,range=(0,boxH),weights=weightsLIPIDS/vbin)[0]
            density_lipids_center += np.histogram(clipids_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=weightsLIPIDS/vbin)[0] # ELECTRON WEIGHTS 
            density_waters_center += np.histogram(cwaters_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=weightsWATERS/vbin)[0] # ELECTRON WEIGHTS

            
            #self.wght=np.ones(self.u.atoms.names.shape[0])
            #clipids = self.u.select_atoms(self.lipids)
            #self.lipid_wght=np.ones(clipids.atoms.names.shape[0])
            #cwaters = self.u.select_atoms(self.waters)
            #self.lipid_wght=np.ones(cwaters.atoms.names.shape[0])
            
            #density_z_centered += np.histogram(crds,bins=self.nbin,range=(-boxH/2,boxH/2),weights=np.ones(self.u.atoms.names.shape[0]))[0]
            #density_z_no_center += np.histogram(crds_no_center,bins=self.nbin,range=(0,boxH),weights=np.ones(self.u.atoms.names.shape[0]))[0]
            #density_lipids_center += np.histogram(clipids_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=np.ones(clipids.atoms.names.shape[0]))[0] # ELECTRON WEIGHTS 
            #density_waters_center += np.histogram(cwaters_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=np.ones(cwaters.atoms.names.shape[0]))[0] # ELECTRON WEIGHTS

            
            #print(self.wght)
            
        print("Calculating the density takes {:10.6f} s".format(time.time()-start_time))


        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        """ Normalizing the profiles """
        density_z_centered /= (frame) 
        density_z_no_center /= frame
        density_lipids_center /= frame
        density_waters_center /= frame
        

        """ Symmetrizing profile if necessary """
        #if args.symmetrize :
        #    density_z_centered += density_z_centered[::-1]
        #    density_z_centered /=2


        #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        """ Post-processign data and writing to file """
        density_data = np.vstack((x,density_z_centered)).transpose()
        density_data_no_center = np.vstack((x,density_z_no_center)).transpose()
        density_lipids_center = np.vstack((x,density_lipids_center)).transpose()
        density_waters_center = np.vstack((x,density_waters_center)).transpose()        
        
        """Post-processing of FF data from individual runs""" # running FF
        
        #fa_run = np.average(fa,axis=0)
        #fb_run = np.average(fb,axis=0)
        #fa_err = np.std(fa,axis=0)/np.sqrt(frame+1)
        #fb_err = np.std(fb,axis=0) /np.sqrt(frame+1)
        
        #fourrier_result= np.sqrt(np.multiply(fa_run,fa_run)+np.multiply(fb_run,fb_run))
        #calculate error (fa*fa_err + fb*fb_err)/ sqrt(fa^2+fb^2)
        #fourrier_error= np.multiply(1/(np.multiply(fa_run,fa_err)+np.multiply(fb_run,fb_err)),fourrier_result)
        #fourrier_data= np.vstack((FF_range*0.1*0.01,fourrier_result,fourrier_error)).transpose()
        
    
        
        
       
        """Get the indexes of the final density data where all the time steps contribute
        In other words, take the coordinates of the smalest box from the simulation"""
        final_FF_start=int(np.round(self.nbin/2-min_z/d/2))+1
        final_FF_end=int(np.round(self.nbin/2+min_z/d/2))-1
        
  
        FF_range = np.linspace(0,999,1000)
        fa_aver, fb_aver = self.fourier(density_data[final_FF_start:final_FF_end,1],density_data[final_FF_end,0]-density_data[final_FF_start,0],FF_range,density_data[1,0]-density_data[0,0])
        
        """Plot density profiles from the average density with minimal box"""
        #self.plot_density(density_data[final_FF_start:final_FF_end,:])
        #self.plot_density(density_data)
        fourrier_result2= np.sqrt(np.multiply(fa_aver,fa_aver)+np.multiply(fb_aver,fb_aver))
        fourrier_data2 = np.vstack((FF_range*0.1*0.01,fourrier_result2)).transpose()
        #self.plot_fourier_final_run(fourrier_data,fourrier_data2)
        
        #print(fourrier_data)
        #print(fourrier_data2)
        
        """Save data into files"""
        #minimum box size density
        with open(str(self.output)+"TotalDensity.dat", 'wb') as f:
            np.savetxt(f, density_data[final_FF_start+1:final_FF_end-1,:],fmt='%8.4f  %.8f')
            
        with open(str(self.output)+"LipidDensity_no_center.dat", 'wb') as f:
            np.savetxt(f, density_data_no_center[final_FF_start+1:final_FF_end-1,:],fmt='%8.4f  %.8f')

        with open(str(self.output)+"LipidDensity.dat", 'wb') as f:
            np.savetxt(f, density_lipids_center[final_FF_start+1:final_FF_end-1,:],fmt='%8.4f  %.8f')

        with open(str(self.output)+"WaterDensity.dat", 'wb') as f:
            np.savetxt(f, density_waters_center[final_FF_start+1:final_FF_end-1,:],fmt='%8.4f  %.8f')

            
            
        # density of the whole thing
        #with open(self.output, 'wb') as f:
        #    np.savetxt(f, density_data,fmt='%8.4f  %.8f')
            
        #with open(str(self.output)+".fourierFromEveryFrame", 'wb') as f:
        #    np.savetxt(f, fourrier_data,fmt='%8.4f  %.8f %.8f')
        
        #this is the important file form factors
        with open(str(self.output)+"FormFactor.dat", 'wb') as f:
            np.savetxt(f, fourrier_data2,fmt='%8.4f  %.8f')
            
        """ write output in json """
        
        with open(str(self.output)+"TotalDensity.json", 'w') as f:
            json.dump(density_data[final_FF_start+1:final_FF_end-1,:],f, cls=NumpyArrayEncoder)

        with open(str(self.output)+"LipidDensity.json", 'w') as f:
            json.dump(density_lipids_center[final_FF_start+1:final_FF_end-1,:],f, cls=NumpyArrayEncoder)

        with open(str(self.output)+"WaterDensity.json", 'w') as f:
            json.dump(density_waters_center[final_FF_start+1:final_FF_end-1,:],f, cls=NumpyArrayEncoder)
                     
        with open(str(self.output)+"FormFactor.json", 'w') as f:
            json.dump(fourrier_data2,f,cls=NumpyArrayEncoder)
            

  
        
        
    def plot_density(self,data):
        #data=np.loadtxt(self.output)
        plt.figure(figsize=(15, 6))
        plt.plot(data[:,0],data[:,1])
        plt.xlabel("Membrane normal [nm]")
        plt.show()
        

    
    def plot_fourier_final(self,data):
        plt.figure(figsize=(15, 6))
        plt.plot(data[:,0],data[:,1])
        plt.xlabel("q [A]")
        plt.show()
        
    def plot_fourier_final_run(self,data,data2):
        plt.figure(figsize=(15, 6))
        plt.plot(data[:,0],data[:,1])
        plt.errorbar(data[:,0],data[:,1],data[:,2])
        plt.plot(data2[:,0],data2[:,1])
        plt.xlabel("q [A]")
        plt.show()
        
    
    def fourier(self,ff_density,box_z,FF_range,d_ff):
        """Calculates fourier transform of ff_density in the FF_range"""
        #calculate a "height" of a bin for FF puroposes; in this case the number of bins is constant and the 
        #bin width changes
        
        
        """Creates the direct space coordinates"""
        #the calculations are stable with rounding (and others) errors in the direct space coordinates
        ff_x = np.linspace(-box_z/2,box_z/2,ff_density.shape[0]+1)[:-1] + box_z/(2*ff_density.shape[0])
        ff_x2 = np.linspace(-box_z/2,box_z/2,ff_density.shape[0]+1)[:-1] + d_ff/2
             
        
        k=0
        bulk=0
        while k*d_ff<0.33:
            bulk+=ff_density[k]+ff_density[-k-1]
            #print("The densities from the left and right ends are: {} {}".format(ff_density[k],ff_density[-k-1]))
            k+=1
        bulk/=(2*k)
        #print("The bulk density is {}".format(bulk))

            
        fa=np.zeros(FF_range.shape[0])
        fb=np.zeros(FF_range.shape[0])

        
        for j in range (0,ff_density.shape[0]):
            fa+=(ff_density[j]-bulk)*np.cos(FF_range*ff_x[j]*0.01)*d_ff
            fb+=(ff_density[j]-bulk)*np.sin(FF_range*ff_x[j]*0.01)*d_ff
            
 

        return fa, fb
 
               


# In[84]:


#new format
#FormFactor('pops_k.gro',"pops_k.xtc",200,5,'pops_k_test_density',"resname POPS","name P","number")
#FormFactor('prod.tpr',"centered.xtc",200,5,'new_test_frame40')

