
# coding: utf-8

# # Database test

# ## Simulation description file

# In[1]:
import sys
import importlib
import re
from random import randint
import argparse
import yaml

from datetime import date

# ## General Imports

# In[2]:


# Working with files and directories
import os

#For quering webs
import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

# From time monitoring
from tqdm import tqdm

import socket

# Python program to find SHA256 hash string of a file
import hashlib

# For dealing with excel and cvs 
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.max_colwidth', 1000)

#To make real independent copies of lists
from copy import deepcopy

import MDAnalysis
from MDAnalysis import Universe

#for calculating order parameters
from OrderParameter import *
import warnings
#from corrtimes import *
import subprocess
import mdtraj
import json
import sys

#for building hydrogens to united atom simulations 
import buildH_calcOP_test


# Download link
def download_link(doi, file):
    if "zenodo" in doi.lower():
        zenodo_entry_number = doi.split(".")[2]
        return 'https://zenodo.org/record/' + zenodo_entry_number + '/files/' + file
    else:
        print ("DOI provided: {0}".format(doi))
        print ("Repository not validated. Please upload the data for example to zenodo.org")
        return ""

#parse input yaml file
parser = argparse.ArgumentParser(description="")
parser.add_argument("-f","--file", help="Input file in yaml "
                        "format.")
args = parser.parse_args()
input_path = "./" + args.file

sim = {}

#open input file for reading and writing
with open(input_path) as yaml_file:
    sim = yaml.load(yaml_file, Loader=yaml.FullLoader)
yaml_file.close()

print(sim)


# ## Directories

# In[3]:


#dir_wrk = os.path.normpath(dir_wrk)
dir_wrk = sim['DIR_WRK']
#mapping_dir =os.path.normpath(mapping_dir)
#print(dir_wrk)
#print(mapping_dir)


# ## Check that the DOI link is valid

# In[4]:


DOI_url = 'https://doi.org/' + sim['DOI']
print(DOI_url)

try:
    response = urllib.request.urlopen(DOI_url)
    print("Status of the DOI link: {0}".format(response.msg))
except HTTPError as e:
    print(DOI_url)
    print('The server couldn\'t fulfill the request.')
    print('Error code: ', e.code)
    user_information = ""
    print('The code will not proceed, please fix DOI')
except URLError as e:
    print(DOI_url)
    print('We failed to reach a server.')
    print('Reason: ', e.reason)
    user_information = ""
    print('The code will not proceed, please fix DOI')
else:
    pass


# ## Read input description

# In[5]:


#bNMRLIPIDS = False #Check if the link contains NMRLIPIDS metadata
#nsims =0 # Counter number of simulations in a submission
#sims = [] #Array with the dictionary containing the information of a simulation

#for line in inFile.user_information.split("\n"):
#    if line.strip() == "":
#        continue
#    if "#NMRLIPIDS BEGIN" in line:
#        bNMRLIPIDS = True
#        continue
#    if "#NMRLIPIDS END" in line:
#        bNMRLIPIDS = False
#        continue
#    if "@SIM" in line:
#        #sims.append({"ID" : nsims, "STATUS" : 0})
#        sims.append({"ID" : nsims})
#        nsims += 1
#        sims[-1]["DOI"] = inFile.DOI
#        continue
#    if not bNMRLIPIDS:
#        continue
#    if line.strip()[0] == "@":
#        key, value = line.split("=")
#        sims[-1][key.strip('@')] = value
#print(nsims)
#print(sims)      
        


# ### Dictionaries

# In[6]:


#Molecules dictionary

lipids_dict = {
            'POPC' : {"REQUIRED": False,
                             "TYPE": "string",
                         },
            'POPG' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'POPS' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'POPE' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CHOL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'DHMDMAB' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                }

molecules_dict = {
                
            'POT' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SOD' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CLA' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'CAL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
            'SOL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                }


               
molecule_numbers_dict = {
            'NPOPC' : {"REQUIRED": False,
                              "TYPE": "array",
                          },
            'NPOPG' : {"REQUIRED": False,
                            "TYPE" : "array",
                         },
            'NPOPS' : {"REQUIRED": False,
                            "TYPE" : "array",
                        },
            'NPOPE' : {"REQUIRED": False,
                       "TYPE" : "array",
                        },
    
            'NCHOL' : {"REQUIRED": False,
                       "TYPE" : "array",
                     },
            'DHMDMAB' : {"REQUIRED": False,
                            "TYPE" : "array",
                        },
    
            'NPOT' : {"REQUIRED": False,
                            "TYPE" : "integer",
                        },
                'NSOD' : {"REQUIRED": False,
                            "TYPE" : "integer",
                        },
                'NCLA' : {"REQUIRED": False,
                            "TYPE" : "integer",
                        },
                'NCAL' : {"REQUIRED": False,
                             "TYPE" : "integer",
                         },
                'NSOL' : {"REQIRED": False,
                            "TYPE" : "integer",
                        },
    
                }

molecule_ff_dict = {
                'FFPOPC' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
                'FFPOPG' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
                'FFPOPS' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
                'FFPOPE' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
                'FFCHOL' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
                'FFDHMDMAB' : {"REQUIRED": False,
                                "TYPE": "string",
                           },
		'FFPOT' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                'FFSOD' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                'FFCLA' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                'FFCAL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
                'FFSOL' : {"REQUIRED": False,
                            "TYPE" : "string",
                        },
               }

                
                

# Gromacs
gromacs_dict = {
               'INI' : {"REQUIRED": False,
                        "TYPE" : "files",
                        "EXTENSION" : ("gro", "pdb",),
                       }, # Could be not needed in the future (tpr)
               'MDP' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("mdp",),
                       }, # Could be not needed in the future (tpr)
               'TRJ' : {"REQUIRED": True,
                        "TYPE" : "files",
                        "EXTENSION" : ("xtc","trr",),
                       },
               'TPR' : {"REQUIRED": True,
                        "TYPE" : "file",
                        "EXTENSION" : ("tpr",),
                       },
               'CPT' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("cpt",),
                       },
               'TOP' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("top",),
                       },
               'ITP' : {"REQUIRED": False,
                        "TYPE" : "files",
                        "EXTENSION" : ("itp",),
                       },
               'LOG' : {"REQUIRED": False,
                        "TYPE": "file",
                        "EXTENSION" :("log",),
                       },
               'FF'  : {"REQUIRED": False,
                        "TYPE" : "string",
                       },
               'FF_SOURCE' : {"REQUIRED": False,
                              "TYPE" : "string",
                              },
               'FF_DATE' : {"REQUIRED": False,
                            "TYPE" : "date",
                           },
               'DOI' : {"REQUIRED": True,
                            "TYPE" : "string",
                           },
    
               'SYSTEM' : {"REQUIRED": True,
                            "TYPE" : "string",
                           },
            'TEMPERATURE' : {"REQUIRED": False,
                            "TYPE" : "integer",
                            },
             'TRJLENGTH' : {"REQUIRED": False,
                           "TYPE" : "integer",
                           },
            'PREEQTIME' : {"REQUIRED": True,
                          "TYPE" : "integer",
                          },
          'TIMELEFTOUT' : {"REQUIRED":True,
                          "TYPE" : "integer",
                          },
            'UNITEDATOM' : {"REQUIRED": False,
                            "TYPE" : "string",
                         },
            'PUBLICATION' : {"REQUIRED": False,
                             "TYPE" : "string",
                            },
            'AUTHORS_CONTACT' : {"REQUIRED": False,
                                 "TYPE": "string",
                                },
            'SOFTWARE_VERSION' : {"REQUIRED": False,
                                  "TYPE": "string",
                             },
            'MAPPING_DICT' : {"REQUIRED": True,
                               "TYPE" : "dictionary",
                             },
            'DATEOFRUNNING' : {"REQUIRED": False,
                               "TYPE" : "string",
                              },
             'DIR_WRK' : {"REQUIRED": True,
                           "TYPE": "string",
                          },
               }

# Amber
amber_dict = {}

# NAMD
namd_dict = {   
            'TRJ' : { "REQUIRED": True,
                      "TYPE": "files",
                      "EXTENSION": ("dcd"),
                    },
            'INP' : { "REQUIRED": False,
                      "TYPE": "file",
                      "EXTENSION": (".inp"),
                    },
            'LOG' : { "REQUIRED": False,
                      "TYPE": "files",
                      "EXTENSION": ("log"),
                      # can be parsed to get software version etc.
                    },
            'PSF' : { "REQUIRED": False,
                      "TYPE": "file",
                      "EXTENSION": ("psf"),
                    },
            'FF'  :  { "REQUIRED": False,
                      "TYPE" : "string",
                    },
            'FF_SOURCE' : {"REQUIRED": False,
                           "TYPE" : "string",
                              },
            'FF_DATE' : {"REQUIRED": False,
                         "TYPE" : "date",
                        },
            'PDB'  : { "REQUIRED": True,
                    "TYPE": "file",
                    "EXTENSION": "pdb",}
               }
          
# CHARMM
charmm_dict = {}

# OPENMM
openmm_dict = {
               'INI' : {"REQUIRED": False,
                        "TYPE" : "files",
                        "EXTENSION" : ("gro", "pdb",),
                       }, # Could be not needed in the future (tpr)
               'MDP' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("mdp",),
                       }, # Could be not needed in the future (tpr)
               'TRJ' : {"REQUIRED": True,
                        "TYPE" : "files",
                        "EXTENSION" : ("xtc","trr",),
                       },
               'PDB' : {"REQUIRED": True,
                        "TYPE" : "file",
                        "EXTENSION" : ("pdb",),
                       },
               'TPR' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("tpr",),
                       },
               'CPT' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("cpt",),
                       },
               'TOP' : {"REQUIRED": False,
                        "TYPE" : "file",
                        "EXTENSION" : ("top",),
                       },
               'ITP' : {"REQUIRED": False,
                        "TYPE" : "files",
                        "EXTENSION" : ("itp",),
                       },
               'FF'  : {"REQUIRED": False,
                        "TYPE" : "string",
                       },
               'FF_SOURCE' : {"REQUIRED": False,
                              "TYPE" : "string",
                              },
               'FF_DATE' : {"REQUIRED": False,
                            "TYPE" : "date",
                           },
               'DOI' : {"REQUIRED": True,
                            "TYPE" : "string",
                           },

               'SYSTEM' : {"REQUIRED": True,
                            "TYPE" : "string",
                           },
            'TEMPERATURE' : {"REQUIRED": False,
                            "TYPE" : "integer",
                            },
             'TRJLENGTH' : {"REQUIRED": False,
                           "TYPE" : "integer",
                           },
            'PREEQTIME' : {"REQUIRED": True,
                          "TYPE" : "integer",
                          },
          'TIMELEFTOUT' : {"REQUIRED":True,
                          "TYPE" : "integer",
                          },
            'MAPPING' : {"REQUIRED": True,
                             "TYPE" : "string",
 #                        "EXTENSION": ("txt"),
                             },

               }

# SOFTWARE
software_dict = {
                "GROMACS" : gromacs_dict, 
                "AMBER"   : amber_dict,
                "NAMD"    : namd_dict,
                "CHARMM"  : charmm_dict,
                "OPENMM"  : openmm_dict,
                }

#print(software_dict.keys())


# ### Check software used by the simulation

# In[7]:


#sims_valid_software = []

if sim['SOFTWARE'].upper() in software_dict.keys():
    msg_info = "Simulation uses supported software {0} and will be further processed"
        #print(msg_info.format(sim['ID'], sim['SOFTWARE'].upper()))
    print(msg_info.format(sim['SOFTWARE'].upper()))
#        sims_valid_software.append(sim.copy())
else:
    msg_err="Simulation performed in an UNSUPPORTED software {0} and will NOT be further processed"
    print(msg_err.format(sim["SOFTWARE"].upper()))
    quit()
        
#print(sims_valid_software) 


# ### Check that all entry keys provided for each simulation are valid:

# In[8]:

#print(gromacs_dict.keys())


#sims_valid_entries = []
#for sim in sims_valid_software:
    #print("ID {0}".format(sim["ID"]))
wrong_key_entries = 0
software_dict_name = "{0}_dict".format(sim['SOFTWARE'].lower())
print(sim.items())
for key_sim, value_sim in sim.items():
        #print(key_sim, value_sim)
        #print(key_sim.upper())
    if key_sim.upper() in ("SOFTWARE"):
            #print("NOT REQUIRED")
        continue
    #Anne: check if key is in molecules_dict, molecule_numbers_dict or molecule_ff_dict too
    if (key_sim.upper() not in software_dict[sim['SOFTWARE'].upper()].keys()) and (key_sim.upper() not in molecules_dict.keys()) and (key_sim.upper() not in lipids_dict.keys()) and (key_sim.upper() not in molecule_numbers_dict.keys()) and (key_sim.upper() not in molecule_ff_dict.keys()):
        print ("{0} NOT in {1}".format(key_sim, software_dict_name)) 
        wrong_key_entries += 1
if wrong_key_entries:
    print("Simulation has {0} unknown entry/ies and won't be longer considered, please correct.\n".format(wrong_key_entries))
    quit()
else:
    msg_info = "All entries in simulation are understood and will be further processed\n"
    print(msg_info)
#        sims_valid_entries.append(sim.copy())
#print(sims_valid_entries)


# ### Process entries with files information to contain file names in arrays

# In[9]:


#sims_files_to_array = deepcopy(sims_valid_entries)

#for sim in sims_files_to_array:
 #   print("ID {0}".format(sim["ID"]), flush=True)
software_sim = software_dict[sim['SOFTWARE'].upper()]
for key_sim, value_sim in sim.items():
    try:
        entry_type = software_sim[key_sim]['TYPE']
        if "file" in entry_type:
            if isinstance(value_sim, list): continue  
            files_list = []
            print("{0} added to list".format(value_sim))
            # Place filenames into arrays
            for file_provided in value_sim.split(";"):
                files_list.append([file_provided.strip()])
            sim[key_sim] = files_list
    except: #It is notmal that fails for "ID" and "SOFTWARE"
        continue
#print(sims_files_to_array)
#print(sims_valid_entries)


# ### Check for multiple files in entries that can only contain one

# In[10]:


#sims_valid_file_entries = []
#for sim in sims_files_to_array:
#    print("ID {0}".format(sim["ID"]), flush=True)
files_issues = 0
software_sim = software_dict[sim['SOFTWARE'].upper()]
for key_sim, value_sim in sim.items():
    try:
        entry_type = software_sim[key_sim]['TYPE']
        if entry_type == "file"  and len(value_sim) > 1:
            print("Multiple values found in {0} and only one allowed (Please correct):\n {1}".format(key_sim,value_sim))
            files_issues += 1
    except: #It is notmal that fails for "ID" and "SOFTWARE"
        continue
if files_issues:
    print("Sim will be no longer processed")
    quit()
else:
    print("ok")
#        sims_valid_file_entries.append(sim.copy())
#print(sims_valid_file_entries)


# ### Check if the submitted simulation has rssion has all required files and information

# In[11]:


missing_required_keys = 0
for key, value in software_dict[sim['SOFTWARE'].upper()].items():
    if value["REQUIRED"]:
        try:
            sim[key]
        except:
            print("Entry not found: {0} {1}".format(key, value))
            missing_required_keys += 1
if missing_required_keys:
    print("{0} missing required entry/ies, please correct.".format(missing_required_keys))
    print("Entry will not be further processed.\n")
    quit()
else:
    print("All required entries present.\n")
    #sims_required_entries.append(sim.copy())


# ### Check status links

# In[12]:

# In[13]:      


wrong_links = 0
software_sim = software_dict[sim['SOFTWARE'].upper()]
for key_sim, value_sim in sim.items():
    #print("key_sim = {0} => value_sim = {1}".format(key_sim, value_sim))
    try:
        entry_type = software_sim[key_sim]['TYPE']
        extension_type = software_sim[key_sim]['EXTENSION']
        #print("entry_type = {0}".format(entry_type))
        if "file" in entry_type and "txt" not in extension_type:
            for file_provided in value_sim:
                #print("File={0}".format(file_provided[0]))
                file_url = download_link(DOI, file_provided[0])
                if file_url == "":
                        
                    wrong_links += 1
                    continue
                try:
                    if key_sim == 'INI':
                        continue
                    response = urllib.request.urlopen(file_url)
                    #print("Status of the DOI link: {0}".format(response.msg))
                except HTTPError as e:
                    print("\nkey={0} => file={1}".format(key_sim, file_provided[0]))
                    print(file_url)
                    print('The server couldn\'t fulfill the request.')
                    print('Error code: ', e.code)
                    wrong_links += 1
                except URLError as e:
                    print(key_sim, file_provided[0])
                    print(file_url)
                    print('We failed to reach a server.')
                    print('Reason: ', e.reason)
                    wrong_links += 1
                else:
                    pass
    except: #It is notmal that fails for "ID" and "SOFTWARE"
        continue
if wrong_links:
    print("{0} link/s failed, please correct.".format(wrong_links))
    print("Entry will not be further processed.\n")
    quit()
else:
     print("All links work.\n")
     #sims_working_links.append(sim.copy())
#print(sims_working_links)


# ## Download files from links

# In[14]:

socket.setdefaulttimeout(30)

download_failed = False

# Create temporary directory where to download files and analyze them
dir_tmp = os.path.join(dir_wrk, "tmp_6-" + str(randint(100000, 999999)))

if (not os.path.isdir(dir_tmp)): 
    os.mkdir(dir_tmp)

file_sizes={}
software_sim = software_dict[sim['SOFTWARE'].upper()]
dir_sim = dir_tmp
DOI = sim['DOI']
if (not os.path.isdir(dir_sim)): 
    os.mkdir(dir_sim)
for key_sim, value_sim in sim.items():
    #print("key_sim = {0} => value_sim = {1}".format(key_sim, value_sim))
    try:
        entry_type = software_sim[key_sim]['TYPE']
        extension_type = software_sim[key_sim]['EXTENSION']
        print("entry_type = {0}".format(entry_type))
        if "file" in entry_type and "txt" not in extension_type:
            for file_provided in tqdm(value_sim, desc = key_sim):
                file_url = download_link(DOI, file_provided[0])
                file_name = os.path.join(dir_sim, file_provided[0])
                #get the size of the file to be downloaded
                url_size = urllib.request.urlopen(file_url).length
                if (not os.path.isfile(file_name)):
                    print("downloading")
                    try:
                        response = urllib.request.urlretrieve(file_url, file_name)
                    except ContentTooShortError:
                        download_failed = True
                        print("Content too short error.")
                    except HTTPError as e:
                        download_failed = True
                        print(e)
                    except URLError as ue:
                        download_failed = True
                        print("failed to download")
                    except socket.timeout as se:
                        download_failed = True
                        print("socket time out")
                    except Exception as ee:
                        download_failed = True
                        print(ee)
                #check if the file is fully downloaded
                size = os.path.getsize(file_name)
                print("size of the file " + file_provided[0] + " to be downloaded: " + str(url_size))
                print("size of the file " + file_provided[0] + " after download: " + str(size) )
                if url_size != size:
                    print("Download of the file " + file_provided[0] + " was interrupted.")
                    quit()
    except:#It is normal that fails for "ID" and "SOFTWARE"
        continue

if download_failed:
    print("One of the downloads failed. Terminating the script.")
    quit()

#print(sim)



# ## Calculate hash downloaded files

# In[15]:


#dir_tmp = os.path.join(dir_wrk, "tmp/")
sim_hashes = deepcopy(sim)

#for sim in sims_hashes:
# print("ID {0}".format(sim["ID"]), flush=True)
software_sim = software_dict[sim['SOFTWARE'].upper()]
# dir_sim = os.path.join(dir_tmp, str(sim["ID"])) 
    
#list_containing the sha1 sums for all required files
sha1_list_requied = []
    
# Make empty dataframe with the desired columns
df_files = pd.DataFrame(columns=['NAME','TYPE','REQUIRED','HASH'])
    
for key_sim, value_sim in sim_hashes.items():
        #print("key_sim = {0} => value_sim = {1}".format(key_sim, value_sim))
    try:
        entry_type = software_sim[key_sim]['TYPE']
            #print("entry_type = {0}".format(entry_type))
        if "file" in entry_type:
            files_list = []
            for file_provided in value_sim:
                file_name = os.path.join(dir_sim, file_provided[0]) 
                sha1_hash = hashlib.sha1()
                with open(file_name,"rb") as f:
                        # Read and update hash string value in blocks of 4K
                    for byte_block in iter(lambda: f.read(4096),b""):
                        sha1_hash.update(byte_block)
                        #print(file_provided, sha256_hash.hexdigest())
                    df_files = df_files.append({
                        "NAME":file_provided[0],
                        "TYPE":key_sim,
                        "REQUIRED": software_dict[sim_hashes['SOFTWARE'].upper()][key_sim]['REQUIRED'],
                        "HASH":sha1_hash.hexdigest(),
                    }, ignore_index=True)
                files_list.append([file_provided[0], sha1_hash.hexdigest()])
                #Find the keys of the required files to calculate the master_hash 
            if software_dict[sim_hashes['SOFTWARE'].upper()][key_sim]['REQUIRED'] == True:
                sha1_list_requied.append(sha1_hash.hexdigest())
            sim_hashes[key_sim] = files_list #Problematic
    except: #It is notmal that fails for "ID" and "SOFTWARE"
        continue
print(df_files)
print("\n{0}\n".format(sha1_list_requied))      
# Calculate the hash of a file contaning the hashes of each of the required files
# This should be always invariant as it will be used unique identifier for a simualtion
# Note order the hashes of the required files before calculating the hash (That means that the required files cannot change)
print(sim_hashes)


# In[16]:

#in case of a united atom simulation make a dictionary of united atom names 
if sim.get('UNITEDATOM'):
    unitedAtoms = sim['UNITEDATOM'].split(',')
    unitedAtomsDic = {}
    for i in range(0, int(len(unitedAtoms)/2)):
        lipid = unitedAtoms[2*i]
        UAlipid = unitedAtoms[2*i+1]
        unitedAtomsDic[lipid]=UAlipid
    sim['UADICTIONARY'] = unitedAtomsDic
        
#    print(sim)

# ## Read molecule numbers into dictionary

# In[22]:

#Anne:Read molecule numbers from tpr or gro file.
#Calculates numbers of lipid molecules in each leaflet. This is done by checking on which side of the centre 
#of mass the membrane each the centre of mass of a lipid molecule is.
#If a lipid molecule is split so that headgroup and tails are their own residues, the centre of mass of the
#headgroup is used in the calculation.
################################################################################################################

tpr = str(dir_tmp) + '/' + str(sim.get('TPR')).translate({ord(c): None for c in "']["})
trj = str(dir_tmp) + '/' + str(sim.get('TRJ')).translate({ord(c): None for c in "']["})

# OTHER SOFTWARES THAN GROMACS!!!!

#if sim['SOFTWARE'] == gromacs:
gro = str(dir_tmp) + '/conf.gro'
    #make gro file
os.system('echo System | gmx trjconv -f '+trj+' -s '+tpr+' -dump 0 -o ' +gro)

    # add gro into dictionary for later use
sim['GRO'] = gro
    
leaflet1 = 0 #total number of lipids in upper leaflet
leaflet2 = 0 #total number of lipids in lower leaflet
    
u = Universe(gro)
lipids = []

# select lipids 
for key_mol in lipids_dict:
    print(key_mol)
    selection = ""
    if key_mol in sim['MAPPING_DICT'].keys():
        m_file = sim['MAPPING_DICT'][key_mol]
        with open('./mapping_files/'+m_file,"r") as f:
            for line in f:
                if len(line.split()) > 2 and "Individual atoms" not in line:
                    selection = selection + "(resname " + line.split()[2] + " and name " + line.split()[1] + ") or "
                elif "Individual atoms" in line:
                    continue
                else:
                    selection = "resname " + sim.get(key_mol)
                    print(selection)
                    break
    selection = selection.rstrip(' or ')
    #print("selection    " + selection)
    molecules = u.select_atoms(selection)
    print("molecules")
    print(molecules)
    if molecules.n_residues > 0:
        lipids.append(u.select_atoms(selection))
        print(lipids) 
# join all the selected the lipids together to make a selection of the entire membrane and calculate the
# z component of the centre of mass of the membrane
membrane = u.select_atoms("")
R_membrane_z = 0
if lipids!= []:
    for i in range(0,len(lipids)):
        membrane = membrane + lipids[i]
    print("membrane") 
    print(membrane)  
    R_membrane_z = membrane.center_of_mass()[2]
print(R_membrane_z)
    
#####number of each lipid per leaflet
        
for key_mol in lipids_dict:
    leaflet1 = 0 
    leaflet2 = 0 
        
    selection = ""
    if key_mol in sim['MAPPING_DICT'].keys():
        m_file = sim['MAPPING_DICT'][key_mol]
        with open('./mapping_files/'+m_file,"r") as f:
            for line in f:
                if len(line.split()) > 2 and "Individual atoms" not in line:
                    selection = selection + "resname " + line.split()[2] + " and name " + line.split()[1] + " or "
                elif "Individual atoms" in line:
                    continue
                else:
                    selection = "resname " + sim.get(key_mol)
                    break
    selection = selection.rstrip(' or ')
#   print(selection)
    molecules = u.select_atoms(selection)
    print(molecules.residues)
    x = 'N' + key_mol
    if molecules.n_residues > 0:
        for mol in molecules.residues:
            R = mol.atoms.center_of_mass()
                
            if R[2] - R_membrane_z > 0:
                leaflet1 = leaflet1 + 1
              # print('layer1  ' + str(leaflet1))
            elif R[2] - R_membrane_z < 0:
                leaflet2 = leaflet2 +1
              # print('layer2  ' + str(leaflet2))
    sim[x] = [leaflet1, leaflet2] 

        
#print("upper leaflet: " + str(leaflet1))
#print("lower leaflet: " + str(leaflet2))

###########################################################################################        
#numbers of other molecules
for key_mol in molecules_dict:
    value_mol = sim.get(key_mol)
    print(value_mol)
    x = 'N' + key_mol
    sim[x] = u.select_atoms("resname " + value_mol ).n_residues
        

# In[23]:


#Anne: Read temperature and trajectory length from tpr file

dt = 0
nsteps = 0
nstxout = 0

file1 = str(dir_tmp) + '/tpr.txt'

os.system('echo System | gmx dump -s '+ tpr + ' > '+file1)
    
    

with open(file1, 'rt') as tpr_info:
    for line in tpr_info:
        if 'ref-t' in line:
            sim['TEMPERATURE']=line.split()[1]
    
mol = Universe(tpr, trj)
Nframes=len(mol.trajectory)
timestep = mol.trajectory.dt
 
trj_length = Nframes * timestep
   
sim['TRJLENGTH'] = trj_length
    
print(sim)

#####DATE OF RUNNING#####
today = date.today().strftime("%d/%m/%Y")
print(today)
sim['DATEOFRUNNING'] = today

print(sim)


# # Save to databank

# In[24]:

# Batuhan: Creating a nested directory structure as discussed on the Issue here https://github.com/NMRLipids/NMRlipidsVIpolarizableFFs/issues/3
    
head_dir = sim_hashes.get('TPR')[0][1][0:3]
sub_dir1 = sim_hashes.get('TPR')[0][1][3:6]
sub_dir2 = sim_hashes.get('TPR')[0][1]
sub_dir3 = sim_hashes.get('TRJ')[0][1]
    
os.system('mkdir ../../Data/Simulations/' + str(head_dir))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2) + '/' + str(sub_dir3))
    
DATAdir = '../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2) + '/' + str(sub_dir3)
#    data_directory[str(ID)] = DATAdir
    
    # SAMULI: I am writin now in txt, but using json might be better in the future
   # outfileDICT=open(str(DATAdir)+'/README.md','w')
   # outfileDICT=str(dir_wrk)+'/tmp/'+str(ID)+'/'+ 'README.json'
    
   # with open(outfileDICT, 'w') as f:
   #     json.dump(sim, f)
   

    #!cp {str(dir_wrk)}'/tmp/'{str(ID)}'/README.json' {data_directory.get(str(ID))}  
# dictionary saved in yaml format
outfileDICT=str(dir_tmp)+ '/README.yaml'
    
    #with open(outfileDICT, 'w') as f:
    #    yaml.dump(sim,f, sort_keys=False)
       
    #    os.system('cp ' + str(dir_tmp) + '/README.yaml ' + DATAdir)
    #outfileDICT.write(str(sim))
    #    outfileDICT.close()
   


# # Analysis starts here
# 

# In[ ]:


#!cp corr_ftios_ind.sh {dir_wrk}
#for sim in sims_working_links:
trj=sim.get('TRJ')
tpr=sim.get('TPR')
   # ID=sim.get('ID')
software=sim.get('SOFTWARE')
EQtime=float(sim.get('TIMELEFTOUT'))*1000
unitedAtom = sim.get('UNITEDATOM')
    
ext=str(trj)[-6:-3] # getting the trajectory extension
print(ext)
    # BATUHAN: Adding a few lines to convert the trajectory into .xtc using MDTRAJ
    #          We will need users to install MDTRAJ in their system so that we can convert other trajectories into xtc

if ext != "xtc":
        
    print("converting the trajectory into xtc")
        
    pdb = sim.get('PDB')
    output_traj = str(dir_tmp) + '/' + 'tmp_converted.xtc'
    input_traj = str(dir_tmp) + '/' + trj[0][0]
    input_pdb = str(dir_tmp) + '/' + pdb[0][0]
      
    if os.path.isfile(output_traj): # when we're done with the converted trajectory we can simply remove it
        os.system('rm {output_traj}')
        
    os.system('echo System | mdconvert {input_traj} -o {output_traj} -t {input_pdb} --force # force overwrite')
        
        # SAMULI: this xtcwhole does not necessarily have molecules as whole. Only if {input_traj} has.
    xtcwhole = str(dir_tmp) + '/' + 'tmp_converted.xtc'
    tpr=input_pdb
        
    print("trajectory conversion is completed")
        
else:
    
    xtc = str(dir_tmp) + '/' + str(trj[0][0])  
    tpr = str(dir_tmp) + '/' + str(tpr[0][0])
    xtcwhole=str(dir_tmp) + '/whole.xtc'
    os.system('echo System | gmx trjconv -f ' + xtc + ' -s ' + tpr + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))
   

    
print("Calculating order parameters")
    
if unitedAtom:
    for key in sim['UADICTIONARY']:
    #construct order parameter definition file for CH bonds from mapping file
        def_file = open(str(dir_tmp) + '/' + key + '.def', 'w')

        mapping_file = sim['MAPPING_DICT'][key]
        previous_line = ""
            
        with open('./mapping_files/'+mapping_file, "r") as f:
            for line in f.readlines():
                if not line.startswith("#"):
                    regexp1_H = re.compile(r'M_[A-Z0-9]*C[0-9]*H[0-9]*_M')
                    regexp2_H = re.compile(r'M_G[0-9]*H[0-9]*_M')
                    regexp1_C = re.compile(r'M_[A-Z0-9]*C[0-9]*_M')
                    regexp2_C = re.compile(r'M_G[0-9]_M')

                    if regexp1_C.search(line) or regexp2_C.search(line):
                        atomC = line.split()
                        atomH = []
                    elif regexp1_H.search(line) or regexp2_H.search(line):
                        atomH = line.split()
                    else:
                        atomC = []
                        atomH = []

                    if atomH:
                        items = [atomC[1], atomH[1], atomC[0], atomH[0]]
                        def_line = items[3] + " " + key + " " + items[0] + " " + items[1] + "\n"
                        if def_line != previous_line:
                            def_file.write(def_line)
                            print(def_line)
                            previous_line = def_line
        def_file.close()

     #Add hydrogens to trajectory and calculate order parameters with buildH
        ordPfile = str(dir_tmp) + '/' + key + '_order_parameter.dat' 
        topfile = sim.get('GRO')
        deffile = str(dir_tmp) + '/' + key + '.def' 
        lipidname = sim['UADICTIONARY'][key]
        #    print(lipidname)
        buildH_calcOP_test.main(topfile,lipidname,deffile,xtcwhole,ordPfile)
else:
    for key in sim['MAPPING_DICT']:    
        mapping_file = sim['MAPPING_DICT'][key]
        OrdParam=find_OP('./mapping_files/'+mapping_file,tpr,xtcwhole,key)

        outfile=open(str(dir_tmp) + '/' + key + 'OrderParameters.dat','w')
        line1="Atom     Average OP     OP stem"+'\n'
        outfile.write(line1)
    
        data = {}
        outfile2=str(dir_tmp) + '/' + key + 'OrderParameters.json' 
        
        for i,op in enumerate(OrdParam):
            resops =op.get_op_res
            (op.avg, op.std, op.stem) =op.get_avg_std_stem_OP
            line2=str(op.name)+" "+str(op.avg)+" "+str(op.stem)+'\n'
            outfile.write(line2)
    
            data[str(op.name)]=[]
            data[str(op.name)].append(op.get_avg_std_stem_OP)
        
        with open(outfile2, 'w') as f:
            json.dump(data,f)

        outfile.close()

        os.system('cp ' + str(dir_tmp) + '/' + key + 'OrderParameters.dat ' + DATAdir)
        os.system('cp ' +str(dir_tmp) + '/' + key + 'OrderParameters.json ' + DATAdir)
    
print("Done calculating order parameters.")


