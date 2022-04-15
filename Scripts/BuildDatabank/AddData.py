# This is the code that generates databank indexing 

#IMPORTING LIBRARIES

import sys
import importlib
import re
from random import randint
import argparse
import yaml

from datetime import date

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

import json
import sys

#for building hydrogens to united atom simulations 


import openmm_parser

#import databank dictionaries
from databankLibrary import lipids_dict, molecules_dict, molecule_ff_dict, gromacs_dict, amber_dict, namd_dict, charmm_dict, openmm_dict, software_dict


# Download link
from databankLibrary import download_link


#parse input yaml file
parser = argparse.ArgumentParser(description="")
parser.add_argument("-f","--file", help="Input file in yaml "
                        "format.")
args = parser.parse_args()
input_path = "./" + args.file

#load input yaml file into empty dictionary
sim = {}

#open input file for reading and writing
with open(input_path) as yaml_file:
    sim = yaml.load(yaml_file, Loader=yaml.FullLoader)
yaml_file.close()

# Show the input read
print("\n Input read from " + input_path + " file:")
print(yaml.dump(sim))

# Working directory
dir_wrk = sim['DIR_WRK']

all_molecules = []
for key in lipids_dict:
    all_molecules.append(key)
for key in molecules_dict:
    all_molecules.append(key)


# Checking that the DOI link is valid

DOI_url = 'https://doi.org/' + sim['DOI']
print("Data will be downloaded from: " + DOI_url)

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


# ### Check software used by the simulation

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

wrong_key_entries = 0
software_dict_name = "{0}_dict".format(sim['SOFTWARE'].lower())
#print(sim.items())
for key_sim, value_sim in sim.items():
        #print(key_sim, value_sim)
        #print(key_sim.upper())
    if key_sim.upper() in ("SOFTWARE"):
            #print("NOT REQUIRED")
        continue
    #Anne: check if key is in molecules_dict, molecule_numbers_dict or molecule_ff_dict too
    if (key_sim.upper() not in software_dict[sim['SOFTWARE'].upper()].keys()) and (key_sim.upper() not in molecules_dict.keys()) and (key_sim.upper() not in lipids_dict.keys()) and (key_sim.upper() not in molecule_ff_dict.keys()):
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


# PLEASE CLARIFY THIS COMMENT
# ### Process entries with files information to contain file names in arrays

software_sim = software_dict[sim['SOFTWARE'].upper()]
for key_sim, value_sim in sim.items():
    try:
        entry_type = software_sim[key_sim]['TYPE']
        if "file" in entry_type:
            if isinstance(value_sim, list): continue  
            files_list = []
            #print(value_sim)
            #print("{0} will be downloaded".format(value_sim))
            print(value_sim + " will be downloaded")
            # Place filenames into arrays
            for file_provided in value_sim.split(";"):
                files_list.append([file_provided.strip()])
            sim[key_sim] = files_list
    except: #It is notmal that fails for "ID" and "SOFTWARE"
        continue
#print(sims_files_to_array)
#print(sims_valid_entries)


# PLEASE CLARIFY THIS COMMENT
# ### Check for multiple files in entries that can only contain one

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
    print("Files are ok")
#        sims_valid_file_entries.append(sim.copy())
#print(sims_valid_file_entries)


# PLEASE CLARIFY THIS COMMENT
# ### Check if the submitted simulation has rssion has all required files and information


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
    print("All required dictionary entries are present.\n")



# ### Check status links

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

print("Starting to download data from " + DOI_url)

socket.setdefaulttimeout(15)

download_failed = False

# Create temporary directory where to download files and analyze them
dir_tmp = os.path.join(dir_wrk, "tmp_6-" + str(randint(100000, 999999)))

print("The data will be processed in directory path " + dir_tmp)

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
        #print("entry_type = {0}".format(entry_type))
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


# ## Calculate hash of downloaded files


#dir_tmp = os.path.join(dir_wrk, "tmp/")
sim_hashes = deepcopy(sim)

software_sim = software_dict[sim['SOFTWARE'].upper()]
    
#list_containing the sha1 sums for all required files
sha1_list_requied = []
    
# Make empty dataframe with the desired columns
df_files = pd.DataFrame(columns=['NAME','TYPE','REQUIRED','HASH'],dtype=object)
    
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

print("\n Summary of downloaded files: ")
print(df_files)
#print("\n{0}\n".format(sha1_list_requied))      

# Calculate the hash of a file contaning the hashes of each of the required files
# This should be always invariant as it will be used unique identifier for a simualtion
# Note order the hashes of the required files before calculating the hash (That means that the required files cannot change)
#print(sim_hashes)



#Anne:Read molecule numbers from tpr or gro file.
#Calculates numbers of lipid molecules in each leaflet. This is done by checking on which side of the centre 
#of mass the membrane each the centre of mass of a lipid molecule is.
#If a lipid molecule is split so that headgroup and tails are their own residues, the centre of mass of the
#headgroup is used in the calculation.
################################################################################################################

print("\n Calculating the numbers of lipid molecules in each leaflet based on the center of mass of the membrane and lipids. \n If a lipid molecule is split to multiple residues, the centre of mass of the headgroup is used.")

top = ''
traj = ''

# OTHER SOFTWARES THAN GROMACS!!!!
if sim['SOFTWARE'] == 'gromacs':
    top = str(dir_tmp) + '/' + sim['TPR'][0][0]
    traj = str(dir_tmp) + '/' + sim['TRJ'][0][0]
elif sim['SOFTWARE'] == 'openMM':
    traj = str(dir_tmp) + '/' + sim['TRJ'][0][0]
    top = str(dir_tmp) + '/' + sim['PDB'][0][0]
    


leaflet1 = 0 #total number of lipids in upper leaflet
leaflet2 = 0 #total number of lipids in lower leaflet
    
#u = Universe(top, traj)
#u.atoms.write(dir_tmp+'/frame0.gro', frames=u.trajectory[[0]]) #write first frame into gro file

gro = str(dir_tmp) + '/frame0.gro'

try:
    u = Universe(top, traj)
    u.atoms.write(gro, frames=u.trajectory[[0]]) #write first frame into gro file
except:
    #conf = str(dir_tmp) + '/conf.gro'
    print("Generating frame0.gro with Gromacs because MDAnalysis cannot read tpr version")
    os.system('echo System | gmx trjconv -s '+ top + ' -f '+ traj + ' -dump 0 -o ' + gro)
    u = Universe(gro, traj)
    u.atoms.write(gro, frames=u.trajectory[[0]]) #write first frame into gro file


try:
    groFORu0 = str(dir_tmp) + '/' + sim['GRO'][0][0]
    print(groFORu0)
except:
    groFORu0 = gro
    
u0 = Universe(groFORu0)
lipids = []

# select lipids 
for key_mol in lipids_dict:
    print("Calculating number of " + key_mol + " lipids")
    selection = ""
    if key_mol in sim['COMPOSITION'].keys():
       m_file = sim['COMPOSITION'][key_mol]['MAPPING']
       mapping_dict = {}
       with open('./mapping_files/'+m_file,"r") as yaml_file:
           mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
       yaml_file.close()
       for key in mapping_dict.keys():
           if 'RESIDUE' in mapping_dict[key].keys():
               selection = selection + "(resname " + mapping_dict[key]['RESIDUE'] + " and name " + mapping_dict[key]['ATOMNAME'] + ") or "
           else:      
               selection = "resname " + sim['COMPOSITION'][key_mol]['NAME']
               break
#       with open('./mapping_files/'+m_file,"r") as f:
#           for line in f:
#               if len(line.split()) > 2 and "Individual atoms" not in line:
#                   selection = selection + "(resname " + line.split()[2] + " and name " + line.split()[1] + ") or "
#               elif "Individual atoms" in line:
#                   continue
#               else:
#                   selection = "resname " + sim['COMPOSITION'][key_mol]['NAME']
#                   #print(selection)
#                   break
    selection = selection.rstrip(' or ')
    #print("selection    " + selection)
    molecules = u0.select_atoms(selection)
    #print("molecules")
    #print(molecules)
    if molecules.n_residues > 0:
        lipids.append(u0.select_atoms(selection))
        #print(lipids) 
# join all the selected the lipids together to make a selection of the entire membrane and calculate the
# z component of the centre of mass of the membrane
membrane = u0.select_atoms("")
R_membrane_z = 0
if lipids!= []:
    for i in range(0,len(lipids)):
        membrane = membrane + lipids[i]
        #print("membrane") 
        #print(membrane)  
    R_membrane_z = membrane.center_of_mass()[2]
print("Center of the mass of the membrane " + str(R_membrane_z))
    
#####number of each lipid per leaflet
        
for key_mol in lipids_dict:
    leaflet1 = 0 
    leaflet2 = 0 
        
    selection = ""
    if key_mol in sim['COMPOSITION'].keys():
        m_file = sim['COMPOSITION'][key_mol]['MAPPING']
        with open('./mapping_files/'+m_file,"r") as yaml_file:
           mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
        yaml_file.close()
        for key in mapping_dict.keys():
           if 'RESIDUE' in mapping_dict[key].keys():
               selection = selection + "resname " + mapping_dict[key]['RESIDUE'] + " and name " + mapping_dict[key]['ATOMNAME'] + " or "
           else:      
               selection = "resname " + sim['COMPOSITION'][key_mol]['NAME']
               break
        
#        with open('./mapping_files/'+m_file,"r") as f:
#            for line in f:
#                if len(line.split()) > 2 and "Individual atoms" not in line:
#                    selection = selection + "resname " + line.split()[2] + " and name " + line.split()[1] + " or "
#                elif "Individual atoms" in line:
#                    continue
#                else:
#                    selection = "resname " + sim['COMPOSITION'][key_mol]['NAME']
#                    break
    selection = selection.rstrip(' or ')
    print(selection)
    molecules = u0.select_atoms(selection)
    print(molecules.residues)

    if molecules.n_residues > 0:
        for mol in molecules.residues:
            R = mol.atoms.center_of_mass()
                
            if R[2] - R_membrane_z > 0:
                leaflet1 = leaflet1 + 1
                # print('layer1  ' + str(leaflet1))
            elif R[2] - R_membrane_z < 0:
                leaflet2 = leaflet2 +1
                  # print('layer2  ' + str(leaflet2))
    try:              
        sim['COMPOSITION'][key_mol]['COUNT'] = [leaflet1, leaflet2] 
    except KeyError:
        continue
    else:
        print("Number of " + key_mol  + " in upper leaflet: " + str(leaflet1))
        print("Number of " + key_mol  + " in lower leaflet: " + str(leaflet2))

###########################################################################################        
#numbers of other molecules
for key_mol in molecules_dict:
    try:
        mol_name = sim['COMPOSITION'][key_mol]['NAME']
    except KeyError:
        continue
    else:
        mol_number = u0.select_atoms("resname " + mol_name).n_residues
        sim['COMPOSITION'][key_mol]['COUNT'] = mol_number
        print("Number of " + key_mol  + ": " + str(sim['COMPOSITION'][key_mol]['COUNT']))   

#Anne: Read trajectory size and length 

sim['TRAJECTORY_SIZE'] = os.path.getsize(traj)

dt = 0
nsteps = 0
nstxout = 0

Nframes=len(u.trajectory)
timestep = u.trajectory.dt
 
trj_length = Nframes * timestep
   
sim['TRJLENGTH'] = trj_length

#Read temperature from tpr
if sim['SOFTWARE'] == 'gromacs':
    file1 = str(dir_tmp) + '/tpr.txt'

    print("Exporting information with gmx dump")                         #need to get temperature from trajectory not tpr !!!
    os.system('echo System | gmx dump -s '+ top + ' > '+file1)

    with open(file1, 'rt') as tpr_info:
        for line in tpr_info:
            if 'ref-t' in line:
                sim['TEMPERATURE']=float(line.split()[1])
#read temperature from xml or inp
elif sim['SOFTWARE'] == 'openMM':
#Use parser written by batuhan to read inp and xml files
    for key in ['INP','XML']:
        try:
            file1 = str(dir_tmp) + '/' + sim[key][0][0]
        except KeyError:
            print(key + ' file does not exist')
            continue
        else:
            type = key.lower()
            sim['TEMPERATURE'] = openmm_parser.openmmParser(file1,type).temperature
            break

print("Parameters read from input files:")
print("TEMPERATURE: " + str(sim['TEMPERATURE']))
print("LENGTH OF THE TRAJECTORY: " + str(sim['TRJLENGTH']))


## Check that the number of atoms between data and README.yaml match

number_of_atomsTRJ = len(u.atoms)

number_of_atoms = 0
for key_mol in all_molecules:
    mapping_dict = {}
    try:
        mapping_file = './mapping_files/'+sim['COMPOSITION'][key_mol]['MAPPING']
    except:
        continue
    else:
        with open(mapping_file,"r") as yaml_file:
           mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
        yaml_file.close()
    if sim.get('UNITEDATOM_DICT') and not 'SOL' in key_mol:
        mapping_file_length = 0
        
        for key in mapping_dict.keys():
            if 'H' in key:
                continue
            else:
                mapping_file_length += 1
    else:
        mapping_file_length = len(mapping_dict.keys())
         
    try: 
        number_of_atoms += np.sum(sim['COMPOSITION'][key_mol]['COUNT']) * mapping_file_length
    except:
        continue
#    if sim.get('UNITEDATOM_DICT') and not 'SOL' in key_mol:
#        lines = open(mapping_file).readlines(  )
#        mapping_file_length = 0
#        for line in lines:
#            if 'H' in line:
#                continue
#            else:
#                mapping_file_length += 1
#    else:
#        mapping_file_length = len(open(mapping_file).readlines(  ))
#    try:
#        number_of_atoms += np.sum(sim['COMPOSITION'][key_mol]['COUNT']) * mapping_file_length
#    except:
#        continue
        

if number_of_atoms != number_of_atomsTRJ:
    stop =  input("Number of atoms in trajectory (" +str(number_of_atomsTRJ) + ") and README.yaml (" + str(number_of_atoms) +") do no match. Check the mapping files and molecule names. \n If you know what you are doing, you can still continue the running the script. Do you want to (y/n)?")
    if stop == "n":
        os._exit("Interrupted because atomnumbers did not match")
    if stop == "y":
        print("Progressed even thought that atom numbers did not match. CHECK RESULTS MANUALLY!")

sim['NUMBER_OF_ATOMS'] = number_of_atomsTRJ
print("Number of atoms in the system: " + str(sim['NUMBER_OF_ATOMS']))


#####DATE OF RUNNING#####
today = date.today().strftime("%d/%m/%Y")
#print(today)
sim['DATEOFRUNNING'] = today

print("Date of adding to the databank: " + sim['DATEOFRUNNING'])

# Type of system is currently hard coded because only lipid bilayers are currently added.
# When we go for other systems, this will be given by user.
sim['TYPEOFSYSTEM'] = 'lipid bilayer'

# BATUHAN: add openmm parser #
# # Save to databank


# Batuhan: Creating a nested directory structure as discussed on the Issue here https://github.com/NMRLipids/NMRlipidsVIpolarizableFFs/issues/3
    
if sim['SOFTWARE'] == 'gromacs':
    head_dir = sim_hashes.get('TPR')[0][1][0:3]
    sub_dir1 = sim_hashes.get('TPR')[0][1][3:6]
    sub_dir2 = sim_hashes.get('TPR')[0][1]
    sub_dir3 = sim_hashes.get('TRJ')[0][1] 
elif sim['SOFTWARE'] == 'openMM':    
    head_dir = sim_hashes.get('TRJ')[0][1][0:3]
    sub_dir1 = sim_hashes.get('TRJ')[0][1][3:6]
    sub_dir2 = sim_hashes.get('TRJ')[0][1]
    sub_dir3 = sim_hashes.get('TRJ')[0][1]
    
print("Creating databank directories.")

os.system('mkdir ../../Data/Simulations/' + str(head_dir))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2))
os.system('mkdir ../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2) + '/' + str(sub_dir3))
    
DATAdir = '../../Data/Simulations/' + str(head_dir) + '/' + str(sub_dir1) + '/' + str(sub_dir2) + '/' + str(sub_dir3)
#    data_directory[str(ID)] = DATAdir
 
 #copy simulation trajectory and top files to DATAdir
os.system('cp '+ traj + ' ' + DATAdir)
os.system('cp '+ top + ' ' + DATAdir) 

# dictionary saved in yaml format
outfileDICT=str(dir_tmp)+ '/README.yaml'

with open(outfileDICT, 'w') as f:
    yaml.dump(sim,f, sort_keys=False)
       
    os.system('cp ' + str(dir_tmp) + '/README.yaml ' + DATAdir)
 #   outfileDICT.write(str(sim))
#outfileDICT.close()
   
print('\033[1m' + "\n Writing the README.yaml dictionary to " + DATAdir + "\n" + '\033[0m')





