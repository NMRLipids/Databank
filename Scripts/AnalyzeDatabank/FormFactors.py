import os
import sys
import numpy as np
import json
import matplotlib.pyplot as plt
import MDAnalysis
import urllib.request
import yaml

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

# From time monitoring
from tqdm import tqdm

import socket

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, lipids_dict, databank
import form_factor


#path = '../../Data/Simulations/102/dfa/102dfa390326acd775347b1d703ccfa85e033cd8/1ad80d5d253e861f01d8f16e0d1590c6ef962916/'
path = '../../Data/Simulations/'

db_data = databank(path)
systems = db_data.get_systems()

with open("system_already_run_17_11_21.out","w") as f:
    f.write("Listing analysed systems \n \n ")   
    
     
for system in systems:
    #download trajectory and gro files
    system_path = system['path']

    print('Analyzing', system_path)

    try:
        if system['UNITEDATOM_DICT']:
            #continue
            print('United atom simulation')
    except:
        pass
    
    output_name = ""
    
    trj_name = system['path'] + system['TRJ'][0][0]
    tpr_name = system['path'] + system['TPR'][0][0]
    trj_url = download_link(system['DOI'], system['TRJ'][0][0])
    tpr_url = download_link(system['DOI'], system['TPR'][0][0])

    socket.setdefaulttimeout(15)

    download_failed = False

    
    # make a function like this
    if (not os.path.isfile(tpr_name)):
        try:
            url_size = urllib.request.urlopen(tpr_url).length
            response = urllib.request.urlretrieve(tpr_url, tpr_name)
            
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
        else:
            size = os.path.getsize(tpr_name)
            print("size of the file "  + " to be downloaded: " + str(url_size))
            print("size of the file " + " after download: " + str(size) )
            if url_size != size:
                print("Download of the file "  + " was interrupted.")
               


                        
    if (not os.path.isfile(trj_name)) and download_failed == False:
        
        try:
            url_size = urllib.request.urlopen(trj_url).length
            response = urllib.request.urlretrieve(trj_url, trj_name)
       
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
        else:
            #check if the file is fully downloaded
            size = os.path.getsize(trj_name)
            print("size of the file " +  " to be downloaded: " + str(url_size))
            print("size of the file " +  " after download: " + str(size) )
            if url_size != size:
                print("Download of the file "  + " was interrupted.")
            print("Download was a success!!!!!!")
     
    if download_failed == True:
        with open("system_already_run_17_11_21.out","a") as f:
            f.write(system_path+" Download failed !!!\n")    
        
    if download_failed == False:
        with open("system_already_run_17_11_21.out","a") as f:
            f.write(system_path+" download successfull \n")   
        if (not os.path.isfile(system_path + "/FormFactor.json")):
            try:
                form_factor.FormFactor(system_path, tpr_name, trj_name, 200, output_name,  system)
            except:
                print(system_path,' Failed')
