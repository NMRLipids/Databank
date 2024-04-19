#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 13:39:43 2024

@author: fabs
"""

import os
import sys
import json
import MDAnalysis
import urllib.request

### Initializing the databank

### This is the NMRlipids databank repository path
databankPath = '../../'

sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
from databankLibrary import ( initialize_databank, download_link )

systems = initialize_databank(databankPath)


for system in systems:
    if "EDR" in system.keys():
        try:
            outfilename = databankPath + '/Data/Simulations/' +  system['path'] + 'surft.json'
            
            if os.path.isfile(outfilename):
                continue
            
            print('Analyzing: ', system['path'] )
            
            systemPath = os.path.dirname(os.path.realpath(__file__)) + '/../../Data/Simulations/' + system['path'] 
            doi = system.get('DOI')
            edr = system.get('EDR')
            edr_name = systemPath + system.get('EDR')[0][0]
            edr_url = download_link(doi, edr[0][0])
            
            # Download the EDR file if not available
            if (not os.path.isfile(edr_name)):
                print('Downloading the EDR file to ', system['path'])
                response = urllib.request.urlretrieve(edr_url, edr_name)
            
            edata = MDAnalysis.auxiliary.EDR.EDRReader( edr_name )
            
            surft = { time: value for time, value in zip( edata.data_dict["Time"][ edata.data_dict["Time"] >= system['TIMELEFTOUT']*1000 ],
                                                          edata.data_dict['#Surf*SurfTen'][ edata.data_dict["Time"] >= system['TIMELEFTOUT']*1000 ]/2 ) } 
        
            with open(outfilename, 'w') as f:
                json.dump(surft,f)
        
        except:
            print ('It was not possible to compute the surface tension of the system.')
