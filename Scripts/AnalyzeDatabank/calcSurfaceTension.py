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
            tpr = system.get('TPR')
            edr_name = systemPath + system.get('EDR')[0][0]
            edr_url = download_link(doi, edr[0][0])
            tpr_name = systemPath + system.get('TPR')[0][0]
            tpr_url = download_link(doi, tpr[0][0])
            
            # Download the EDR file if not available
            if (not os.path.isfile(edr_name)):
                print('Downloading the EDR file to ', system['path'])
                response = urllib.request.urlretrieve(edr_url, edr_name)
                
            # Download the TPR file if not available
            if (not os.path.isfile(tpr_name)):
                print('Downloading the TPR file to ', system['path'])
                response = urllib.request.urlretrieve(tpr_url, tpr_name)
                
            # Generate a PDB file from the TPR
            os.system( f'gmx editconf -f {tpr_name} -o {systemPath}conf.pdb >/dev/null 2>&1' )
            
            # Get the size of the box in the Z dimension
            u = MDAnalysis.Universe( tpr_name, systemPath + 'conf.pdb' )
            Lz = u.dimensions[2]
            
            edata = MDAnalysis.auxiliary.EDR.EDRReader( edr_name )
            
            frames = edata.data_dict["Time"] >= system['TIMELEFTOUT']*1000
            
            Pxx = edata.data_dict['Pres-XX'][ frames ]
            Pyy = edata.data_dict['Pres-YY'][ frames ]
            Pzz = edata.data_dict['Pres-ZZ'][ frames ]
            
            surft = { time: value for time, value in zip( edata.data_dict["Time"][ frames ],
                                                          Lz*(Pzz-0.5*(Pxx+Pyy))/2 * 0.01 ) } 
        
            with open(outfilename, 'w') as f:
                json.dump(surft,f)
                
            os.remove( systemPath + 'conf.pdb' )
        
        except:
            print ('It was not possible to compute the surface tension of the system.')
