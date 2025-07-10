#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 14:00:26 2024

@author: fabs
TODO: Unite with createPDBs.py
"""

import os
import sys
import yaml
import urllib.request

# where to deploy PDBs
outpath_prefix = os.environ['NMLDB_DEPLOY_PATH']

FAILS=[]

if len( sys.argv ) > 1:
    files = sys.argv[1:]
    
    for file in files:
        try:
            path = f"./Databank/Data/Simulations/{file}/"
            outpath = os.path.join(outpath_prefix, file)
            if not os.path.isdir( outpath ):
                os.system( f"mkdir -p {outpath}" )
            readme = path + "README.yaml"
            data = yaml.load( open( readme, "r" ), Loader=yaml.FullLoader)
            if ( not os.path.isfile( path + data["TPR"][0][0] ) ):
                response = urllib.request.urlretrieve( "https://zenodo.org/record/" + data["DOI"].split(".")[2] + "/files/" + data["TPR"][0][0], 
                                                      path + data["TPR"][0][0] )    
                
            if 'gromacs' in data['SOFTWARE']:
                if os.path.isfile( f"{outpath}conf.pdb" ):
                    os.remove( f"{outpath}conf.pdb" )
                if os.path.isfile( f"{outpath}conf.pdb.gz" ):
                    os.remove( f"{outpath}conf.pdb.gz" )
                os.system( f'gmx editconf -f {path}{data["TPR"][0][0]} -o {outpath}conf.pdb' )
                os.remove( path + data["TPR"][0][0] )
                os.system( f"gzip {outpath}conf.pdb" )
        except:
            FAILS.append( file )
    
    if FAILS:
        print( f"It was not possible to generate a PDB for:" )
        for file in FAILS:
            print(file)
    
else:
    print("No systems provided")

