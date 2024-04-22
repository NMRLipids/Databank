#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:48:05 2024

@author: fabs
"""
import os
import re
import sys
import json
import yaml
import functools 
import MDAnalysis 
import multiprocess

import numpy as np


Parallel = False
if len(sys.argv)-1:
    if sys.argv[0]:
        print( "Analysis will be parallelized" )
        Parallel = True

if __name__ == "__main__":
    
    # ! H= -3.790e-5
    neutron_dictionary={"H":-3.790e-5,"He":3.260e-5,"Li":-1.900e-5, "Be":-7.790e-5,"B":  5.300e-5, "C": 6.646e-5,"N":  9.36e-5,"O": 5.803e-5,
                        "F": 5.654e-5,"Ne":4.566e-5,"Na": 3.630e-5, "Mg": 5.375e-5,"Al": 3.499e-5, "Si":1.491e-5,"P":  5.13e-5,"S": 2.847e-5,
                        "Cl":9.577e-5,"Ar":1.909e-5,"K": 3.670e-5,  "Ca":4.700e-5,"Sc": 12.290e-5,"Ti":-34.38e-5,"V": -0.3824e-5,"Cr":3.635e-5,
                        "Mn":-3.73e-5,"Fe":9.45e-5, "Co":2.49e-5,   "Ni":1.03e-5,"Cu":7.718e-5,"Zn":5.680e-5,"Ga": 7.288e-5, "Ge": 8.185e-5,"Vi": 0, 
                        "Cs":5.420e-5,"D":  6.671e-5,"":0}

    databankPath = '../../'

    import sys
    sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
    from databankLibrary import ( initialize_databank, lipids_dict )

    systems = initialize_databank(databankPath)
    #systems = [ system for system in systems if system["path"]=="0c7/33d/0c733dce5d09ee8b1f66ba847183f7f2329ba4a2/1a60d0f1be4b69b95bcf23fbe207c20e96e2d341/" ]

    
    def getNeutrons(mapping_name):
        # Adapted from FormFactor.py
        name1 = re.sub(r'[0-9]*_M$', '', mapping_name) # removes numbers and '_M' from the end of string
        name2 = re.sub(r'^M_([A-Z]{1,2}[0-9]{1,4})*','',name1) # removes M_X12Y23... sequence from the beginning
        if name2 == 'G': # G is a glycerol carbon so change G to C
            name2 = 'C'
        try:
            n = neutron_dictionary[name2]
        except KeyError:
            print (f'ERROR: This mapping name cannot be read by our rules: {mapping_name}', file=sys.stderr)
            print ('Consider changing naming in your mapping file.')
            sys.exit(1)
        return n
    
    
    for system in systems:
        
        try:
            if "TYPEOFSYSTEM" in system.keys():
                if system["TYPEOFSYSTEM"] == "lipid monolayer":
                    
                    path = databankPath + "Data/Simulations/" + system["path"] + "/"
                    
                    outfilename = databankPath + '/Data/Simulations/' +  system['path'] + 'nsdens.json'
            
                    if os.path.isfile(outfilename):
                        continue    
                    
                    print( "Analyzing system: ", system["path"] )

            
                    u = MDAnalysis.Universe( path + system["TPR"][0][0],
                                             path + system["TRJ"][0][0],
                                             in_memory=True if Parallel else False )
                    
                    mapping = { mol: yaml.load( open( "../BuildDatabank/mapping_files/" + system["COMPOSITION"][mol]["MAPPING"], "r" ), 
                                                 Loader=yaml.FullLoader ) for mol in system["COMPOSITION"] }
                    
                    translator = { mol: { mapping[mol][j]["ATOMNAME"]: j for j in mapping[mol] } for mol in mapping.keys() }
                    
                ### Lipids
                    # List of lipids in the system
                    lipids = list( set( system["COMPOSITION"].keys() ) & set( lipids_dict.keys() ) )
                    
                    # Split leaflets
                    COG = u.select_atoms( f"resname {' '.join( lipids )}" ).center_of_geometry()[2]
                    
                    UL = u.select_atoms( f"resname {' '.join( lipids )} and prop z>{COG}" )
                    LL = u.select_atoms( f"resname {' '.join( lipids )} and prop z<{COG}" )
                    
                    # Selection with all lipids
                    Lipids = UL + LL
                    
                    # Universal names of the lipid atoms
                    Lipid_names = [ translator[lipid][name] for lipid, name in zip( Lipids.resnames, Lipids.names ) ]
                    
                    # Neutron scattering length associated to each lipid atom
                    Lipid_N = np.vectorize( getNeutrons )( Lipid_names )
                    
                    
                ### Water
                    # A selection with the water
                    Waters = u.select_atoms( f"resname {system['COMPOSITION']['SOL']['NAME']}" )
                
                    # Universal names of the water atoms
                    Water_names = [ translator["SOL"][name] for name in Waters.names  ]
                    
                    # Neutron scattering length associated to each water atom
                    Water_N = np.vectorize( getNeutrons )( Water_names * 2 )
                
                
                ### Ions
                    ions = list( set( system["COMPOSITION"].keys() ) - set( lipids + ["SOL"] ) )
                    
                    if ions:
                        # A selection with the ions
                        Ions = u.select_atoms( f"resname {' '.join(ions)}" )
                        
                        # Universal names of the ions
                        Ions_names = [ translator[ion][name] for ion, name in zip( Ions.resnames, Ions.names ) ]
                        
                        # Neutron scattering length associated to each ion
                        Ions_N = np.vectorize( getNeutrons )( Ions_names * 2 )
                        
                    
                ##### Run
                    # Skip the initial equilibration time
                    frames_to_skip = int( system['TIMELEFTOUT']*1000 / ( u.trajectory[1].time-u.trajectory[0].time ) )
                    
                    def getNSDens( ts, bins=None, first=False ):
                        
                        # Set the current timestep
                        u.trajectory[ts]
                        
                        ref = { lipid: dict(zip(translator[lipid].values(), translator[lipid].keys()))["M_G2_M"] for lipid in lipids if "M_G2_M" in translator[lipid].values() }
                        
                        # Split leaflets
                        COG = u.select_atoms( f"resname {' '.join( lipids )}" ).center_of_geometry()[2]
                        
                        UL = u.select_atoms( f"resname {' '.join( lipids )} and prop z>{COG}" )
                        LL = u.select_atoms( f"resname {' '.join( lipids )} and prop z<{COG}" )
                        
        
                        # Put the coordinates relatively to the limit of the leaflet
                        ref_UL = np.mean( UL.select_atoms( " or ".join( [ f"( resname {resname} and name {name} )" for resname, name in ref.items() ] ) ).positions[:,2] ) # min( UL.positions[:,2] )
                        ref_LL = np.mean( LL.select_atoms( " or ".join( [ f"( resname {resname} and name {name} )" for resname, name in ref.items() ] ) ).positions[:,2] ) # max( LL.positions[:,2] )
                        
                        Coords_lip = np.hstack( [  ref_UL - UL.positions[:,2],
                                                   LL.positions[:,2] - ref_LL  ] )   
                        
                        Waters = u.select_atoms( f"resname {system['COMPOSITION']['SOL']['NAME']}" )
                        
                        # Water coodinates (duplicated and relative to both leaflets)
                        Coords_wat = np.hstack( [  ref_UL - Waters.positions[:,2],
                                                   Waters.positions[:,2] - ref_LL ] )  
                        
                        # Ion coordinates
                        if ions:
                            Ions = u.select_atoms( f"resname {' '.join(ions)}" )
                            
                            Coords_ion = np.hstack( [  ref_UL - Ions.positions[:,2],
                                                       Ions.positions[:,2] - ref_LL ] )  
                        
                        # Generate a slicing of the space
                        if first:
                            bins = np.linspace( min( Coords_lip ), ( max( UL.positions[:,2] ) - COG ) + min( Coords_lip ), num=50 )
                            
                        Dens_l = np.zeros( len(bins) )
                        Dens_w = np.zeros( len(bins) )
                        Dens_i = np.zeros( len(bins) )
                                        
                
                        # Assign each atom to its slice
                        Slices_lip = np.digitize( Coords_lip, bins ) 
                        Slices_wat = np.digitize( Coords_wat, bins ) 
                        if ions:
                            Slices_ion = np.digitize( Coords_ion, bins ) 
                            
                        for i in range( len( bins ) ):
                            Dens_l[i] += sum( Lipid_N[ Slices_lip == i ] )
                            Dens_w[i] += sum( Water_N[ Slices_wat == i ] )
                            if ions: 
                                Dens_i[i] += sum( Ions_N[ Slices_ion == i ] )
                        
                        if first:
                            return Dens_l, Dens_w, Dens_i, bins
                        else:
                            return Dens_l, Dens_w, Dens_i
                    
                    
                    Dens_l, Dens_w, Dens_i, bins = getNSDens( frames_to_skip,first=True )
                    
                    if Parallel:
                        pool = multiprocess.Pool( )
                        run_per_frame = functools.partial( getNSDens, bins=bins, first=False )
                        dens = pool.map( run_per_frame, np.arange(frames_to_skip+1,len(u.trajectory)), )
                    else:
                        dens = []
                        for ts in range(frames_to_skip+1,len(u.trajectory)):
                            dens.append( (lambda x: getNSDens(x,bins,False))(ts) )
                    
                    dens = np.sum(dens,axis=0)
                    
                    Dens_l += dens[0]
                    Dens_w += dens[1]
                    Dens_i += dens[2]
                
                    # Divide by the number of frames, the volume of the slice and 2, because of counting twice
                    Dens_l /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                    Dens_w /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                    Dens_i /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                
                    # Total electron density
                    Dens = Dens_l + Dens_w + Dens_i
                    
                    # Compute XRR 
                    qz_range = np.arange(0,0.5,0.005)
                    
                    NSDensData = { z: [ dens_l, dens_w, dens_i ] for z, dens_l, dens_w, dens_i in zip( bins, Dens_l, Dens_w, Dens_i ) }
                    
                    with open( outfilename, 'w') as f:
                        json.dump( NSDensData, f )
                
                    if Parallel:
                        del u
                
        except:
            print ('It was not possible to compute the neutron scattering density of the system.')
        

