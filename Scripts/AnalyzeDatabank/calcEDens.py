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
import traceback

import numpy as np


Parallel = False
if len(sys.argv)-1:
    if sys.argv[0]:
        print( "Analysis will be parallelized" )
        Parallel = True




if __name__ == "__main__":
    
    # It is not the same as in FormFactor.py !
    electron_dictionary={"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
                       "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,
                       "K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,
                         "Cu":29,"Zn":30,"Ga":31,"Ge":32,"Vi":0,"Cs":55,"D":0,"":0}

    databankPath = '../../'

    import sys
    sys.path.insert(1, databankPath + '/Scripts/BuildDatabank/')
    from databankLibrary import ( initialize_databank, lipids_dict )

    systems = initialize_databank(databankPath)
    
    def getElectrons(mapping_name):
        # From FormFactor.py
        name1 = re.sub(r'[0-9]*_M$', '', mapping_name) # removes numbers and '_M' from the end of string
        name2 = re.sub(r'^M_([A-Z]{1,2}[0-9]{1,4})*','',name1) # removes M_X12Y23... sequence from the beginning
        if name2 == 'G': # G is a glycerol carbon so change G to C
            name2 = 'C'
        try:
            el = electron_dictionary[name2]
        except KeyError:
            print (f'ERROR: This mapping name cannot be read by our rules: {mapping_name}', file=sys.stderr)
            print ('Consider changing naming in your mapping file.')
            sys.exit(1)
        return el
    
    
    for system in systems:
        
        try:
            if "TYPEOFSYSTEM" in system.keys():
                if system["TYPEOFSYSTEM"] == "lipid monolayer":
                    
                    path = databankPath + "Data/Simulations/" + system["path"] + "/"
                    
                    outfilename = databankPath + '/Data/Simulations/' +  system['path'] + 'edens.json'
            
                    if os.path.isfile(outfilename):
                        continue
                    
                    print( "Analyzing system: ", system["path"] )
                    
                    u = MDAnalysis.Universe( path + system["TPR"][0][0],
                                             path + system["TRJ"][0][0],
                                             in_memory=True if Parallel else False )
                    
                    mapping = {
                        mol: yaml.load(open("../BuildDatabank/mapping_files/" + 
                                            system["COMPOSITION"][mol]["MAPPING"], "r" ), 
                                       Loader=yaml.FullLoader) 
                             for mol in system["COMPOSITION"] 
                    }
                    
                    translator = {
                        mol: {mapping[mol][j]["ATOMNAME"]: j for j in mapping[mol]} 
                             for mol in mapping.keys() 
                    }
                    
                    
                ### Lipids
                    # List of lipids in the system
                    lipids = list( set( system["COMPOSITION"].keys() ) & set( lipids_dict.keys() ) )
                    print(lipids)
                    
                    # Split leaflets
                    COG = u.select_atoms( f"resname {' '.join( lipids )}" ).center_of_geometry()[2]
                    
                    UL = u.select_atoms( f"resname {' '.join( lipids )} and prop z>{COG}" )
                    LL = u.select_atoms( f"resname {' '.join( lipids )} and prop z<{COG}" )
                    
                    # Selection with all lipids
                    Lipids = UL + LL
                    
                    # resname - name pair array
                    resnm_anm_array = list(zip( Lipids.resnames, Lipids.names ))
                    # Universal names of the lipid atoms
                    Lipid_names = [ translator[lipid][name] for lipid, name in resnm_anm_array ]
                    
                    # Atomic number of the lipid atoms
                    Lipid_Z = np.vectorize( getElectrons )( Lipid_names )
                    
                    # Start United-atom processing
                    try:
                        uaDict = system["UNITEDATOM_DICT"]
                    except:
                        uaDict = None
                    if type(uaDict) is not dict:
                        uaDict = None
                    # Dictionary for adding additional electrons
                    add_el = {'CH': 1, 'CH2': 2, 'CH3': 3, 'nop': 0}
                    if uaDict is not None:
                        print("Processing United-Atom topology..")
                        for lipid in lipids:
                            if lipid not in uaDict:
                                continue
                            print(f"Reading {uaDict[lipid]}..", end='')
                            UAlipidjsonNAME = './lipid_json_buildH/' + uaDict[lipid] + '.json'
                            with open(UAlipidjsonNAME) as json_file:
                                UAlipidjson = json.load(json_file)
                            i = 0
                            print("done.")
                            for res, name in resnm_anm_array:
                                print(res,name)
                                if res == lipid:
                                    if name in UAlipidjson:
                                        ua_typ = UAlipidjson[name][0]
                                    else:
                                        ua_typ = 'nop'
                                    Lipid_Z[i] += add_el[ua_typ]
                                i += 1
                    # debug msg
                    print(Lipid_Z)
                    ### Water
                    # A selection with the water
                    Waters = u.select_atoms( f"resname {system['COMPOSITION']['SOL']['NAME']}" )
                    print(translator["SOL"]) 
                    # Universal names of the water atoms
                    Water_names = [ translator["SOL"][name] for name in Waters.names  ]
                    
                    # Atomic number of the lipid atoms
                    Water_Z = np.vectorize( getElectrons )( Water_names * 2 )


                    ### Ions
                    ions = list( set( system["COMPOSITION"].keys() ) - set( lipids + ["SOL"] ) )
                    
                    if ions:
                        # A selection with the ions
                        Ions = u.select_atoms( f"resname {' '.join(ions)}" )
                        
                        # Universal names of the ions
                        Ions_names = [ translator[ion][name] for ion, name in zip( Ions.resnames, Ions.names ) ]
                        
                        # Atomic number of the lipid atoms
                        Ions_Z = np.vectorize( getElectrons )( Ions_names * 2 )
                        
                    
                ##### Run
                    
                    # Skip the initial equilibration time
                    frames_to_skip = int( system['TIMELEFTOUT']*1000 / ( u.trajectory[1].time-u.trajectory[0].time ) )
                
                    
                    def getElDens( ts, bins=None, first=False ):
                        
                        # Set the current timestep
                        u.trajectory[ts]
                        
                        ref = { lipid: dict(zip(translator[lipid].values(), translator[lipid].keys()))["M_G2_M"] 
                                       for lipid in lipids if "M_G2_M" in translator[lipid].values() }
                        
                        # Split leaflets
                        COG = u.select_atoms( f"resname {' '.join( lipids )}" ).center_of_geometry()[2]
                        
                        UL = u.select_atoms( f"resname {' '.join( lipids )} and prop z>{COG}" )
                        LL = u.select_atoms( f"resname {' '.join( lipids )} and prop z<{COG}" )
                        
                        # Selection with all lipids
                        Lipids = UL + LL
                        
                        # Put the coordinates relatively to the limit of the leaflet
                        ref_UL = np.mean( UL.select_atoms( " or ".join( [ f"( resname {resname} and name {name} )" 
                                                           for resname, name in ref.items() ] ) ).positions[:,2] ) # min( UL.positions[:,2] )
                        ref_LL = np.mean( LL.select_atoms( " or ".join( [ f"( resname {resname} and name {name} )" 
                                                           for resname, name in ref.items() ] ) ).positions[:,2] ) # max( LL.positions[:,2] )
                        
                        Coords_lip = np.hstack( [  ref_UL - UL.positions[:,2],
                                                   LL.positions[:,2] - ref_LL  ] )   
                        
                        # The lipids will be associated to negative values of the position
                        coef = -1 if np.sign( np.mean(Coords_lip) )>0 else 1
                        Coords_lip *= coef
                        
                        Waters = u.select_atoms( f"resname {system['COMPOSITION']['SOL']['NAME']}" )
                        
                        # Water coodinates (duplicated and relative to both leaflets)
                        Coords_wat = np.hstack( [  ref_UL - Waters.positions[:,2],
                                                   Waters.positions[:,2] - ref_LL ] )  
                        Coords_wat *= coef
                        Coords_wat += u.dimensions[2]*(Coords_wat<-(max(Coords_lip)-min(Coords_lip)))
                        
                        # Ion coordinates
                        if ions:
                            Ions = u.select_atoms( f"resname {' '.join(ions)}" )
                            
                            Coords_ion = np.hstack( [  ref_UL - Ions.positions[:,2],
                                                       Ions.positions[:,2] - ref_LL ] )  
                            
                            Coords_ion *= coef
                            Coords_ion += u.dimensions[2]*(Coords_ion<-(max(Coords_lip)-min(Coords_lip)))
                        
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
                            Dens_l[i] += sum( -Lipids[ Slices_lip == i ].charges + Lipid_Z[ Slices_lip == i ] )
                            Dens_w[i] += sum( -(Waters+Waters)[ Slices_wat == i ].charges + Water_Z[ Slices_wat == i ] )
                            if ions: 
                                Dens_i[i] += sum( -(Ions+Ions)[ Slices_ion == i ].charges + Ions_Z[ Slices_ion == i ] )
                        
                        if first:
                            return Dens_l, Dens_w, Dens_i, bins
                        else:
                            return Dens_l, Dens_w, Dens_i
                    
                    
                    Dens_l, Dens_w, Dens_i, bins = getElDens( frames_to_skip,first=True )
                    
                    if Parallel:
                        pool = multiprocess.Pool( )
                        run_per_frame = functools.partial( getElDens, bins=bins, first=False )
                        dens = pool.map( run_per_frame, np.arange(frames_to_skip+1,len(u.trajectory)), )
                    else:
                        dens = []
                        for ts in range(frames_to_skip+1,len(u.trajectory)):
                            dens.append( (lambda x: getElDens(x,bins,False))(ts) )
                    
                    dens = np.sum(dens,axis=0)
                    
                    Dens_l += dens[0]
                    Dens_w += dens[1]
                    Dens_i += dens[2]
                
                    # Divide by the number of frames, the volume of the slice and 2, because of counting twice
                    Dens_l /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                    Dens_w /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                    Dens_i /= 2 * len( u.trajectory[frames_to_skip:] ) * ( u.trajectory[-1].dimensions[0] * u.trajectory[-1].dimensions[1] * (bins[1]-bins[0]) )
                    
                    EDensData = { z: [ dens_l, dens_w, dens_i ] for z, dens_l, dens_w, dens_i in zip( bins, Dens_l, Dens_w, Dens_i ) }
                    
                    with open( outfilename, 'w') as f:
                        json.dump( EDensData, f )
                    
                    if Parallel:
                        del u
                
        except Exception as e:
            print ('It was not possible to compute the electron density of the system.')
            print(traceback.format_exc())
        

