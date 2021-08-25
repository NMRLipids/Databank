# code to generate transition matrix
import yaml
import os
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis import Universe

import urllib.request
from urllib.error import URLError,HTTPError,ContentTooShortError

sys.path.insert(1, '../BuildDatabank/')
from databankLibrary import download_link, read_mapping_file, lipids_dict

def getLipids(readme):
        lipids = []
        for key in lipids_dict.keys():
            try:
                if readme['N'+key] != [0,0]: 
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
        
def average_box_z(universe, readme):
    EQtime=float(readme['TIMELEFTOUT'])*1000
    stop = len(universe.trajectory)
    dt = universe.trajectory.dt
    start = int(EQtime*dt) 
    z = []
    for ts in universe.trajectory[start:stop]:
        z_ts = universe.dimensions[2]
        z.append(z_ts)
    avg_box_z = sum(z)/len(universe.trajectory)
    return avg_box_z

def whichBin(w_molecule,bins, H):   
    w_com_z = w_molecule.center_of_mass()[2]   # z-axis coordinate
   # print(w_com_z)
    # frame i
    for i, b_i in enumerate(bins):
         #if periodic boundary is crossed translate z coordinate of the centre of mass of the water molecule by the lenght of the box 
         if H <= w_com_z/H:
             print(str(H) + ' <= ' + str(w_com_z/H))
             w_com_z = w_com_z - H
             
         if w_com_z/H < 0:
             print(str(w_com_z/H))
             w_com_z = w_com_z + H
            
         # find the index of the bin where the water molecule is
         if i == 0 and w_com_z/H < (b_i + spacing/2):
              in_bin = i
         elif i == H and (b_i - spacing / 2) <= w_com_z/H:
             in_bin = i
         else:
             if (b_i - spacing / 2) <= w_com_z/H < (b_i + spacing/2):
                 in_bin = i
             #else:
                 #print(str(b_i - spacing / 2) + ' <= ' + str(w_com_z) + ' < ' + str(b_i + spacing/2))
                 
         
                 
 #   print(str(in_bin))
    return in_bin                  
    
    
    
def membraneCentreOfMass(universe, readme):
    lipids = []
    for key_mol in lipids_dict:
        print("Calculating number of " + key_mol + " lipids")
        selection = ""
        if key_mol in readme['MAPPING_DICT'].keys():
            m_file = readme['MAPPING_DICT'][key_mol]
            with open('./mapping_files/'+m_file,"r") as f:
                for line in f:
                    if len(line.split()) > 2 and "Individual atoms" not in line:
                        selection = selection + "(resname " + line.split()[2] + " and name " + line.split()[1] + ") or "
                    elif "Individual atoms" in line:
                       continue
                    else:
                        selection = "resname " + sim.get(key_mol)
                        #print(selection)
                        break
        selection = selection.rstrip(' or ')
        #print("selection    " + selection)
        molecules = universe.select_atoms(selection)
        #print("molecules")
        #print(molecules)
        if molecules.n_residues > 0:
            lipids.append(universe.select_atoms(selection))
            #print(lipids) 
    # join all the selected the lipids together to make a selection of the entire membrane and calculate the
    # z component of the centre of mass of the membrane
    membrane = universe.select_atoms("")
    R_membrane_z = 0
    if lipids!= []:
        for i in range(0,len(lipids)):
            membrane = membrane + lipids[i]
            #print("membrane") 
            #print(membrane)  
    R_membrane = membrane.center_of_mass()
    return R_membrane
    


# thickness of the membrane calculated based on the average distance of phosphate group to membrane centre of mass
def membraneThickness(universe,readme,membraneCoM):
    lipids = getLipids(readme)
    all_phosphates = []
    for lipid in lipids:
        phosphate1 = ['M_G3O1_M', 'M_G3P2_M', 'M_G3P2O1_M', 'M_G3P2O2_M', 'M_G3O3_M']
        phosphate2 = []
        mapping_file = '../BuildDatabank/mapping_files/'+readme['MAPPING_DICT'][lipid]
        #translate atom names
     #   for i,a in enumerate(phosphate1):
     #       phosphate2[i] = read_mapping_file(mapping_file, a)
        atom1, atom2, atom3, atom4, atom5 = [read_mapping_file(mapping_file, a) for a in phosphate1]
        selection = u.select_atoms("resname " + lipid + " and name {} {} {} {} {}".format(atom1,atom2,atom3,atom4,atom5))
        all_phosphates.append(selection)
    
    distances = []    
    for sel in all_phosphates:
        for group in sel:
            dist = np.abs(group.center_of_mass()[2] - R_membrane[2])
            distances.append(dist)
    h = np.sum(distances) / len(distances)
    
    return h*2
        
        
        
#normal to the centre of mass of the bilayer ?????
#z = 0
 
#average length of box along z 




    

# [b_i - spacing/2, b_i + spacing/2)

# centre of mass of a water molecule moves from [b_i - spacing/2, b_i + spacing/2) to [b_j - spacing/2, b_j + spacing/2)

# matrix T_k contains number of water molecules in each bin b_i in frame k
# matrix T_k+1 contains number of water molecules in each bin b



##################################################################################################################################################################
# make directory for saving transition matrices
matrix_path = "../../Data/TransitionMatrix/"
os.system('mkdir ' + matrix_path)
        
for subdir, dirs, files in os.walk(r'../../Data/Simulations/f40/bb6/f40bb6ab5d44402be07059e8df74b5a8200f031e/6774168dfec0a5a7377c8a46341eba603f320cf7'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readme = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[4:8])
               # print(indexingPath)
                #download files to
                sub_dirs = indexingPath.split("/")
                os.system('mkdir ../../Data/TransitionMatrix/' + sub_dirs[0])
                os.system('mkdir ../../Data/TransitionMatrix/' + sub_dirs[0] + '/' + sub_dirs[1])
                os.system('mkdir ../../Data/TransitionMatrix/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2])    
                os.system('mkdir ../../Data/TransitionMatrix/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2] + '/' + sub_dirs[3])
                simulationDir = "../../Data/Simulations/" + indexingPath
                dataDIR = matrix_path + indexingPath
                 
                doi = readme['DOI']
                trj=readme['TRJ'][0][0]
                tpr=readme['TPR'][0][0]
                trj_name = simulationDir + '/' + trj
                tpr_name = simulationDir + '/' + tpr
                trj_url = download_link(doi, trj)
                tpr_url = download_link(doi, tpr)
                
                #Download tpr and xtc files to a tmp directory in dir wrk defined in readme file / or should these be put into the same directory as README.yaml ???
                if (not os.path.isfile(tpr_name)):
                    response = urllib.request.urlretrieve(tpr_url, tpr_name)
                        
                if (not os.path.isfile(trj_name)):
                    response = urllib.request.urlretrieve(trj_url, trj_name)

                EQtime=float(readme['TIMELEFTOUT'])*1000
                
                xtcwhole=str(simulationDir) + '/whole.xtc'
                if (not os.path.isfile(xtcwhole)):
                    print("Make molecules whole in the trajectory")
                    os.system('echo System | gmx trjconv -f ' + xtc + ' -s ' + tpr + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime))

                #mda universe
                u = Universe(tpr_name,xtcwhole)
                
                # length of simulation box in z-axis direction at ts 
          
                
                # simulation domain is divided into n bins HOW TO DIVIDE and HOW MANY BINS??
                # try 100 for now
                nbins = 100
                


                # list of matrices to save matrix of each frame
                matrix_list = []
                
                end = len(u.trajectory)
                dt = u.trajectory.dt
                start = int(EQtime*dt)
          
                water = u.select_atoms("resname " + readme['SOL'] )
                w_molecules = water.split('residue') 

                bins_previous = []
                
                for ts in u.trajectory[start:end]:
                    bins_fr = []
                    print('frame :' + str(ts.frame))
                    
                    # length of box in z-axis direction
                    H = u.dimensions[2]
                    # divide simulation box into n bins - 100 bins for now
                    spacing = 1 / nbins
                    b = []
                    b0 = 0
                    for i in range(0,nbins + 1):
                        b.append(b0)
                        b0 += spacing

                   # print(len(w_molecules)) 
                    transitions = np.zeros((len(w_molecules),2),dtype=int)
                    
                    for j,w in enumerate(w_molecules):  
                        bins_fr.append(whichBin(w,b,H))
                    
                    if (ts.frame == 0):
                        bins_previous = bins_fr.copy()
                    else:
                        transitions[:,0] = bins_previous    # eiiii nää on samat
                        transitions[:,1] = bins_fr
                       # print(transitions)

                        M_fr = np.zeros((nbins+1,nbins+1),dtype=int)
                        # count transitions from bin to bin per consecutive frames into matrix M_fr
                        for i in range(0,len(transitions[:,0])):
                            k = transitions[i,0]
                          #  print('k: ' + str(k))
                            l = transitions[i,1]
                          #  print('l: ' + str(l)) 
                            M_fr[k,l] += 1 
                        #print(M_fr)
                        matrix_list.append(M_fr)
                        # save the locations of water molecules of the current frame
                        bins_previous = bins_fr.copy()
                    
  
                    
                # sum matrices  
                t_matrix = matrix_list[0]            
                for m in matrix_list[1:end]:
                    A = t_matrix + m
                    t_matrix = A 
                    
                print(t_matrix)  
                      
                # matrix is written to a txt file
                matrix_file = dataDIR + "/transition_matrix_bins" + str(nbins) + ".txt"
                
                # lagtime 
                lt = str(0) #???
                #dn = ???
                avg_box_z = average_box_z(u, readme)
                bins_to_txt_file = np.linspace(-0.5*avg_box_z, 0.5*avg_box_z, num=100, endpoint=True)
                #bin edges z = 0 at the center of mass of the membrane
                edges = np.array2string(bins_to_txt_file,separator=' ')
                
                characters = '[]\n'
                edges_cleaned = edges
                for c in characters:
                    edges_cleaned = edges_cleaned.replace(c, '')
                    
                    
                header =  'lt    ' + lt + '\n' + 'dt    ' + dt + '\n' + 'dn   ??? \n' + 'edges   ' + edges_cleaned + '\n'

                np.savetxt(matrix_file,t_matrix,header=header,  fmt='%d')



             

             
                    
                    
                    
                    
                    
