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

from databankLibrary import lipids_dict


# In[2]:

# to write data in numpy arrays into json file
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)
        
def getLipids(readme, molecules=lipids_dict.keys()):
    lipids = 'resname '

    for key in readme['COMPOSITION'].keys():
        if key in molecules:
            mapping_dict = {}
            mapping_file = readme['COMPOSITION'][key]['MAPPING']
            with open('./mapping_files/'+mapping_file,"r") as yaml_file:
                mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
            yaml_file.close()
            
            switch = 0
            for key in mapping_dict:
                try:
                    res = mapping_dict[key]['RESIDUE']
                except KeyError:
                    switch = 1
                    continue
                else:
                    if res not in lipids:
                        lipids = lipids + res + ' or resname '
            
            if switch == 1:
                lipids = lipids + readme['COMPOSITION'][key]['NAME'] + ' or resname '
            
            
 #           with open('./mapping_files/'+m_file,"r") as f:
 #               for line in f:
 #                   if len(line.split()) > 2 and "Individual atoms" not in line:
 #                       if line.split()[2] not in lipids:
 #                           lipids = lipids + line.split()[2] + ' or resname '
 #                   elif "Individual atoms" in line:
 #                       continue
 #                   else:
 #                       lipids = lipids + readme['COMPOSITION'][key]['NAME'] + ' or resname '
                        break    
    lipids = lipids[:-12]
    print(lipids)
    return lipids

def getWater(readme, molecules=lipids_dict.keys()):
    waters = 'resname ' + readme['COMPOSITION']['SOL']['NAME']
    return waters


#for testing purposes to safe time for loading trajectory and creating dictonary
#the electron.dat will be decripted in the future as the values will be used from the universal mapping file

def temporary_mapping_dictionary(readme):
    mapping_dictonary={}
    for key1 in readme['COMPOSITION'].keys(): #
        
        key2 = readme['COMPOSITION'][key]['NAME']
        mapping_dict ={}
    
        with open('./mapping_files/'+mapping_file,"r") as yaml_file:
            mapping_dict = yaml.load(yaml_file, Loader=yaml.FullLoader)
        yaml_file.close()
        
        # put mapping files into a larger dictionary where molecule name is key and the 
        # dictionary in the correponding mapping file is the value
        mapping_dictionary[key2]= mapping_dict #try if this works
        
    return mapping_dictonary


electron_dictionary={"H":1,"He":2,"Li":3,"Be":4,"B":5,"C":6,"N":7,"O":8,"F":9,"Ne":10,
                   "Na":11,"Mg":12,"Al":13,"Si":14,"P":15,"S":16,"Cl":17,"Ar":18,
                   "K":19,"Ca":20,"Sc":21,"Ti":22,"V":23,"Cr":24,"Mn":25,"Fe":26,"Co":27,"Ni":28,
                     "Cu":29,"Zn":30,"Ga":31,"Ge":32,"Vi":0}
	


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
            self.u = mda.Universe(self.conf,self.traj)
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
        
        self.calculate_weight()
        

        self.calculate_density()
        
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
            #with open ('electrons.dat') as f:
            #    rawdata = f.read().split('\n')
             #   lines= rawdata[0:len(rawdata)-1]
             #   tmp_atoms = np.asarray([l.split()[0:2] for l in lines])
             #   electrons = dict([ (elem[0],int(elem[1])) for elem in tmp_atoms])
        
            #weights to calculate the electron density
            wght=np.zeros(self.u.atoms.names.shape[0])
           ########3 
          
                    
            ##########
            #Assign electron weights to atoms of the system
#RICKY
#            for i in range(0,self.u.atoms.names.shape[0]):
#                name = self.u.atoms.names[i] #names of atoms (names in gro or similar file)
#                residue = self.u.atoms.resnames[i] #residue names (names in gro or similar file)
#                temporary_mapping = temporary_mapping_dictionary(self.readme) #dictionary made of mapping files
#                #try:
#                #
#                m_name = temporary_mapping[][] 
#                #except KeyError:
#                #    #print(m_name)
#                #    print(residue)
#                #    print(name)
#                #else:
#                name1 = re.sub(r'M_[0-9]*','',m_name[::-1])
#                name2 = re.sub(r'M_([A-Z]{1,2}[0-9]{1,4})*','',name1[::-1]) #name of the atom extracted from mapping name
#                if name2 == 'G':
#                    name2 = 'C'
#                #print(name2)                
#                wght[i]=electron_dictionary[name2] #get number of electrons
#               # print(wght[i])
#            self.wght=wght
#        if self.density_type=="number":
#            self.wght=np.ones(self.u.atoms.names.shape[0])
#        if self.density_type=="mass":
#            self.wght=self.u.atoms.masses
            
#######ANNE
        # combine all mapping files of the system to a single nested dictionary
            temporary_mapping = temporary_mapping_dictionary(self.readme) #dictionary made of mapping files
        
            for molecule, mapping_name in temporary_mapping:
                name1 = re.sub(r'M_[0-9]*','',mapping_name[::-1]) # removes numbers and '_M' from the end of string and reverses the string
                name2 = re.sub(r'M_([A-Z]{1,2}[0-9]{1,4})*','',name1[::-1]) # name2 is the atom and electrons are assigned to this atom
            
                if name2 == 'G': # G is carbon so change G to C
                    name2 = 'C'
                
                #write number of electrons to the mapping dictionary of the whole system 
                #temporary_mapping[molecule][mapping_name]['ELECTRONS'] = electron_dictionary[name2]
                wght[i]=electron_dictionary[name2] #get number of electrons
                
        if self.density_type=="number":
            self.wght=np.ones(self.u.atoms.names.shape[0])
        if self.density_type=="mass":
            self.wght=self.u.atoms.masses           
        
        print(self.lipids)
        print(self.u.atoms.resnames)
            
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
            lipid_wght=np.ones(c.atoms.names.shape[0])
            
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

            
            """shif the coordinates so that the center in z-dimention is in 0; 
            #divide by 10 to get the coordinates in nm, since now the crds are only the z coordinates"""
            crds = (self.u.atoms.positions[:,2] - box_z/2)/10
            clipids_center = (clipids.atoms.positions[:,2] - box_z/2)/10
            cwaters_center = (cwaters.atoms.positions[:,2] - box_z/2)/10



            water_wght=np.ones(cwaters.atoms.names.shape[0])
            
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
          
 
         

            """calculates the total density profile; keep for now"""
            density_z_centered += np.histogram(crds,bins=self.nbin,range=(-boxH/2,boxH/2),weights=self.wght/vbin)[0]
            density_z_no_center += np.histogram(crds_no_center,bins=self.nbin,range=(0,boxH),weights=lipid_wght/vbin)[0]
            density_lipids_center += np.histogram(clipids_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=lipid_wght/vbin)[0]
            density_waters_center += np.histogram(cwaters_center,bins=self.nbin,range=(-boxH/2,boxH/2),weights=water_wght/vbin)[0]

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

