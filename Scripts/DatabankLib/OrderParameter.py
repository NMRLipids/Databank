
#!/usr/bin/env python
"""
 calculation of order parameters of lipid bilayers
 from a MD trajectory
 
 meant for use with NMRlipids projects
 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/02/02
------------------------------------------------------------
 input: Order parameter definitions
        gro and xtc file (or equivalents)
 output: order parameters (2 textfiles)
--------------------------------------------------------
"""

# coding: utf-8

import MDAnalysis as mda
import numpy as np
import math
import os, sys
import warnings
import re

from optparse import OptionParser
from collections import OrderedDict
from operator import add
from linecache import getline


#k_b = 0.0083144621  #kJ/Mol*K
#f_conc=55430  # factor for calculating concentrations of salts from numbers of ions/waters; in mM/L

bond_len_max=1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq=bond_len_max**2

#%%
class OrderParameter:
    """
    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """
    def __init__(self, resname, atom_A_name, atom_B_name, M_atom_A_name, M_atom_B_name, *args):  #removed name, resname from arguments
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
#        self.name = name             # name of the order parameter, a label
        self.resname = resname       # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        self.M_atAname = M_atom_A_name
        self.M_atBname = M_atom_B_name
        self.name = M_atom_A_name + " " + M_atom_B_name # generic name of atom A and atom B
        for field in self.__dict__:
            if not isinstance(field, str):
                warnings.warn("provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur.")#.format(field)
            else:
                if not field.strip():
                    raise RuntimeError("provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions.")#.format(field)
        # extra optional arguments allow setting avg,std values -- suitable for reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
            self.stem = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
            self.stem = None
        else:
            warnings.warn("Number of optional positional arguments is {len}, not 2 or 0. Args: {args}\nWrong file format?")
        self.traj = []  # for storing OPs
        self.selection = []
        

    def calc_OP(self, atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
       # print(atoms)
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()
	
        if d2>bond_len_max_sq:
            at1=atoms[0].name
            at2=atoms[1].name
            resnr=atoms[0].resid
            d=math.sqrt(d2)
            warnings.warn("Atomic distance for atoms \
            {at1} and {at2} in residue no. {resnr} is suspiciously \
            long: {d}!\nPBC removed???")
        cos2 = vec[2]**2/d2
        S = 0.5*(3.0*cos2-1.0)
        return S
        


    @property
    def get_avg_std_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))
    @property
    def get_avg_std_stem_OP(self):
        """
        Provides average and stddev of all OPs in self.traj
        """
        std=np.std(self.traj)
        # convert to numpy array
        return (np.mean(self.traj),std,std/np.sqrt(len(self.traj)-1))
    @property
    def get_op_res(self):
        """
	Provides average and stddev of all OPs in self.traj
        """

	# convert to numpy array
        return self.traj


#Anne: Added calc_angle from https://github.com/NMRLipids/MATCH/blob/master/scripts/calcOrderParameters.py
#Anne: removed self from function arguments
def calc_angle(atoms, com):
        """
        calculates the angle between the vector and z-axis in degrees
        no PBC check!
        Calculates the center of mass of the selected atoms to invert bottom leaflet vector
        """
        vec = atoms[1].position - atoms[0].position
        d = math.sqrt(np.square(vec).sum())
        cos = vec[2]/d
        # values for the bottom leaflet are inverted so that 
        # they have the same nomenclature as the top leaflet
        cos *= math.copysign(1.0, atoms[0].position[2]-com)
        try:
            angle = math.degrees(math.acos(cos))
        except ValueError:
            if abs(cos)>=1.0:
                print("Cosine is too large = {} --> truncating it to +/-1.0".format(cos))
                cos = math.copysign(1.0, cos)
                angle = math.degrees(math.acos(cos))
        return angle


def calc_z_dim(gro):
        u = mda.Universe(gro)
        z = u.dimensions[2]
        return z


def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.
    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use in trajectory
    for op in ordPars:
        # selection = pairs of atoms, split-by residues
        selection = mol.select_atoms("resname {rnm} and name {atA} {atB}".format(
                                    rnm=op.resname, atA=op.atAname, atB=op.atBname)
                                    ).atoms.split("residue")
    #    print(op.resname + " " + op.atAname + " " + op.atBname)
        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                #print(res.resnames, res.resids)
                for atom in res.atoms:
                  # print(atom.name, atom.id)
                   atA=op.atAname
                   atB=op.atBname
                   nat=res.n_atoms
                   warnings.warn("Selection >> name {atA} {atB} << \
                   contains {nat} atoms, but should contain exactly 2!")
        op.selection = selection

    # go through trajectory frame-by-frame
        Nres=len(op.selection)

    Nframes=len(mol.trajectory)
    for op in ordPars:
        op.traj= [0]*Nres
#        op.traj=[0]*Nres
#        if len(op.traj) == 0:
#            print(op.atAname + " " + op.atBname)
#        print("op.traj length " + str(len(op.traj)))

    for frame in mol.trajectory:
        for op in ordPars:            #.values():
            for i in range(0,Nres):
#             for i, op in enumerate(ordPars,1):
                residue=op.selection[i]
#                print(residue)
                S = op.calc_OP(residue)
#                print(S)
#                    print(op.atAname + " " + op.atBname)
#                    print(i)
#                op.traj.append(S/Nframes)
                op.traj[i]=op.traj[i]+S/Nframes
            #op.traj.append(np.mean(tmp))

        #print "--", mol.atoms[0].position
#    for op in ordPars:
#        op.traj=op.traj/Nframes

def parse_op_input(fname,lipid_name):
    ordPars = []
    atomC = []
    atomH = []
    resname = lipid_name

#    try:
#    with open(fname,"r") as f:
#        for ind, line in enumerate(f,1):
#            if line.startswith("#Whole molecules"):
#                print(getline(f.name, ind + 1).split())
#                resname = getline(f.name, ind + 1).split()[1]
#                break
    
    with open(fname, "r") as f:
        for line in f.readlines():
            if not line.startswith("#"):
                regexp1_H = re.compile(r'M_[A-Z0-9]*C[0-9]*H[0-9]*_M')
                regexp2_H = re.compile(r'M_G[0-9]*H[0-9]*_M')
                regexp3_H = re.compile(r'M_C[0-9]*H[0-9]*_M')
                regexp1_C = re.compile(r'M_[A-Z0-9]*C[0-9]*_M')
                regexp2_C = re.compile(r'M_G[0-9]_M')
                regexp3_C = re.compile(r'M_C[0-9]_M')

                if regexp1_C.search(line) or regexp2_C.search(line) or regexp3_C.search(line):
                    atomC = line.split()
                    atomH = []
                elif regexp1_H.search(line) or regexp2_H.search(line) or regexp3_H.search(line):
                    atomH = line.split()
                    if len(line.split())> 2:
                        resname = line.split()[2]
                else:
                    atomC = []
                    atomH = []

                if atomH:
                    items = [atomC[1], atomH[1], atomC[0], atomH[0]]
#                    print(resname + " " + atomH[1])
                    op = OrderParameter(resname, items[0], items[1], items[2], items[3])
                    ordPars.append(op)
#    except:
#        inpf=opts.inp_fname
#        raise RuntimeError("Couldn't read input file >> {inpf} <<")
    return ordPars


#def parse_op_input(fname):
#    """
#    parses input file with Order Parameter definitions
#    file format is as follows:
#    OP_name    resname    atom1    atom2
#    (flexible cols)
#    fname : string
#        input file name
#    returns : dictionary
#        with OrderParameters class instances
#    """
#    ordPars = OrderedDict()
#    try:
#        with open(fname,"r") as f:
#            for line in f.readlines():
#                if not line.startswith("#"):
#                    items = line.split()
#
#                    ordPars[items[0]] = OrderParameter(*items)
#
#    except:
#        inpf=opts.inp_fname
#        raise RuntimeError("Couldn't read input file >> {inpf} <<")
#    return ordPars



#%%

def find_OP(inp_fname, top_fname, traj_fname,lipid_name):
    ordPars = parse_op_input(inp_fname,lipid_name)

 #   for file_name in os.listdir(os.getcwd()):
 #       if file_name.startswith(traj_fname):
 #           trajs.append(file_name)

    read_trajs_calc_OPs(ordPars, top_fname, traj_fname)

    return ordPars

#


def read_trj_PN_angles(molname,atoms, top_fname, traj_fname, gro_fname):

    mol = mda.Universe(top_fname, traj_fname)
#    z_dim = calc_z_dim(gro_fname)

    selection = mol.select_atoms("resname " + molname + " and (name " + atoms[0] + ")", "resname " + molname + " and (name " + atoms[1] + ")").atoms.split("residue")
    com = mol.select_atoms("resname " + molname + " and (name " + atoms[0] + " or name " + atoms[1] + ")").center_of_mass()
    
    Nres=len(selection)
    Nframes=len(mol.trajectory)
    angles = np.zeros((Nres,Nframes))

#    sumsResAngles = [0]*Nres
    resAverageAngles = [0]*Nres
    resSTDerror = [0]*Nres
    j = 0

    for frame in mol.trajectory:
        for i in range(0,Nres):
            residue = selection[i]
            #print(residue)
            angles[i,j] = calc_angle(residue, com[2])
        j=j+1
#        sumsResAngles = map(add,sumsResAngles, frameAngles) 

#    resAverageAngles = [x / Nframes for x in sumsResAngles]
    for i in range(0,Nres):
        resAverageAngles[i] = sum(angles[i,:]) / Nframes
        resSTDerror[i] = np.std(angles[i,:])

    totalAverage = sum(resAverageAngles) / Nres
    totalSTDerror = np.std(resAverageAngles) / np.sqrt(Nres)

# standard errors
    



#    for i in range(0,Nres):
#        residue = selection[i]
#        PNangles = []
#        for frame in mol.trajectory:
#            PNangle = calc_angle(residue)
#            PNangles.append(PNangle)
#
#        averageAngle = sum(PNangles) / Nframes
#        resAveragePNangles.append(averageAngle)

    return angles, resAverageAngles, totalAverage, totalSTDerror

#####Dihedrals#####
#Anne: Functions needed for dihedral analysis

class Dihedrals:

    def __init__(self, resname, atom_1_name, atom_2_name, atom_3_name, atom_4_name, *args):
        self.resname = resname       # name of residue atoms are in
        self.at1name = atom_1_name
        self.at2name = atom_2_name
        self.at3name = atom_3_name
        self.at4name = atom_4_name
        self.name = atom_1_name + " " + atom_2_name + " " + atom_3_name + " " + atom_4_name  # generic name of atom A and atom B

        self.traj = []
        self.selection = []

#def calcDihedrals(atom1, atom2, atom3, atom4):
    
#    return

def parseDihedralInput(mapping_file):
#reads all heavy atoms from mapping file and forms groups of four consecutive atoms for dihedral calculation
#not ready yet!!!!
    G1 = []
    G2 = []
    G3 = []
    resname = ""

    regexp_H = re.compile(r'M_[A-Z0-9]*H[0-9]*_M')
    regexp_1 = re.compile(r'M_G[0-9]*[A-Z]*[0-9]*_M')
#    try:
    with open(mapping_file,"r") as f:
        for ind, line in enumerate(f,1):
            if line.startswith("#Whole molecules"):
#                print(getline(f.name, ind + 1).split())
                resname = getline(f.name, ind + 1).split()[1]
                break
    with open(mapping_file, "r") as f:
        for line in f.readlines():
            if not line.startswith("#") and not regexp_H.search(line):
                    if line.startswith("M_G1"):
                        atom = line.split()
                        G1.append(atom[0])
                    elif line.startswith("M_G2"):
                        atom = line.split()
                        G2.append(atom[0])
                    elif line.startswith("M_G3"):
                        atom = line.split()
                        G3.append(atom[0])

    G1.reverse()
    G2.reverse()

    G1G3 = G1 + G3
    G2G3 = G2 + G3

    dihAtoms = []

#make groups of 4 atoms from G1 tail
    dihGroups1 = []

    for atom in G1G3:
        if regexp_1.search(atom):
            dihAtoms.append(atom)

            if len(dihAtoms) == 4:
                if "M_G3" in dihAtoms[0]: #do not read groups of atoms containing only head group atoms yet
                    break
                dihGroups1.append(dihAtoms)
                #print(dihAtoms)
#                print(dihGroups)
                dihAtoms.pop(0) #remove first atom from dihAtoms
        else:
            dihAtoms.append(atom)
            if len(dihAtoms) == 4:
                #print(dihAtoms)
                dihGroups1.append(dihAtoms)
#                print(dihGroups)
                dihAtoms.pop(3) #remove last atom from dihAtoms
#empty dihAtoms list
    dihAtoms = []

#make group of 4 atoms from G2 tail and G3 head group
    dihGroups2 = []

    for atom in G2G3:
        if regexp_1.search(atom):
            dihAtoms.append(atom)

            if len(dihAtoms) == 4:
                dihGroups2.append(dihAtoms)
                #print(dihAtoms)
#                print(dihGroups)
                dihAtoms.pop(0) #remove first atom from dihAtoms
                
        else:
            dihAtoms.append(atom)
            if len(dihAtoms) == 4:
                dihGroups2.append(dihAtoms)
               # print(dihAtoms)
#                print(dihGroups)
                dihAtoms.pop(3) #remove last atom from dihAtoms


    return dihGroups1 + dihGroups2


#def readTrjCalcDih(trj, tpr, mapping_file):
