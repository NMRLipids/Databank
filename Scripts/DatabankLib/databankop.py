"""
    :meta private:
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

import sys
import re
import MDAnalysis as mda
import numpy as np
import warnings  # TODO: should we change to NMRlipids' logger?
from tqdm import tqdm

from DatabankLib.databankLibrary import loadMappingFile

bond_len_max = 1.5  # in A, max distance between atoms for reasonable OP calculation
bond_len_max_sq = bond_len_max**2


class OrderParameter:
    """
    :meta private:

    Class for storing&manipulating
    order parameter (OP) related metadata (definition, name, ...)
    and OP trajectories
    and methods to evaluate OPs.
    """

    def __init__(
        self, resname, atom_A_name, atom_B_name, M_atom_A_name, M_atom_B_name, *args
    ):
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        self.resname = resname  # name of residue atoms are in
        self.atAname = atom_A_name
        self.atBname = atom_B_name
        self.M_atAname = M_atom_A_name
        self.M_atBname = M_atom_B_name
        self.name = (
            M_atom_A_name + " " + M_atom_B_name
        )  # generic name of atom A and atom B
        for field in self.__dict__:
            if not isinstance(field, str):
                warnings.warn(
                    "provided name >> {} << is not a string! \n \
                Unexpected behaviour might occur."
                )  # .format(field)
            else:
                if not field.strip():
                    raise RuntimeError(
                        "provided name >> {} << is empty! \n \
                    Cannot use empty names for atoms and OP definitions."
                    )  # .format(field)
        # extra optional arguments allow setting avg,std values -- suitable for
        # reading-in results of this script
        if len(args) == 0:
            self.avg = None
            self.std = None
            self.stem = None
        elif len(args) == 2:
            self.avg = args[0]
            self.std = args[1]
            self.stem = None
        else:
            warnings.warn(
                f"Number of optional positional arguments is {len}, not 2 or 0."
                f" Args: {args}\nWrong file format?"
            )
        self.traj = []  # for storing OPs
        self.selection = []

    @staticmethod
    def calc_OP(atoms):
        """
        calculates Order Parameter according to equation
        S = 1/2 * (3*cos(theta)^2 -1)
        """
        vec = atoms[1].position - atoms[0].position
        d2 = np.square(vec).sum()

        if d2 > bond_len_max_sq:
            at1 = atoms[0].name
            at2 = atoms[1].name
            resnr = atoms[0].resid
            d = np.sqrt(d2)
            warnings.warn(
                f"Atomic distance for atoms"
                f"{at1} and {at2} in residue no. {resnr} is suspiciously "
                f"long: {d}!\nPBC removed???"
            )
        cos2 = vec[2] ** 2 / d2
        S = 0.5 * (3.0 * cos2 - 1.0)
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
        std = np.std(self.traj)
        # convert to numpy array
        return (np.mean(self.traj), std, std / np.sqrt(len(self.traj) - 1))

    @property
    def get_op_res(self):
        """
        Provides average and stddev of all OPs in self.traj
        """

        # convert to numpy array
        return self.traj


def read_trajs_calc_OPs(ordPars, top, trajs):
    """
    :meta private:

    procedure that
    creates MDAnalysis Universe with top,
    reads in trajectories trajs and then
    goes through every frame and
    evaluates each Order Parameter "S" from the list of OPs ordPars.
    ordPars : list of OrderParameter class
       each item in this list describes an Order parameter to be calculated in the
       trajectory
    top : str
        filename of a top file (e.g. conf.gro)
    trajs : list of strings
        filenames of trajectories
    """
    # read-in topology and trajectory
    mol = mda.Universe(top, trajs)

    # make atom selections for each OP and store it as its attribute for later use
    # in trajectory
    c = -1
    improperOPs = []
    for op in ordPars:
        c += 1
        # selection = pairs of atoms, split-by residues
        selStr = "resname {rnm} and name {atA} {atB}".format(
                rnm=op.resname, atA=op.atAname, atB=op.atBname)
        selection = mol.select_atoms(selStr).atoms.split("residue")
        if len(selection) == 0:
            warnings.warn(
                f"Selection is empty: [{selStr}]. "
                f"Check carefully residue name and names in the mapping file.")
            improperOPs.append(c)
            continue

        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                atA = op.atAname
                atB = op.atBname
                nat = res.n_atoms
                print(atA, atB, nat)
                warnings.warn(
                    f"Selection >> name {atA} {atB} << "
                    f"contains {nat} atoms, but should contain exactly 2!")
                improperOPs.append(c)
                continue

        op.selection = selection

    # remove OPs, which are incorrect
    improperOPs.sort(reverse=True)
    for i in improperOPs:
        del ordPars[i]

    # go through trajectory frame-by-frame
    Nframes = len(mol.trajectory)
    for op in ordPars:
        # print(op.selection)
        Nres = len(op.selection)
        op.traj = [0] * Nres

    for frame in tqdm(mol.trajectory):
        for op in ordPars:  # .values():
            Nres = len(op.selection)
            for i in range(0, Nres):
                residue = op.selection[i]
                S = OrderParameter.calc_OP(residue)
                op.traj[i] = op.traj[i] + S / Nframes


def parse_op_input(mapping_file: str, lipid_name: str):
    """Form pair-list of all OPs which should be calculated

    Args:
        mapping_file (str): mapping file name
        lipid_name (str): lipid name (residue name)

    Returns:
        array: List of OrderParameter instances
    """
    ordPars = []
    atomC = []
    atomH = []
    resname = lipid_name
    mapping_dict = loadMappingFile(mapping_file)

    regexp1_H = re.compile(r"M_[A-Z0-9]*C[0-9]*H[0-9]*_M")
    regexp2_H = re.compile(r"M_G[0-9]*H[0-9]*_M")
    regexp3_H = re.compile(r"M_C[0-9]*H[0-9]*_M")
    regexp1_C = re.compile(r"M_[A-Z0-9]*C[0-9]*_M")
    regexp2_C = re.compile(r"M_G[0-9]{1,2}_M")
    regexp3_C = re.compile(r"M_C[0-9]{1,2}_M")

    for mapping_key in mapping_dict.keys():
        if (
            regexp1_C.search(mapping_key)
            or regexp2_C.search(mapping_key)
            or regexp3_C.search(mapping_key)
        ):
            atomC = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
            try:
                resname = mapping_dict[mapping_key]["RESIDUE"]
            except (KeyError, TypeError):
                pass
            atomH = []
        elif (
            regexp1_H.search(mapping_key)
            or regexp2_H.search(mapping_key)
            or regexp3_H.search(mapping_key)
        ):
            atomH = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
        else:
            atomC = []
            atomH = []

        if atomH and not len(atomC):
            print(
                f"Cannot define carbon for the hydrogen {atomH[0]} ({atomH[1]})",
                file=sys.stderr)
            continue
        if atomH and len(atomC):
            items = [atomC[1], atomH[1], atomC[0], atomH[0]]
            op = OrderParameter(resname, items[0], items[1], items[2], items[3])
            ordPars.append(op)
    return ordPars


def find_OP(inp_fname: str, top_fname: str, traj_fname: str, lipid_name: str):
    """Externally used funcion for computing OP values.

    Args:
        inp_fname (_type_): mapping file
        top_fname (_type_): TPR file
        traj_fname (_type_): TRAJ file
        lipid_name (_type_): lipid name (residue name)

    Returns:
        ordPars: list of OrderParameter instances
    """
    ordPars = parse_op_input(inp_fname, lipid_name)
    read_trajs_calc_OPs(ordPars, top_fname, traj_fname)
    return ordPars
