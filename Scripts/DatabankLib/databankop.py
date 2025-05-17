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
        self, resname, atom_name_a, atom_name_b,
        univ_atom_name_a, univ_atom_name_b, *args
    ):
        """
        it doesn't matter which comes first,
        atom A or B, for OP calculation.
        """
        self.resname = resname  # name of residue atoms are in
        self.aname_a = atom_name_a
        self.aname_b = atom_name_b
        # TODO: consider removing these two variables
        self.m_aname_a = univ_atom_name_a
        self.m_aname_b = univ_atom_name_b
        self.name = (
            univ_atom_name_a + " " + univ_atom_name_b
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
    def calc_OP(atoms):  # noqa: N802 (API)
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
        val = 0.5 * (3.0 * cos2 - 1.0)
        return val

    @property
    def get_avg_std_OP(self):  # noqa: N802 (API)
        """
        Provides average and stddev of all OPs in self.traj
        """
        # convert to numpy array
        return (np.mean(self.traj), np.std(self.traj))

    @property
    def get_avg_std_stem_OP(self):  # noqa: N802 (API)
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


def read_trajs_calc_OPs(  # noqa: N802 (API)
        op_obj_list: list[OrderParameter], top, trajs):
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
    improper_ops = []
    for op in op_obj_list:
        c += 1
        # selection = pairs of atoms, split-by residues
        sel_str = "resname {rnm} and name {atA} {atB}".format(
                rnm=op.resname, atA=op.aname_a, atB=op.aname_b)
        selection = mol.select_atoms(sel_str).atoms.split("residue")
        if len(selection) == 0:
            warnings.warn(
                f"Selection is empty: [{sel_str}]. "
                f"Check carefully residue name and names in the mapping file.")
            improper_ops.append(c)
            continue

        for res in selection:
            # check if we have only 2 atoms (A & B) selected
            if res.n_atoms != 2:
                at_a = op.atAname
                at_b = op.atBname
                nat = res.n_atoms
                print(at_a, at_b, nat)
                warnings.warn(
                    f"Selection >> name {at_a} {at_b} << "
                    f"contains {nat} atoms, but should contain exactly 2!")
                improper_ops.append(c)
                continue

        op.selection = selection

    # remove OPs, which are incorrect
    improper_ops.sort(reverse=True)
    for i in improper_ops:
        del op_obj_list[i]

    # go through trajectory frame-by-frame
    n_frames = len(mol.trajectory)
    for op in op_obj_list:
        # print(op.selection)
        n_res = len(op.selection)
        op.traj = [0] * n_res

    for frame in tqdm(mol.trajectory):
        for op in op_obj_list:  # .values():
            n_res = len(op.selection)
            for i in range(0, n_res):
                residue = op.selection[i]
                opval = OrderParameter.calc_OP(residue)
                op.traj[i] = op.traj[i] + opval / n_frames


def parse_op_input(mapping_dict: dict, lipid_resname: str):
    """Form pair-list of all OPs which should be calculated

    Args:
        mapping_dict (dict): mapping dictionary
        lipid_resname (str): lipid name (residue name)

    Returns:
        array: List of OrderParameter instances
    """
    opvals = []
    atom_c = []
    atom_h = []
    resname = lipid_resname

    regexp1_h = re.compile(r"M_[A-Z0-9]*C[0-9]*H[0-9]*_M")
    regexp2_h = re.compile(r"M_G[0-9]*H[0-9]*_M")
    regexp3_h = re.compile(r"M_C[0-9]*H[0-9]*_M")
    regexp1_c = re.compile(r"M_[A-Z0-9]*C[0-9]*_M")
    regexp2_c = re.compile(r"M_G[0-9]{1,2}_M")
    regexp3_c = re.compile(r"M_C[0-9]{1,2}_M")

    for mapping_key in mapping_dict:
        if (
            regexp1_c.search(mapping_key)
            or regexp2_c.search(mapping_key)
            or regexp3_c.search(mapping_key)
        ):
            atom_c = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
            try:
                resname = mapping_dict[mapping_key]["RESIDUE"]
            except (KeyError, TypeError):
                pass
            atom_h = []
        elif (
            regexp1_h.search(mapping_key)
            or regexp2_h.search(mapping_key)
            or regexp3_h.search(mapping_key)
        ):
            atom_h = [mapping_key, mapping_dict[mapping_key]["ATOMNAME"]]
        else:
            atom_c = []
            atom_h = []

        if atom_h and not len(atom_c):
            print(
                f"Cannot define carbon for the hydrogen {atom_h[0]} ({atom_h[1]})",
                file=sys.stderr)
            continue
        if atom_h and len(atom_c):
            items = [atom_c[1], atom_h[1], atom_c[0], atom_h[0]]
            op = OrderParameter(resname, items[0], items[1], items[2], items[3])
            opvals.append(op)
    return opvals


def find_OP(  # noqa: N802 (API)
        mdict: dict, top_fname: str, traj_fname: str, lipid_name: str):
    """Externally used funcion for computing OP values.

    Args:
        mdict (_type_): mapping dict
        top_fname (_type_): TPR file
        traj_fname (_type_): TRAJ file
        lipid_name (_type_): lipid name (residue name)

    Returns:
        ordPars: list of OrderParameter instances
    """
    op_pairs = parse_op_input(mdict, lipid_name)
    read_trajs_calc_OPs(op_pairs, top_fname, traj_fname)
    return op_pairs
