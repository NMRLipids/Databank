"""
    :meta private:
 calculation of order parameters of lipid bilayers
 from a MD trajectory

 meant for use with NMRlipids projects
 ------------------------------------------------------------
 Made by Joe,  Last edit 2017/02/02
 Refactored for performance by Gemini
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
import warnings
from tqdm import tqdm

# Maximum bond length for a C-H bond to be considered reasonable.
bond_len_max = 1.5  # in Angstrom
bond_len_max_sq = bond_len_max**2


class OrderParameter:
    """
    :meta private:

    Class for storing and manipulating order parameter (OP) related metadata
    (definition, name, etc.), OP trajectories, and methods to evaluate OPs.
    """

    def __init__(
        self,
        resname,
        atom_name_a,
        atom_name_b,
        univ_atom_name_a,
        univ_atom_name_b,
        *args,
    ):
        """
        Initializes the OrderParameter object. It doesn't matter which atom
        (A or B) comes first for the OP calculation.

        Args:
            resname (str): Name of the residue the atoms are in.
            atom_name_a (str): Name of the first atom in the topology.
            atom_name_b (str): Name of the second atom in the topology.
            univ_atom_name_a (str): Generic/mapping name for atom A.
            univ_atom_name_b (str): Generic/mapping name for atom B.
        """
        self.resname = resname
        self.aname_a = atom_name_a
        self.aname_b = atom_name_b
        self.m_aname_a = univ_atom_name_a
        self.m_aname_b = univ_atom_name_b
        self.name = f"{univ_atom_name_a} {univ_atom_name_b}"

        for field_name, field_value in self.__dict__.items():
            if isinstance(field_value, str):
                if not field_value.strip():
                    raise RuntimeError(
                        f"Provided name for field '{field_name}' is empty! "
                        "Cannot use empty names for atoms and OP definitions."
                    )
            else:
                warnings.warn(
                    f"Provided value for '{field_name}' is not a string: {field_value}. "
                    "Unexpected behaviour might occur."
                )

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
                f"Number of optional positional arguments is {len(args)}, not 0 or 2. "
                f"Args: {args}\nWrong file format?"
            )

        self.traj = []  # For storing final OP results.
        self.selection = []  # List of AtomGroups, one for each residue.
        self.atomgroup = None  # A single AtomGroup containing all atoms for this OP.

    @property
    def get_avg_std_stem_OP(self):  # noqa: N802 (API compliance)
        """Provides average, stddev, and standard error of the mean of OPs."""
        std = np.std(self.traj)
        n = len(self.traj)
        stem = std / np.sqrt(n - 1) if n > 1 else 0
        return np.mean(self.traj), std, stem


def read_trajs_calc_OPs(
    op_obj_list: list[OrderParameter], top: str, trajs: list[str]
):  # noqa: N802 (API)
    """
    :meta private:

    Creates an MDAnalysis Universe, reads trajectories, and calculates
    Order Parameters ("S") for each definition in `op_obj_list`. This version
    is optimized for single-core performance using vectorized calculations.
    """
    # --- 1. Setup Universe and Atom Selections ---
    mol = mda.Universe(top, trajs)
    improper_ops = []

    for i, op in enumerate(op_obj_list):
        sel_str = f"resname {op.resname} and name {op.aname_a} {op.aname_b}"
        selection_by_residue = mol.select_atoms(sel_str).split("residue")

        if not selection_by_residue:
            warnings.warn(
                f"Selection is empty: [{sel_str}]. "
                "Check residue and atom names in the mapping file.",
                UserWarning,
            )
            improper_ops.append(i)
            continue

        # Validate that each residue selection contains exactly two atoms
        valid_selection = []
        for res in selection_by_residue:
            if res.n_atoms != 2:
                warnings.warn(
                    f"Selection 'name {op.aname_a} {op.aname_b}' in residue "
                    f"{res.resids[0]} contains {res.n_atoms} atoms, but should be 2. "
                    "This residue will be skipped.",
                    UserWarning,
                )
            else:
                valid_selection.append(res)

        if not valid_selection:
            warnings.warn(
                f"No valid atom pairs found for selection: [{sel_str}]", UserWarning
            )
            improper_ops.append(i)
            continue

        op.selection = valid_selection
        # Create a single, combined AtomGroup for efficient, vectorized access
        op.atomgroup = mda.AtomGroup([atom for res in op.selection for atom in res])

    # Remove OP definitions that resulted in invalid selections
    for i in sorted(improper_ops, reverse=True):
        del op_obj_list[i]

    n_frames = len(mol.trajectory)
    for op in op_obj_list:
        n_res = len(op.selection)
        # We accumulate sums here, so initialize with zeros
        op.traj = np.zeros(n_res, dtype=np.float64)

    print("Processing trajectory with optimized single-core engine...")
    for _ in tqdm(mol.trajectory, total=n_frames, unit="frame"):
        for op in op_obj_list:
            if op.atomgroup is None or len(op.atomgroup) == 0:
                continue

            # Get all atom positions for this OP in one go
            # Shape: (n_residues * 2, 3)
            positions = op.atomgroup.positions

            # Reshape to easily access atom pairs
            # Shape: (n_residues, 2, 3) where axis 1 is [atom_A, atom_B]
            positions_reshaped = positions.reshape(-1, 2, 3)

            # Calculate vectors between atom pairs for all residues at once
            vec = positions_reshaped[:, 1, :] - positions_reshaped[:, 0, :]

            # Calculate squared distance for all residues
            d2 = np.sum(vec**2, axis=1)

            # Create a mask for valid bond lengths to avoid unnecessary calculations
            # and warnings for atoms that are too far apart (e.g., due to PBC issues).
            valid_mask = d2 <= bond_len_max_sq

            # Initialize cos2 array. We only compute for valid pairs.
            cos2 = np.zeros_like(d2)

            # Safely calculate cosine-squared of the angle with the z-axis
            # for all valid vectors simultaneously.
            # np.divide handles potential division by zero if d2 is 0.
            d2_valid = d2[valid_mask]
            vec_valid = vec[valid_mask]
            cos2[valid_mask] = np.divide(
                vec_valid[:, 2] ** 2, d2_valid, where=d2_valid != 0
            )

            # Calculate order parameters for all residues
            op_values = 0.5 * (3.0 * cos2 - 1.0)

            # Add the results for the current frame to the running sum.
            # We only add the valid ones, others remain 0 for this frame.
            op.traj += op_values

    # Average the accumulated sums over all frames
    for op in op_obj_list:
        if n_frames > 0:
            op.traj /= n_frames
        # Convert back to a list to maintain original API behavior
        op.traj = op.traj.tolist()


def parse_op_input(mapping_dict: dict, lipid_resname: str):
    """
    Parses a mapping dictionary to form a list of C-H pairs for OP calculation.

    Args:
        mapping_dict (dict): The mapping dictionary.
        lipid_resname (str): The default lipid residue name.

    Returns:
        list[OrderParameter]: A list of OrderParameter instances.
    """
    opvals = []
    atom_c = []
    atom_h = []
    resname = lipid_resname

    # Regex to identify carbon and hydrogen atoms from mapping keys
    regexp_h = re.compile(r"M_([A-Z0-9]*C[0-9]*|G[0-9]*|C[0-9]*)H[0-9]*_M")
    regexp_c = re.compile(r"M_([A-Z0-9]*C[0-9]*|G[0-9]{1,2}|C[0-9]{1,2})_M")

    for mapping_key, value in mapping_dict.items():
        if not isinstance(value, dict) or "ATOMNAME" not in value:
            continue

        if regexp_c.search(mapping_key) and not regexp_h.search(mapping_key):
            atom_c = [mapping_key, value["ATOMNAME"]]
            # Use residue name from mapping if available, otherwise use default
            resname = value.get("RESIDUE", lipid_resname)
            atom_h = []  # Reset hydrogen atom
        elif regexp_h.search(mapping_key):
            atom_h = [mapping_key, value["ATOMNAME"]]
        else:
            atom_c, atom_h = [], []  # Reset for non-matching keys

        if atom_h and not atom_c:
            warnings.warn(
                f"Cannot define carbon for the hydrogen {atom_h[0]} ({atom_h[1]})",
                UserWarning,
            )
            continue

        # If both a carbon and a hydrogen have been found, create the pair
        if atom_h and atom_c:
            op = OrderParameter(resname, atom_c[1], atom_h[1], atom_c[0], atom_h[0])
            opvals.append(op)
            # Important: Reset hydrogen to look for the next one for the same carbon
            atom_h = []

    return opvals


def find_OP(
    mdict: dict, top_fname: str, traj_fname: str, lipid_name: str
):  # noqa: N802 (API)
    """
    Externally used function for computing OP values.

    Args:
        mdict (dict): The mapping dictionary.
        top_fname (str): Filename of the topology file (e.g., .gro, .tpr).
        traj_fname (str or list[str]): Filename(s) of the trajectory file(s).
        lipid_name (str): The residue name of the lipid.

    Returns:
        list[OrderParameter]: A list of OrderParameter instances with calculated data.
    """
    op_pairs = parse_op_input(mdict, lipid_name)
    if not isinstance(traj_fname, list):
        traj_fname = [traj_fname]
    read_trajs_calc_OPs(op_pairs, top_fname, traj_fname)
    return op_pairs
