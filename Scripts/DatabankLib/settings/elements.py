"""
Elements assigneer functions. Helps Molecule to define brutto formula
from existing mapping files and element guesser.
"""

from DatabankLib.core import System
import MDAnalysis as mda
import periodictable
import re


def uname2element(mapping_name):
    # removes numbers and '_M' from the end of string
    name1 = re.sub(r"[0-9]*_M$", "", mapping_name)
    # removes M_X12Y23... sequence from the beginning
    name2 = re.sub(r"^M_([A-Z]{1,2}[0-9]{1,4})*", "", name1)

    # name2 is the atom and electrons are assigned to this atom

    if name2 == "G":  # G is a glycerol carbon so change G to C
        name2 = "C"

    try:
        el = getattr(periodictable, name2).number
    except AttributeError as e:
        raise KeyError("This mapping name cannot be read by our rules: {mapping_name}") from e

    return name2

def guess_elements(system: System, u: mda.Universe) -> None:
    for _mol,comp in system['COMPOSITION'].items():
        # get the name of molecule used in simulation files
        mol_simu_name = comp['NAME']
        mol = system.content[_mol]

        # if lipid is split to multiple residues
        _moltype = f"moltype {comp['MOLTYPE']} and " if 'MOLTYPE' in comp else ""

        for uname,props in mol.mapping_dict.items():
            try: 
                res = str(props['RESIDUE'])
            except KeyError:
                res = ""
            _residue = f"resname {mol_simu_name if res == '' else res} and "
            _name = f"name {props['ATOMNAME']}"
            selstr = _moltype + _residue + _name
            selection = u.select_atoms(selstr)
            assert (selection.n_atoms > 0), \
                f"No atoms found for selection: {selstr}"
            elname = uname2element(uname)
            selection.atoms.elements = selection.n_atoms*[elname]
        # end mapping loop
    # end molecules loop
    return


