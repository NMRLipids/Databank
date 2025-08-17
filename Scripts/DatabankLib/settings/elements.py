"""
Elements assigneer functions. Helps Molecule to define brutto formula
from existing mapping files and element guesser.
"""

import json
import os
from DatabankLib import NMLDB_DATA_PATH
from DatabankLib.core import System
import MDAnalysis as mda
import periodictable
import re


def uname2element(mapping_name: str):
    """Guess element name from universal atom name.

    :param mapping_name: name of the mapping, e.g. 'M_C1', 'M_G'
    :raises KeyError: if element cannot be determined from mapping name
    :return: element name, e.g. 'C', 'N', ...
    """
    # removes numbers and '_M' from the end of string
    name1 = re.sub(r"[0-9]*_M$", "", mapping_name)
    # removes M_X12Y23... sequence from the beginning
    name2 = re.sub(r"^M_([A-Z]{1,2}[0-9]{1,4})*", "", name1)

    # name2 is the atom and electrons are assigned to this atom

    if name2 == "G":  # G is a glycerol carbon so change G to C
        name2 = "C"

    try:
        if name2 != "X":
            # X is a placeholder for unknown element, we use 'X' as element name
            _ = getattr(periodictable, name2).number
    except AttributeError as e:
        raise KeyError(
            "This mapping name cannot be read by our rules: {mapping_name}") from e
    else:
        return name2


def guess_elements(system: System, u: mda.Universe) -> None:
    """Assigns elements to atoms in the MDAnalysis universe based on the system's
    composition.

    :param system: Databank System object
    :type system: System
    :param u: MDAnalysis Universe object
    :type u: mda.Universe
    """
    for _mol, comp in system['COMPOSITION'].items():
        # get the name of molecule used in simulation files
        mol_simu_name = comp['NAME']
        mol = system.content[_mol]

        # if lipid is split to multiple residues
        _moltype = f"moltype {comp['MOLTYPE']} and " if 'MOLTYPE' in comp else ""

        try:
            ua_dict_f = os.path.join(
                NMLDB_DATA_PATH, "lipid_json_buildH",
                system["UNITEDATOM_DICT"][_mol]+".json")
            with open(ua_dict_f, 'r') as f:
                ua_dict = json.load(f)
        except (KeyError, TypeError):
            ua_dict = False

        for uname, props in mol.mapping_dict.items():
            try:
                res = str(props['RESIDUE'])
            except KeyError:
                res = ""
            _residue = f"resname {mol_simu_name if res == '' else res} and "
            _name = f"name {props['ATOMNAME']}"
            selstr = _moltype + _residue + _name
            selection = u.select_atoms(selstr)
            if ua_dict:
                if props['ATOMNAME'] in ua_dict:
                    elname = ua_dict[props['ATOMNAME']][0]
                    if elname == "CH":
                        elname = "CH1"
                    if elname == "CHdoublebond":
                        elname = "CH1"
                else:
                    if selection.n_atoms == 0:
                        continue  # this is an implicit atom
                    else:
                        elname = uname2element(uname)
            else:
                # if not united-atom, ALL atoms must be present
                assert selection.n_atoms > 0, \
                    f"Selection '{selstr}' did not match any atoms in the universe."
                elname = uname2element(uname)
            selection.atoms.elements = selection.n_atoms*[elname]
        # end mapping loop
    # end molecules loop
    return
