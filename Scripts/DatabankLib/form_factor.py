"""
Module form_factor serves for bilayer form factor calculations.
"""

import sys
import os
import re
import MDAnalysis as mda
import numpy as np
import json
import time
from pprint import pprint
from tqdm import tqdm

# for testing purposes to safe time for loading trajectory and creating dictonary
# the electron.dat will be decripted in the future as the values will be used from
# the universal mapping file
import periodictable

from DatabankLib.databankLibrary import lipids_dict, loadMappingFile, getLipids
from DatabankLib.jsonEncoders import CompactJSONEncoder
from DatabankLib import NMLDB_ROOT_PATH


# To write data in numpy arrays into json file, we inherit compact JSON
# encoder to make him store 2xN numpy arrays as just a nested list
class NumpyArrayEncoder(CompactJSONEncoder):
    def encode(self, o):
        if isinstance(o, np.ndarray):
            return CompactJSONEncoder.encode(self, o.tolist())
        else:
            return CompactJSONEncoder.encode(self, o)


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

    # Group of static private helper methods

    @staticmethod
    def __listNamePairs(mapping_dict, molecule):
        pairs_dictionary = {}
        residues = [mv['RESIDUE']
                    for _, mv in mapping_dict.items()
                    if 'RESIDUE' in mv]

        if residues:
            pairs_dictionary = dict.fromkeys(residues)
            print(pairs_dictionary.keys())
            for res in pairs_dictionary.keys():
                pairs_dictionary[res] = [
                    [mk, mv['ATOMNAME']]
                    for mk, mv in mapping_dict.items()
                    if res == mv['RESIDUE']
                    ]
        else:
            pairs_dictionary = {
                molecule: [[mn, mv['ATOMNAME']] for mn, mv in mapping_dict.items()]
            }

        return pairs_dictionary

    @staticmethod
    def __filterHbonds(mapping_names):
        re_H = re.compile(r'M_([A-Z]{1,2}[0-9]{1,4})*H[0-4]{1,2}_M')
        filtered = filter(re_H.match, mapping_names)
        return filtered

    @staticmethod
    def __getWater(rmf):
        return ('resname ' + rmf['COMPOSITION']['SOL']['NAME'])

    # -- regular methods --

    def __init__(self, path, conf, traj, nbin, output, readme, density_type="electron"):
        self.path = path
        self.conf = conf
        self.traj = traj
        self.readme = readme
        start_time = time.time()
        try:
            self.u = mda.Universe(self.conf, self.traj)
            print()
        except Exception as e:
            if (self.conf[-3:] == 'gro'):
                raise e
            # assume conf was TPR
            gro = self.path + os.sep + 'conf.gro'
            print("Generating conf.gro because MDAnalysis cannot read tpr version")
            os.system(f"echo System | gmx trjconv "
                      f"-s {self.conf} -f {self.traj} -dump 0 -o {gro}")
            self.conf = gro
            self.u = mda.Universe(self.conf, self.traj)

        print("Loading the trajectory takes {:10.6f} s".format(time.time()-start_time))

        # number of bins
        self.nbin = nbin
        # the totatl box size in [nm] - will be probably removed and tpr box size or
        # transform of final FF will be used instead

        self.output = path + output
        self.lipids = getLipids(readme)
        self.waters = FormFactor.__getWater(readme)

        self.density_type = density_type

        self.system_mapping = self.temporary_mapping_dictionary()

        self.calculate_weight()
        self.calculate_density()

    def temporary_mapping_dictionary(self):
        mapping_dictionary = {}
        for key1 in self.readme["COMPOSITION"].keys():  #

            key2 = self.readme["COMPOSITION"][key1]["NAME"]
            mapping_file = self.readme["COMPOSITION"][key1]["MAPPING"]
            mapping_dict = loadMappingFile(mapping_file)

            # put mapping files into a larger dictionary where molecule name is key and
            # the dictionary in the correponding mapping file is the value
            mapping_dictionary[key2] = mapping_dict  # try if this works

        return mapping_dictionary

    def getElectrons(self, mapping_name):
        # removes numbers and '_M' from the end of string
        name1 = re.sub(r"[0-9]*_M$", "", mapping_name)
        # removes M_X12Y23... sequence from the beginning
        name2 = re.sub(r"^M_([A-Z]{1,2}[0-9]{1,4})*", "", name1)

        # name2 is the atom and electrons are assigned to this atom

        if name2 == "G":  # G is a glycerol carbon so change G to C
            name2 = "C"

        try:
            el = getattr(periodictable, name2).number
        except AttributeError:
            print(
                f"ERROR: This mapping name cannot be read by our rules: {mapping_name}",
                file=sys.stderr,
            )
            print("Consider changing naming in your mapping file.")
            sys.exit(1)

        # TODO Modify default charge modification
        if name2 in ['K', 'Na', 'Cs']:
            el -= 1
        elif name2 in ['Cl']:
            el += 1
        #

        # TODO REMOVE WEIRD HERITAGE!!!
        if name2 == "D":
            el = 0
        # WEIRD HERITAGE!!!

        return el

    def residueElectronsAll(self, molecule):
        """
        Return list of electrons
        """
        print(f"Gathering electron list for {molecule}")
        electrons = []
        # get the name of molecule used in simulation files
        molname = self.readme['COMPOSITION'][molecule]['NAME']

        pairs_residue = FormFactor.__listNamePairs(
            self.system_mapping[molname], molname)

        # if lipid is split to multiple residues
        selection_txt = ""

        try:
            selection_txt = (
                "moltype " +
                self.readme['COMPOSITION'][molecule]['MOLTYPE'] +
                " and ")
        except KeyError:
            pass

        for res in pairs_residue.keys():  # fix second res in explicit_atoms
            # TODO: REPLACE THIS ENTIRE CODE WITH MDANALYSIS
            # when lipids are split into more than one residue in topology with
            # same residue names e.g. PGR,PA,OL and PE, PA,OL then
            # `u.select_atoms("resname " + res)` does not distinguish
            # which PA or PL belongs to POPG or POPE
            selection_txt = selection_txt + "resname " + res
            # explicit atoms of the residue
            explicit_atoms = list(self.u.select_atoms(selection_txt).atoms.names)
            for atom in explicit_atoms:
                # find generic mapping name matching to forcefield atom name
                atom_i1 = None
                for i in range(len(pairs_residue[res])):
                    _cm = pairs_residue[res][i]  # [UNIQUE_NAME, NAME] pair
                    if re.match(_cm[1].replace('+', '\\+'), atom):
                        atom_i1 = i
                        break
                if atom_i1 is None:
                    raise ValueError(f"Atom was not found: {res}:{atom}")
                # # find mapping name
                # index_list = [atom in pairs for pairs in pairs_residue[res]]
                # print (atom, pairs_residue[res])
                # atom_i1 = index_list.index(True)

                mapping_name = pairs_residue[res][atom_i1][0]
                # get number of electrons in an atom i of residue
                e_atom_i = self.getElectrons(mapping_name)
                electrons.append(e_atom_i)

        return electrons

    def assignElectrons(self):
        # check if simulation is united atom or all atom
        try:
            UA = self.readme['UNITEDATOM_DICT']
        except KeyError:
            UA = False
        else:
            UA = True

        weights = []
        l_weights = []  # separate list for lipid electrons
        w_weights = []  # separate list for water electrons

        if UA:  # UNITED ATOM SIMULATIONS
            for molecule1 in self.readme['COMPOSITION'].keys():
                print(molecule1)

                if molecule1 in lipids_dict:
                    molecule2 = self.readme['COMPOSITION'][molecule1]['NAME']
                    electrons = []

                    pairs_residue = FormFactor.__listNamePairs(
                        self.system_mapping[molecule2], molecule2)
                    atomsH = FormFactor.__filterHbonds(
                        self.system_mapping[molecule2].keys())  # list of hydrogen atoms

                    # extract explicit atoms and get the mapping names
                    for res in pairs_residue.keys():
                        # explicit atoms of the residue
                        explicit_atoms = list(self.u.select_atoms(
                            "resname " + res).atoms.names)
                        for atom in explicit_atoms:
                            # print(atom)
                            # find generic mapping name matching to forcefield atom name
                            index_list = [atom in pairs for pairs in pairs_residue[res]]
                            atom_i1 = index_list.index(True)
                            mapping_name = pairs_residue[res][atom_i1][0]
                            # remove _M from the end of mapping name
                            name1 = re.sub(r'_M', '', mapping_name[::1])

                            # count how many times name1 occurs in atomsH i.e. how
                            # many H are bound to it
                            numberH = sum(
                                1 for atomH in atomsH
                                if name1 == re.sub(r'H[1-4]_M', '', atomH[::1])
                                )
                            number_e = self.getElectrons(mapping_name)
                            e_atom_i = number_e + numberH
                            electrons.append(e_atom_i)

                    weights.extend(electrons)
                    l_weights.extend(electrons)
                else:
                    w = self.residueElectronsAll(molecule1)
                    weights.extend(w)
                    if molecule1 == 'SOL':
                        w_weights.extend(w)

        else:  # ALL ATOM SIMULATIONS
            for molecule in self.readme['COMPOSITION'].keys():
                w = self.residueElectronsAll(molecule)
                weights.extend(w)
                if molecule in lipids_dict:
                    l_weights.extend(w)
                elif molecule == 'SOL':
                    w_weights.extend(w)
        # how to do lipids split into three residues??    MAYBE WORKS

        return weights, l_weights, w_weights

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
        start_time = time.time()
        if self.density_type == "electron":
            # calc weights to calculate the electron density
            # with the function that maps the electrons to atoms
            self.wght, self.lipid_wght, self.water_wght = self.assignElectrons()
            print("lenght of wght: " + str(len(self.wght)))
            print("lenght of lipid_wght: " + str(len(self.lipid_wght)))
            print("atoms in system: " + str(len(self.u.atoms.names)))

            if len(self.wght) != len(self.u.atoms.names):
                print("ERROR: Number of atoms mismatch")
                print("ERROR: If lipids are split to many residues add moltype "
                      "of lipids from tpr as 'MOLTYPE' to README.yaml "
                      "under 'COMPOSITION'.")

        if self.density_type == "number":
            self.wght = np.ones(self.u.atoms.names.shape[0])
            clipids = self.u.select_atoms(self.lipids)
            self.lipid_wght = np.ones(clipids.atoms.names.shape[0])
            cwaters = self.u.select_atoms(self.waters)
            self.lipid_wght = np.ones(cwaters.atoms.names.shape[0])

        if self.density_type == "mass":
            self.wght = self.u.atoms.masses
            clipids = self.u.select_atoms(self.lipids)
            self.lipid_wght = np.ones(clipids.atoms.names.shape[0])
            cwaters = self.u.select_atoms(self.waters)
            self.lipid_wght = np.ones(cwaters.atoms.names.shape[0])

        print("Creating the electron mapping dictonary takes {:10.6f} s"
              .format(time.time()-start_time))

    def calculate_density(self):
        c = self.u.select_atoms(self.lipids)
        print(c)  # TODO: remove excessive debug printing!

        box_z = self.u.dimensions[2]  # + 10 if fails
        print(box_z)
        d = box_z/10 / self.nbin  # bin width
        boxH = box_z/10
        print(boxH)
        x = np.linspace(-boxH/2, boxH/2, self.nbin+1)[:-1] + d/2
        density_z_centered = np.zeros(self.nbin)
        density_z_no_center = np.zeros(self.nbin)
        density_lipids_center = np.zeros(self.nbin)
        density_waters_center = np.zeros(self.nbin)

        # Calculte density profiles and FF from individual frames
        start_time = time.time()
        min_z = 10000000
        frame = 0

        ElectronNumbers = {}

        try:
            UA = self.readme['UNITEDATOM_DICT']
        except KeyError:
            UA = False
        else:
            if self.readme['UNITEDATOM_DICT'] is None:
                UA = False
            else:
                UA = True

        for key1 in self.readme['COMPOSITION'].keys():
            mol = self.readme['COMPOSITION'][key1]['NAME']
            print(mol)
            for universalAN in self.system_mapping[mol]:
                AtomName = self.system_mapping[mol][universalAN]['ATOMNAME']
                try:
                    ResName = self.system_mapping[mol][universalAN]['RESIDUE']
                except (KeyError, TypeError):
                    ResName = mol

                if ResName not in ElectronNumbers.keys():
                    ElectronNumbers[ResName] = {}

                if UA and key1 in lipids_dict:
                    UAlipidjsonNAME = os.path.join(
                        NMLDB_ROOT_PATH, 'Scripts', 'DatabankLib',
                        'lipid_json_buildH',
                        self.readme['UNITEDATOM_DICT'][key1] + '.json')

                    with open(UAlipidjsonNAME) as json_file:
                        UAlipidjson = json.load(json_file)

                    numberH = 0
                    try:
                        if UAlipidjson[AtomName][0] == "CH":
                            numberH = 1
                        if UAlipidjson[AtomName][0] == "CH2":
                            numberH = 2
                        if UAlipidjson[AtomName][0] == "CH3":
                            numberH = 3
                    except (KeyError, TypeError):
                        pass

                    ElectronNumbers[ResName][AtomName] = \
                        self.getElectrons(universalAN) + numberH

                else:
                    ElectronNumbers[ResName][AtomName] = self.getElectrons(universalAN)

        print("ElectronNumbers dictionary: ")
        pprint(ElectronNumbers)

        # define selections and dictionaries for profile calculations
        clipids = self.u.select_atoms(self.lipids)
        cwaters = self.u.select_atoms(self.waters)

        print("Generating final electron arrays from frame #1..", end='')
        weightsALL = np.zeros(self.u.atoms.n_atoms)

        for rname, AName2Enum in ElectronNumbers.items():
            for aname, enum in AName2Enum.items():
                # all
                cSel = self.u.select_atoms(f"resname {rname} and name {aname}")
                weightsALL[[_a.index for _a in cSel]] = enum

        weightsLIPIDS = weightsALL[[_a.index for _a in clipids]]
        weightsWATERS = weightsALL[[_a.index for _a in cwaters]]
        print("done.")

        for ts in tqdm(self.u.trajectory, desc="Iterating over trajectory"):
            # count the index of the frame, numbered from 0, used to be used for the
            # density profile averaging
            # -- posible not needed now
            frame += 1  # ts.frame

            # reads the dimension in z-direction
            box_z = ts.dimensions[2]
            if box_z/10 < min_z:
                min_z = box_z/10

            # reads the coordinates of all of the atoms
            crds = self.u.atoms.positions

            crds_no_center = c.atoms.positions[:, 2]/10

            # calculates the center of mass of the selected atoms that the density
            # should be centered around and takes the z-coordinate value
            ctom = c.atoms.center_of_mass()[2]

            # moves the center of mass of the selected centering group into box/2
            crds[:, 2] += box_z/2 - ctom

            # shifts the coordinates in the universe by the value of
            # the center of mass
            self.u.atoms.positions = crds

            # puts the atoms back to the original box dimension; it possibly does not
            # take PBC into account, therefore it may brake some of the water molecules;
            # -- try it, come to the issue later
            self.u.atoms.pack_into_box()

            # shif the coordinates so that the center in z-dimention is in 0;
            # divide by 10 to get the coordinates in nm, since now the crds are
            # only the z coordinates
            crds = (self.u.atoms.positions[:, 2] - box_z/2)/10
            clipids_center = (clipids.atoms.positions[:, 2] - box_z/2)/10
            cwaters_center = (cwaters.atoms.positions[:, 2] - box_z/2)/10

            # calculates the volume of the bin; d- the "height" of a bin;
            # assumes in [nm]
            # ts.dimension[0], ts.dimension[1] - the x and y dimension;
            # in [A] --> devides by 100
            vbin = d*np.prod(ts.dimensions[:2])/100  # needed for density, correct!

            # calculates the total density profile; keep for now
            density_z_centered += np.histogram(
                crds, bins=self.nbin, range=(-boxH/2, boxH/2),
                weights=weightsALL / vbin,
            )[0]
            density_z_no_center += np.histogram(
                crds_no_center, bins=self.nbin, range=(0, boxH),
                weights=weightsLIPIDS / vbin,
            )[0]
            # ELECTRON WEIGHTS
            density_lipids_center += np.histogram(
                clipids_center, bins=self.nbin, range=(-boxH/2, boxH/2),
                weights=weightsLIPIDS / vbin,
            )[0]
            # ELECTRON WEIGHTS
            density_waters_center += np.histogram(
                cwaters_center, bins=self.nbin, range=(-boxH/2, boxH/2),
                weights=weightsWATERS / vbin,
            )[
                0
            ]

        print("Calculating the density takes {:10.6f} s".format(time.time()-start_time))

        """ Normalizing the profiles """
        density_z_centered /= (frame)
        density_z_no_center /= frame
        density_lipids_center /= frame
        density_waters_center /= frame

        # Post-processign data and writing to file
        density_data = np.vstack((x, density_z_centered)).transpose()
        density_data_no_center = np.vstack((x, density_z_no_center)).transpose()
        density_lipids_center = np.vstack((x, density_lipids_center)).transpose()
        density_waters_center = np.vstack((x, density_waters_center)).transpose()

        # Get the indexes of the final density data where all the time steps contribute
        # In other words, take the coordinates of the smalest box from the simulation
        final_FF_start = int(np.round(self.nbin/2 - min_z/d/2)) + 1
        final_FF_end = int(np.round(self.nbin/2 + min_z/d/2)) - 1

        FF_range = np.linspace(0, 999, 1000)
        fa_aver, fb_aver = self.fourier(
            density_data[final_FF_start:final_FF_end, 1],
            density_data[final_FF_end, 0] - density_data[final_FF_start, 0],
            FF_range,
            density_data[1, 0] - density_data[0, 0],
        )

        # Plot density profiles from the average density with minimal box
        fourrier_result2 = np.sqrt(np.multiply(fa_aver, fa_aver) +
                                   np.multiply(fb_aver, fb_aver))
        fourrier_data2 = np.vstack((FF_range*0.1*0.01, fourrier_result2)).transpose()

        # Save data into files
        # minimum box size density
        with open(str(self.output)+"TotalDensity.dat", 'wb') as f:
            np.savetxt(f, density_data[final_FF_start+1:final_FF_end-1, :],
                       fmt='%8.4f  %.8f')

        with open(str(self.output)+"LipidDensity_no_center.dat", 'wb') as f:
            np.savetxt(f, density_data_no_center[final_FF_start+1:final_FF_end-1, :],
                       fmt='%8.4f  %.8f')

        with open(str(self.output)+"LipidDensity.dat", 'wb') as f:
            np.savetxt(f, density_lipids_center[final_FF_start+1:final_FF_end-1, :],
                       fmt='%8.4f  %.8f')

        with open(str(self.output)+"WaterDensity.dat", 'wb') as f:
            np.savetxt(f, density_waters_center[final_FF_start+1:final_FF_end-1, :],
                       fmt='%8.4f  %.8f')

        # this is the important file form factors
        with open(str(self.output)+"FormFactor.dat", 'wb') as f:
            np.savetxt(f, fourrier_data2, fmt='%8.4f  %.8f')

        # write outputs in JSON

        with open(str(self.output)+"TotalDensity.json", 'w') as f:
            json.dump(density_data[final_FF_start+1:final_FF_end-1, :],
                      f, cls=NumpyArrayEncoder)

        with open(str(self.output)+"LipidDensity.json", 'w') as f:
            json.dump(density_lipids_center[final_FF_start+1:final_FF_end-1, :],
                      f, cls=NumpyArrayEncoder)

        with open(str(self.output)+"WaterDensity.json", 'w') as f:
            json.dump(density_waters_center[final_FF_start+1:final_FF_end-1, :],
                      f, cls=NumpyArrayEncoder)

        with open(str(self.output)+"FormFactor.json", 'w') as f:
            json.dump(fourrier_data2, f, cls=NumpyArrayEncoder)
    # end calculate_density

    def fourier(self, ff_density, box_z, FF_range, d_ff):
        """Calculates fourier transform of ff_density in the FF_range.
        It calculates a "height" of a bin for FF puroposes; in this case the number of
        bins is constant and the bin width changes.

        Args:
            ff_density (_type_): TODO: fill
            box_z (_type_): TODO: fill
            FF_range (_type_): TODO: fill
            d_ff (_type_): TODO: fill

        Returns:
            _type_: TODO: fill
        """

        # Creates the direct space coordinates
        # the calculations are stable with rounding (and others) errors in the direct
        # space coordinates
        ff_x = np.linspace(-box_z/2, box_z/2, ff_density.shape[0]+1)[:-1] +\
            box_z/2/ff_density.shape[0]

        k = 0
        bulk = 0
        while k*d_ff < 0.33:
            bulk += ff_density[k] + ff_density[-k-1]
            k += 1
        bulk /= (2*k)

        fa = np.zeros(FF_range.shape[0])
        fb = np.zeros(FF_range.shape[0])

        for j in range(0, ff_density.shape[0]):
            fa += (ff_density[j]-bulk)*np.cos(FF_range*ff_x[j]*0.01)*d_ff
            fb += (ff_density[j]-bulk)*np.sin(FF_range*ff_x[j]*0.01)*d_ff

        return fa, fb
