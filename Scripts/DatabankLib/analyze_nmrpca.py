#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This code aims to calculate the relaxation time based on the
PCAlipids analysis.

For details check:

[Principal Component Analysis of Lipid Molecule Conformational
Changes in Molecular Dynamics Simulations
Pavel Buslaev, Valentin Gordeliy, Sergei Grudinin, and Ivan Gushchin
Journal of Chemical Theory and Computation 2016 12 (3), 1019-1028
DOI: 10.1021/acs.jctc.5b01106 ]

and

[Principal component analysis highlights the influence of temperature,
curvature and cholesterol on conformational dynamics of lipids
Pavel Buslaev, Khalid Mustafin, and Ivan Gushchin
Biochimica et Biophysica Acta (BBA) - Biomembranes
Volume 1862, Issue 7, 1 July 2020, 183253
DOI: 10.1016/j.bbamem.2020.183253]

The main idea is to run PCA on concatenated lipid trajectory,
calculate the autocorrelation times, which are linearly related
to equilibration times of individual lipids. The parameters of linear
transform are calculated based on the trajectories from NMRlipids
databank.

The code was developed by
Dr. Pavel Buslaev, pavel.i.buslaev@jyu.fi
Dr. Samuli Olilla
Patrik Kula
Alexander Kuzmin
"""

import os
import json
import numpy as np
from scipy import signal
import MDAnalysis as mda

from DatabankLib.databankLibrary import lipids_dict, loadMappingFile
from DatabankLib.databankio import resolve_download_file_url, download_resource_from_uri

from MDAnalysis.analysis.base import AnalysisFromFunction
import warnings
# suppress some MDAnalysis warnings issued from mda.analysis.align
warnings.filterwarnings('ignore')
from MDAnalysis.analysis import align  # noqa

# TODO: now there are only regular phospholipids.
# The list should be verified by method authors.
ALLOWLIPIDS = [
    "POPC", "POPG", "POPS", "POPE", "PYPC", "DMPC", "DPPC", "DPPE", "DPPG",
    "DEPC", "DRPC", "DYPC", "DLPC", "DLIPC", "DOPC", "DOPE", "DDOPC", "DOPS",
    "DSPC", "DAPC", "SDPE", "SOPC", "POPI", "SAPI", "SAPI24", "SAPI25",
    "SLPI"
]

mergeCutoff = 2.0
trjSizeCutoff = 5000000000

TAILSN1 = "sn-1"
TAILSN2 = "sn-2"
HEADGRP = "headgroup"

# In case you want to test for a specific coordinate change the flag to yes
TEST = False


class Parser:
    """
    Class Parser is a basic class to work with the trajectory. It has basic utility
    for preparing trajectories:

    1. Checking if simulation has a correct name
    2. Checking if trajectory is downloaded. Downloading, if not
    3. Calling AmberMerger to merge head groups and tails for Amber
    trajectories
    4. Concatenate trajectories
    """

    def __init__(
            self, root, readme, eq_time_fname="eq_times.json",
            path=None, v=True):
        """
        Constructor for Parser:
            root          - path to the directory with trajectories
            readme        - simulation data
            eq_time_fname - name of the output file
            path          - name of a particular trajectory
            v             - verbosity for testing
        """

        # Technical
        self.verbose = v
        self.error = 0

        # Path
        self.root = root
        self.eq_time_fname = eq_time_fname

        # Extracting data from readme
        self.indexingPath = os.path.join(root, readme["path"])
        print('INdexing path:', self.indexingPath)
        if self.verbose:
            print(f"Parser: Processing trajectory {self.indexingPath}")
        self.doi = readme["DOI"]
        self.soft = readme["SOFTWARE"]
        if self.soft == "openMM" or self.soft == "NAMD":
            try:
                self.trj = readme["TRJ"][0][0]
                self.tpr = readme["PDB"][0][0]
            except Exception:
                print("Parser: Did not find trajectory or pdb for openMM or NAMD")
                self.error = 1
        else:
            try:
                self.trj = readme["TRJ"][0][0]
                self.tpr = readme["TPR"][0][0]
            except Exception:
                print("Parser: Did not find trajectory or tpr for GROMACS")
                self.trj = ""
                self.tpr = ""
                self.error = 1
        self.trj_name = os.path.join(self.indexingPath, self.trj)
        self.tpr_name = os.path.join(self.indexingPath, self.tpr)
        self.trjLen = readme["TRJLENGTH"] / 1000  # ns
        self.FF = readme.get("FF")

        self.size = readme["TRAJECTORY_SIZE"]

        self.composition = readme["COMPOSITION"]
        lipids_list = list(lipids_dict.keys())
        self.lipids = []
        # TODO: add to the lipids those lipids that are indicated in .yaml
        for lipid in lipids_list:
            try:
                if self.composition[lipid] != 0:
                    self.lipids.append(lipid)
            except KeyError:
                continue

        self.path = path

    def validatePath(self):
        """
        Path validation. Behaviour depends on the input.
        1. If ther were any errors, validation fails. Currently the only
        error tested is that .xtc and .tpr files are not present. This
        is the case for non-GROMACS trajectories.
        2. If path is not provided, parser is iterating over all trajectories.
        3. If path is provided and current path is equal to the provided one,
        parser reports that it finds the trajectory for the analysis.
        """
        if self.error > 0:
            # This is an error message. Printing even in silent mode
            print(
                "Parser: Can't read TPR/PDB and TRJ from README for " +
                f"{self.indexingPath}")
        if self.verbose and not self.path:
            print(
                "Parser: Iterating over all trajectories. "
                + f"Current trajectory is {self.indexingPath}"
            )
            # return
        # if self.path == self.indexingPath:
        if os.path.isfile(os.path.join(self.indexingPath, self.eq_time_fname)):
            if self.verbose:
                print("Parser: Found file with equilibration data. \n"
                      "Not processing the trajectory")
                return -1
        if self.verbose:
            print(f"Parser: Found trajectory {self.indexingPath}")
        return 0

    def downloadTraj(self):
        """
        Basic trajectory and TPR download.
        """
        print("Downloading")
        if not os.path.isfile(self.tpr_name):
            self.tpr_url = resolve_download_file_url(self.doi, self.tpr)
            # This is a log message. Printing even in silent mode
            print("Parser: Downloading tpr ", self.doi)
            # urllib.request.urlretrieve(self.tpr_url, self.tpr_name)
            download_resource_from_uri(self.tpr_url, self.tpr_name)
        if not os.path.isfile(self.trj_name):
            self.trj_url = resolve_download_file_url(self.doi, self.trj)
            # This is a log message. Printing even in silent mode
            print("Parser: Downloading trj ", self.doi)
            # urllib.request.urlretrieve(self.trj_url, self.trj_name)
            download_resource_from_uri(self.trj_url, self.trj_name)

    def prepareGMXTraj(self):
        """
        Preparing trajectory. If centered trajectory is found, use it. If whole
        trajectory is found, use it. Otherwise, call gmx trjconv to make whole
        trajectory. The selected trajectory is loaded into Universe.
        """
        # Look for centered.xtc
        if os.path.isfile(os.path.join(self.indexingPath, "centered.xtc")):
            print("Parser: Founder centered trajectory: centered.xtc")
            self.trj_name = os.path.join(self.indexingPath, "centered.xtc")
        # Look for whole.xtc
        elif os.path.isfile(os.path.join(self.indexingPath, "whole.xtc")):
            print("Parser: Founder whole trajectory: whole.xtc")
            self.trj_name = os.path.join(self.indexingPath, "whole.xtc")
        # Run gmx trjconv
        else:
            print("Parser: Making trajectory whole")
            trj_out_name = os.path.join(self.indexingPath, "whole.xtc")
            os.system(
                f"echo System | gmx trjconv -f {self.trj_name} -s "
                f"{self.tpr_name} -pbc mol -o {trj_out_name}"
            )
            self.trj_name = trj_out_name

        # Skip frames for large trajectories
        if self.size > trjSizeCutoff:
            trj_out_name = os.path.join(self.indexingPath, "short.xtc")
            os.system(
                f"echo System | gmx trjconv -f {self.trj_name} -s "
                f"{self.tpr_name} -pbc mol -o {trj_out_name} -skip 10"
            )
            self.trj_name = trj_out_name

        # Create .gro file
        self.gro_name = os.path.join(self.indexingPath, "conf.gro")
        if os.path.isfile(self.gro_name):
            print("Parser: Found gro file")
        else:
            print("Parser: Making a gro file from the first frame")
            os.system(
                f"echo System | gmx trjconv -s {self.tpr_name} "
                f" -f {self.trj_name} -dump 0 -o {self.gro_name}"
            )

        # Loading trajectory
        try:
            self.traj = mda.Universe(self.tpr_name, self.trj_name)
        except Exception:
            self.traj = mda.Universe(self.gro_name, self.trj_name)

    # TODO: not tested!!
    def prepareOpenMMTraj(self):
        print("openMM or NAMD")
        if self.size > trjSizeCutoff:
            trj_out_name = os.path.join(self.indexingPath, "short.xtc")
            if os.path.isfile(trj_out_name):
                print("Parser: Short trajectory is found")
            else:
                u = mda.Universe(self.tpr_name, self.trj_name)
                with mda.Writer(trj_out_name,
                                u.select_atoms("all").n_atoms) as W:
                    for ts in u.trajectory[::10]:
                        W.write(u.select_atoms("all"))

            self.trj_name = trj_out_name

        self.traj = mda.Universe(self.tpr_name, self.trj_name)

    def prepareTraj(self):
        if self.soft == "openMM" or self.soft == "NAMD":
            self.prepareOpenMMTraj()
        else:
            self.prepareGMXTraj()

    def concatenateTraj(self):
        """
        Create Concatenator and corresponding concatenated trajectories for
        all lipids available for the trajectory.
        """

        self.concatenated_trajs = []
        for lipid in self.lipids:
            if lipid not in ALLOWLIPIDS:
                # We do not treat cholesterols
                continue
            topology = Topology(
                self.FF,
                self.traj,
                self.composition[lipid]["NAME"],
                self.composition[lipid]["MAPPING"]
            )
            concatenator = Concatenator(topology,
                                        self.traj,
                                        lipid,
                                        self.composition[lipid]["NAME"])
            self.concatenated_trajs.append(concatenator.concatenate())

    def dumpData(self, data):
        """
        Write data to json file
        """
        with open(os.path.join(self.indexingPath,
                               self.eq_time_fname), "w") as f:
            json.dump(data, f)


class Topology:
    """
    Class Topology is a class needed to extract lipid specific data and also to
    merge lipids from Amber trajectories, where lipid is often represented as 3
    residue types. It has the following methods:

    1. Simple constructor, which sets the force field,residue names, trajectory,
    and loads mapping_file
    2. Mapping file loader
    3. Function that outputs atomNames
    4. Checker if the merge is needed. Currently defunct
    5. Runner, that returns lists for head, tail1, tail2 orders for merging

    Currently defunct
    """

    def __init__(self, ff, traj, lipid_resname, mapping_file):
        """
        Constructor for Topology:
            ff            - force field name
            traj          - MDAnalysis trajectory
            lipid_resname - name of lipid
            mapping_file  - path to maping file
        """

        self.ff = ff
        self.lipid_resname = lipid_resname
        self.traj = traj
        self.mapping = loadMappingFile(mapping_file)

    def atomNames(self):
        """
        Extract all names of heavy atoms from the mapping
        """
        atoms = []
        for key in self.mapping:
            atom = self.mapping[key]["ATOMNAME"]
            if atom[0] == "H":
                continue
            atoms.append(atom)
        return atoms

    def isMergeNeeded(self):
        """
        Checker for merge. Currently it checks if the RESIDUE key is in the
        mapping file.
        NOTE: This has to be changed if the content of mapping file changes
        """

        if "RESIDUE" not in self.mapping[list(self.mapping.keys())[0]].keys():
            return False
        resnames = []
        for key in self.mapping:
            atom = self.mapping[key]["ATOMNAME"]
            atomNames = self.atomNames()
            if atom in atomNames:
                resnames.append(self.mapping[key]["RESIDUE"])
        resnames = set(resnames)
        if len(resnames) > 1:
            return resnames
        return False

    def getLipidResnames(self):
        """
        Helper function that gets the residue names for lipid, if merge is needed
        """
        resnames = self.isMergeNeeded()
        if self.isMergeNeeded():
            return resnames
        # How can we end up here?
        # Since currently we only call this function when Merge is needed,
        # this is a place for checking if everything is ok and we can raise
        else:
            return self.lipid_resname

    def assignResnames(self, resnames):
        """
        Helper function that finds head, sn-1 and sn-2 tails
        NOTE: currently there are only lipids with 2 different tails available in
        databank: e.g. POPC or POPG. This leads to different names of tails. This
        won't be the case for DOPC. Currently the algorithm is using that all three
        groups differ
        """
        if self.isMergeNeeded():
            resnameDict = {}
            # First find headgroup
            for key in self.mapping:
                resname = self.mapping[key]["RESIDUE"]
                if self.mapping[key]["FRAGMENT"] == HEADGRP:
                    resnameDict[HEADGRP] = resname
                    break
            for key in self.mapping:
                resname = self.mapping[key]["RESIDUE"]
                if (
                    TAILSN1 not in resnameDict
                    and TAILSN1 == self.mapping[key]["FRAGMENT"]
                    and not resnameDict[HEADGRP] == resname
                ):
                    resnameDict[TAILSN1] = resname
                if (
                    TAILSN2 not in resnameDict
                    and TAILSN2 == self.mapping[key]["FRAGMENT"]
                    and not resnameDict[HEADGRP] == resname
                ):
                    resnameDict[TAILSN2] = resname
                if TAILSN1 in resnameDict and TAILSN2 in resnameDict:
                    break
            # TODO: add check that all resnames from input are in the dict
            return resnameDict
        # How can we end up here?
        else:
            return None

    def runMerger(self):
        """
        Find lipid tails that correspond to a particular head-group.
        NOTE: currently there are only lipids with 2 different tails available in
        databank: e.g. POPC or POPG. This leads to different names of tails. This
        won't be the case for DOPC. Currently the algorithm is using that all three
        groups differ
        TODO: get the correspondence from structure
        """

        resnames = self.getLipidResnames()
        resnameDict = self.assignResnames(resnames)
        head_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resnameDict[HEADGRP]}"
            ).residues
        ]
        sn_1_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resnameDict[TAILSN1]} and "
                + f"around {mergeCutoff} (resname {resnameDict[HEADGRP]} "
                + "and not name H*)"
            ).residues
        ]
        sn_2_residues = [
            r.atoms.select_atoms("not name H* and prop mass > 0.8")
            for r in self.traj.select_atoms(
                f"not name H* and resname {resnameDict[TAILSN2]} and "
                + f"around {mergeCutoff} (resname {resnameDict[HEADGRP]} "
                + "and not name H*)"
            ).residues
        ]
        return head_residues, sn_1_residues, sn_2_residues


class Concatenator:
    """
    Class Concatenator is a class needed to concatenate trajectory for lipid types.
    It has the following methods:
    1. Simple constructor, which sets the topology,residue names, and
    trajectory
    2. concatenateTraj method to do the basic by lipid concatenation
    3. alignTraj method to perform the alignment of the concatenated trajectory
    4. The enveloping concatenate method
    """

    def __init__(self, topology, traj, lipid_name, lipid_resname):
        """
        Constructor for Concatenator:
            topology      - topology for lipid
            traj          - MDAnalysis trajectory
            lipid_name    - lipid name in the databank
            lipid_resname - lipid resname in the trajectory
        """
        self.topology = topology
        self.traj = traj
        self.lipid_name = lipid_name
        self.lipid_resname = lipid_resname

        if self.topology.isMergeNeeded():
            self.headlist, self.tail1list, self.tail2list = \
                self.topology.runMerger()
            # TODO: raise if length of 3 lists is not equal
        else:
            self.headlist = None
            self.tail1list = None
            self.tail2list = None

    def concatenateTraj(self):
        """
        concatenateTraj performs basic trajectory concatination. First, it extracts
        coordinates from trajectory, next, it reshapes the coordinate array, swaps
        time and resid axes to obtain continuous trajectories of individual lipids
        (this is needed for autocorrelation time analysis), and finally merges
        individual lipid trajectories.
        """
        traj = self.traj.trajectory
        n_frames = len(traj)

        heavy_atoms_topology = self.traj.select_atoms(
            f"resname {self.lipid_resname} and not name H*"
        )

        n_atoms_lipid = len(self.topology.atomNames())

        n_lipid = heavy_atoms_topology.n_residues
        if n_lipid * n_atoms_lipid != heavy_atoms_topology.n_atoms:
            n_lipid = heavy_atoms_topology.n_atoms // n_atoms_lipid

        # Get all coordinates n_frames, n_lipid * n_atoms_lipid
        coords = (
            AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                 heavy_atoms_topology)
            .run()
            .results["timeseries"]
        )
        # Adding axis to separate individual lipid trajectories
        coords = coords.reshape((n_frames, n_lipid, n_atoms_lipid, 3))
        # Swapping time frame with lipid axis
        # Now we have continuous lipid trajectory
        coords = np.swapaxes(coords, 0, 1)
        # Reshaping to desired shape
        coords = coords.reshape((n_frames * n_lipid, n_atoms_lipid, 3))

        # Creating new Universe from the concatenated coordinates
        atom_resindex = np.repeat(0, n_atoms_lipid)
        concatenated_traj = mda.Universe.empty(
            n_atoms_lipid, 1, atom_resindex=atom_resindex, trajectory=True
        )
        concatenated_traj.add_TopologyAttr("resname", [self.lipid_resname])
        concatenated_traj.load_new(coords, in_memory=False)

        return concatenated_traj, n_lipid, n_frames * n_lipid

    def concatenateTrajWithMerging(self):
        """
        concatenateTrajWithMerging performs basic trajectory concatination. In
        contrast to basic concatenateTraj it additionally merges splitted lipids.
        First, it creates extracts coordinates from trajectory, next, it reshapes
        the coordinate array, swaps time and resid axes to obtain continuous
        trajectories of individual lipids (this is needed for autocorrelation
        time analysis), and finally merges individual lipid trajectories.
        """
        traj = self.traj.trajectory
        n_frames = len(traj)

        heavy_atoms_topology = \
            self.headlist[0] + \
            self.tail1list[0] + \
            self.tail2list[0]

        for head, sn1, sn2 in zip(
            self.headlist[1:], self.tail1list[1:], self.tail2list[1:]
        ):
            heavy_atoms_topology = heavy_atoms_topology.union(head + sn1 + sn2)

        n_atoms_lipid = len(self.topology.atomNames())

        n_lipid = heavy_atoms_topology.n_atoms // n_atoms_lipid

        # TODO: add check
        # n_atoms_lipid == heavy_atoms_topology.n_atoms

        # Get all coordinates n_frames, n_lipid * n_atoms_lipid
        coords = (
            AnalysisFromFunction(lambda ag: ag.positions.copy(),
                                 heavy_atoms_topology)
            .run()
            .results["timeseries"]
        )
        # Adding axis to separate individual lipid trajectories
        coords = coords.reshape((n_frames, n_lipid, n_atoms_lipid, 3))
        # Swapping time frame with lipid axis
        # Now we have continuous lipid trajectory
        coords = np.swapaxes(coords, 0, 1)
        # Reshaping to desired shape
        coords = coords.reshape((n_frames * n_lipid, n_atoms_lipid, 3))

        # Creating new Universe from the concatenated coordinates
        atom_resindex = np.repeat(0, n_atoms_lipid)
        concatenated_traj = mda.Universe.empty(
            n_atoms_lipid, 1, atom_resindex=atom_resindex, trajectory=True
        )
        concatenated_traj.add_TopologyAttr("resname", [self.lipid_resname])
        concatenated_traj.load_new(coords, in_memory=False)

        return concatenated_traj, n_lipid, n_frames * n_lipid

    """
    alignTraj alignes the concatenated trajectory in two steps: (1) it computes
    average structure after alignment to the first frame, and (2) it alignes
    the structure to the calculated average structure in (1).
    """
    def alignTraj(self, concatenated_traj):
        # Compute average structure after alignment to the first frame
        av = align.AverageStructure(concatenated_traj, ref_frame=0).run()
        # Align to average structure
        align.AlignTraj(concatenated_traj, av.results.universe).run()
        # Compute average structure after second alignment
        coords = (
            AnalysisFromFunction(
                lambda ag: ag.positions.copy(),
                concatenated_traj.select_atoms("all")
            )
            .run()
            .results["timeseries"]
        )
        return coords.reshape(-1, 3), coords.mean(axis=0).reshape(1, -1)

    """
    Simple enveloping function to perform concatenation
    """
    def concatenate(self):
        print(f"Concatenator: Concatenating lipid {self.lipid_name}")
        if not self.topology.isMergeNeeded():
            # Merging is not needed
            concatenated_traj, n_lipid, n_frames = self.concatenateTraj()
        else:
            concatenated_traj, n_lipid, n_frames = \
                self.concatenateTrajWithMerging()
        aligned_traj, av_pos = self.alignTraj(concatenated_traj)

        return aligned_traj, av_pos, n_lipid, n_frames, self.lipid_name


class PCA:
    """
    Class PCA is a class that actually performs PCA. It has the following methods:
    1. Simple constructor, which sets the aligned trajtory, average coordinates,
    number of lipids, number of frames in the concatenated trajectory and
    trajectory length in ns
    2. PCA prepares the trajectory coordinates for analysis and calculates
    principal components
    3. get_proj projects the trajectory on principal components
    4. get_lipid_autocorrelation calculates the autocorrelation timeseries for
    individual lipid
    5. get_autocorrelations calculates the autocorrelation timeseries for
    trajectory
    """

    def __init__(self, aligned_traj, av_pos, n_lipid, n_frames, traj_time):
        """
        Constructor for PCA:
            aligned_traj - np.array with positions of concatenated and aligned
                        trajectory
            av_pos       - np.array with average positions for the lipid
            n_lipid      - number of lipids of particular type in the system
            n_frames     - number of frames in the concatenated trajectory
            traj_time    - trajectory length in ns
        """
        self.aligned_traj = aligned_traj
        self.av_pos = av_pos
        self.n_lipid = n_lipid
        self.n_frames = n_frames
        self.traj_time = traj_time

    def PCA(self):
        """
        PCA calculates the PCA. First the data is centered and then covariance
        matrix is calculated manually.
        """
        # centering of positions relative to the origin
        X = self.aligned_traj.astype(np.float64)
        X = X.reshape(self.n_frames, self.av_pos.shape[1]) - self.av_pos
        # the sum of all coordinates (to calculate mean)
        X_1 = X.sum(axis=0)
        # production of X and X^T
        X_X = np.tensordot(X, X, axes=(0, 0))
        # covariance matrix calculation
        cov_mat = (
            X_X
            - np.dot(X_1.reshape(len(X_1), 1), (X_1.reshape(len(X_1), 1)).T)
            / self.n_frames
        ) / (self.n_frames - 1)
        # eigenvalues and eigenvectors calculation
        eig_vals, eig_vecs = np.linalg.eigh(cov_mat)
        self.eig_vecs = np.flip(eig_vecs, axis=1).T

        return X

    def get_proj(self, cdata):
        """
        Projecting the trajectory on the 1st principal component
        """
        # projection on PC1
        proj = np.tensordot(cdata, self.eig_vecs[0:1], axes=(1, 1)).T
        proj = np.concatenate(np.array(proj), axis=None)

        self.proj = proj

    def get_lipid_autocorrelation(self, data, variance, mean):
        """
        Autocorrelation calculation for individual lipid.
        """
        # Centering the data
        data -= mean
        # Convolve the data
        r = signal.fftconvolve(data, data[::-1], mode="full")[-len(data):]
        # Weight the correlation to get the result in range -1 to 1
        return r / (variance * (np.arange(len(data), 0, -1)))

    def get_autocorrelations(self):
        """
        Autocorrelation calculation for the trajectory.
        """
        variance = self.proj.var()
        mean = self.proj.mean()
        # extract the trajectories for individual lipids
        separate_projs = [
            self.proj[
                len(self.proj) * i // self.n_lipid:len(self.proj)
                * (i + 1) // self.n_lipid]
            for i in range(self.n_lipid)
        ]
        # calculate autocorrelations for individual lipids
        R = np.array(
            [self.get_lipid_autocorrelation(x, variance, mean)
                for x in separate_projs]
        )
        R = R.mean(axis=0)
        T = np.arange(len(R)) * self.traj_time / len(R)
        self.autocorrelation = np.array([T, R]).T


class TimeEstimator:
    """
    Class TimeEstimator is a class that estimates equilbration time from
    autocorrelation data. It includes the following methods:
        1. Simple constructor, which sets the autocorrelation data
        2. Helper method get_nearest_value that finds the interval in the data
        where the autocorrelation falls below specific value
        3. timerelax method calculates the logarithms of autocorrelation and
        calculates the decay time
        4. calculate_time method calculates the estimated equilibration time
    """

    def __init__(self, autocorrelation):
        """
        Constructor for PCA:
            autocorrelation - autocorrelation data
        """
        self.autocorrelation = autocorrelation

    def get_nearest_value(self, iterable, value):
        """
        get_nearest_value method return the indices that frame the particular value
        As an input it gets an array, which is expected to decay. The method tries
        to find and index (index_A), for which the array gets below the cutoff value
        for the first time. Then, for index_A-1 the array was larger than the cutoff
        value for the last time. It means that the function, represented by the
        array data is equal to cutoff somewhere between timesteps that correspond to
        index_A-1 and index_A. It can happen, that this index_A is equal 0. Then,
        method searches for the first occurance, where array is smaller than
        array[0]. The method returns the interval inbetween the array gets below
        cutoff.
        """
        try:
            A = np.where(iterable < value)[0][0]
        except Exception:
            print("TimeEstimator: Autocorrelations do not converge. "
                  "We shift to extrapolation regime.")
            A = np.where(iterable == np.min(iterable))[0][0]
            return A, A - 1

        if A == 0:
            B = np.where(iterable < iterable[A])[0][0]
        else:
            B = A - 1
        return A, B

    def timerelax(self):
        """
        Estimates autocorrelation decay time.
        """
        time = self.autocorrelation[:, 0]
        autocorrelation = self.autocorrelation[:, 1]

        points = []
        j = 1
        while j < len(autocorrelation):
            points.append(j)
            j = int(1.5 * j) + 1
        # Get timeseries for autocorrelation
        T_pic = np.array([time[i] for i in points if autocorrelation[i] > 0.0])
        R_pic = np.array(
            [autocorrelation[i] for i in points if autocorrelation[i] > 0.0]
        )
        # Calculate logs for time and autocorrelation arrays
        # We use log since in log-log scale autocorrelation is closer to linear
        R_log = np.log(R_pic)
        T_log = np.log(T_pic)
        # data interpolation. We are searching for time interval where
        # autocorrelations decay by e. This is the most stable method.
        # Autocorrelations decay by e is equivalent to
        # log(autocorrelation) < - 1
        power = -1
        A, B = self.get_nearest_value(R_log, power)
        # perform interpolation
        a = (T_log[B] - T_log[A]) / (R_log[B] - R_log[A])
        b = T_log[A] - a * R_log[A]
        t_relax1 = a * power + b
        T_relax1 = np.exp(t_relax1)

        return T_relax1

    def calculate_time(self):
        """
        Basic enveloping method that estimates equilibration time from
        autocorrelation decay time. They are linearly connected and the
        coefficient is calculated experimentally.
        """
        # relaxation time at e^1 decay
        te1 = self.timerelax()

        return te1 * 49  # 49 from sample data average of tkss2/tac1
