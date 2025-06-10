"""
:module: DatabankLib.analyze
:description: major analysis methods calculating a property from NMRLipipds
              `system` object
"""

import json
import os
import sys
import re
import traceback
from logging import Logger
import buildh
import urllib.request
import socket

from tqdm import tqdm
import numpy as np
import gc

from DatabankLib import (
    NMLDB_SIMU_PATH, RCODE_ERROR, RCODE_SKIPPED, RCODE_COMPUTED,
    NMLDB_DATA_PATH)
from DatabankLib.core import System
from DatabankLib.settings.molecules import lipids_set
from DatabankLib.settings.engines import get_struc_top_traj_fnames
from DatabankLib.databankLibrary import (
    GetNlipids, system2MDanalysisUniverse)
from DatabankLib.jsonEncoders import CompactJSONEncoder
from DatabankLib.databankio import resolve_download_file_url
from DatabankLib.databankop import find_OP
from DatabankLib.form_factor import FormFactor
from DatabankLib import analyze_nmrpca as nmrpca


def computeNMRPCA(  # noqa: N802 (API)
        system: System, logger: Logger, recompute: bool = False) -> int:
    """Compute eq_times.json using NMR PCA analysis.

    Args:
        system (System): one of systems of the Databank
        recompute (bool, optional): Delete previous apl.json and recompute it if True.
        Defaults to False.
    Returns:
        int success code (RCODE_...)
    """
    print('== Doi: ' + system['DOI'])
    print('== System name: ' + system['SYSTEM'])
    # getting data from databank and preprocessing them
    # Start Parser
    # TODO: 2test|    parser = Parser(NMLDB_SIMU_PATH, readme, eq_time_fname, testTraj)
    parser = nmrpca.Parser(system, 'eq_times.json')
    # Check trajectory
    print(parser._path)
    vpcode = parser.validate_path()
    print("ValidatePath code: ", vpcode)
    if vpcode > 0:
        sys.stderr.write("Some errors in analyze_nmrpca.py::Parser constructor.\n")
        return RCODE_ERROR
    elif vpcode < 0 and not recompute:
        return RCODE_SKIPPED

    if (
        'WARNINGS' in system and
        system['WARNINGS'] is not None and
        'AMBIGUOUS_ATOMNAMES' in system['WARNINGS']
    ):
        return RCODE_SKIPPED

    if (
        'WARNINGS' in system and
        system['WARNINGS'] is not None and
        'SKIP_EQTIMES' in system['WARNINGS']
    ):
        return RCODE_SKIPPED

    # Check if TPR is defined and non-null or ""
    try:
        __ = system['TPR'][0][0]
        _ = __[0]
    except (KeyError, TypeError):
        sys.stderr.write("TPR is required for NMRPCA analysis!")
        return RCODE_ERROR

    # Download files
    parser.download_traj()
    # Prepare trajectory
    parser.prepare_traj()
    # Concatenate trajectory
    parser.concatenate_traj()
    equilibration_times = {}
    # Iterate over trajectories for different lipids
    for lip, traj in parser.concatenated_trajs.items():
        print(f"Main: parsing lipid {lip}")
        # Creat PCA for trajectory
        pca_runner = nmrpca.PCA(traj[0], traj[1], traj[2], traj[3], parser.trj_len)
        # Run PCA
        data = pca_runner.PCA()
        print("Main: PCA: done")
        # Project trajectory on PC1
        pca_runner.get_proj(data)
        print("Main: Projections: done")
        # Calculate autocorrelations
        pca_runner.get_autocorrelations()
        print("Main: Autocorrelations: done")
        # Estimate equilibration time
        te2 = nmrpca.TimeEstimator(pca_runner.autocorrelation).calculate_time()
        equilibration_times[lip] = te2 / parser.trj_len
        print("Main: EQ time: done")

        print(te2 / parser.trj_len)

    parser.dump_data(equilibration_times)
    gc.collect()
    return RCODE_COMPUTED


def computeAPL(  # noqa: N802 (API)
        system: System, logger: Logger, recompute: bool = False) -> int:
    """Generate apl.json analysis file for a system.

    Args:
        system (dict): one of systems of the Databank
        recompute (bool, optional): Delete previous apl.json and recompute it if True.
        Defaults to False.
    Returns:
        int success code (RCODE_...)
    """
    # TODO: reading software and file path for simulations
    # software = system['SOFTWARE']
    path = system['path']

    # this is checking if area per lipid is already calculated for the systems
    outfilename = os.path.join(NMLDB_SIMU_PATH,  path, 'apl.json')
    if os.path.isfile(outfilename):
        if recompute:
            os.unlink(outfilename)
        else:
            return RCODE_SKIPPED

    print('Analyzing: ', path)
    print('Will write into: ', outfilename)

    # calculates the total number of lipids
    n_lipid = GetNlipids(system)

    # makes MDAnalysis universe from the system. This also downloads the data if not
    # yet locally available
    u = system2MDanalysisUniverse(system)

    if u is None:
        print('Generation of MDAnalysis universe failed in folder', path)
        return RCODE_ERROR

    # this calculates the area per lipid as a function of time and stores it
    # in the databank
    apl = {}
    for ts in tqdm(u.trajectory, desc='Scanning the trajectory'):
        if u.trajectory.time >= system['TIMELEFTOUT']*1000:
            dims = u.dimensions
            apl_frame = dims[0]*dims[1]*2/n_lipid
            apl[u.trajectory.time] = apl_frame

    with open(outfilename, 'w') as f:
        json.dump(apl, f, cls=CompactJSONEncoder)

    return RCODE_COMPUTED


def computeTH(  # noqa: N802 (API)
        system: dict, logger: Logger, recompute: bool = False) -> int:
    thick_fn = os.path.join(NMLDB_SIMU_PATH, system['path'], 'thickness.json')
    if os.path.isfile(thick_fn) and not recompute:
        # file exist. Skipping
        return RCODE_SKIPPED

    wat_dens_name = os.path.join(
        NMLDB_SIMU_PATH, system['path'], 'WaterDensity.json')
    lip_dens_name = os.path.join(
        NMLDB_SIMU_PATH, system['path'], 'LipidDensity.json')
    print(lip_dens_name)
    try:
        with open(wat_dens_name) as f:
            wat_dens = json.load(f)
        with open(lip_dens_name) as f:
            lip_dens = json.load(f)
        wd = np.array(wat_dens)
        ld = np.array(lip_dens)
        idx = np.argwhere(np.diff(np.sign(wd[:, 1] - ld[:, 1]))).flatten()
        if len(idx) < 2:
            print("Dehydrated sample! Standard thickness rule doesn't work."
                  " Will extract boxZ.")
            thickness = wd[-1, 0] - wd[0, 0]
        else:
            thickness = wd[idx[1], 0] - wd[idx[0], 0]
        with open(thick_fn, 'w') as f:
            json.dump(thickness, f)
    except Exception as e:
        print('Calculation failed for ' + system['path'])
        print(str(e))
        print(traceback.format_exc())
        return RCODE_ERROR

    return RCODE_COMPUTED


# TODO: implement onlyLipid
def computeOP(  # noqa: N802 (API)
        system: System, logger: Logger, recompute: bool = False) -> int:
    """_summary_

    Args:
        system (dict): _description_
        recompute (bool, optional): _description_. Defaults to False.

    Returns:
        int: _description_
    """
    path = system['path']

    # Check if order parameters are calculated or something in the system prevents
    # order parameter calculation
    for key in system['COMPOSITION']:
        outfilename = os.path.join(NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')
        # print(outfilename)
        if (
            os.path.isfile(outfilename) or
            (
                'WARNINGS' in system.keys() and
                type(system['WARNINGS']) is dict and
                'AMBIGUOUS_ATOMNAMES' in system['WARNINGS'].keys() and
                key in system['WARNINGS']['AMBIGUOUS_ATOMNAMES']
            )
        ):
            file_found = True
        elif key in lipids_set:
            file_found = False
            break

    if file_found and not recompute:
        return RCODE_SKIPPED

    # If order parameters are not calculated and system ok, then continue to
    # calculating order parameters
    print('Analyzing: ', path)

    # Download the simulations and create MDAnalysis universe (the universe is not used
    # here but this script downloads the data)
    _ = system2MDanalysisUniverse(system)

    # Software and time for equilibration period
    software = system['SOFTWARE']
    eq_time = float(system['TIMELEFTOUT'])*1000

    # Check if all or united atom simulation
    try:
        united_atom = system['UNITEDATOM_DICT']
    except (KeyError):
        united_atom = False

    # Check relevant warnings
    g3switch = (
        'WARNINGS' in system and
        type(system['WARNINGS']) is dict and
        'GROMACS_VERSION' in system['WARNINGS'] and
        system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3'
    )
    trconv_command = 'trjconv' if g3switch else 'gmx trjconv'

    cur_path = os.path.join(NMLDB_SIMU_PATH, path)

    try:
        struc_fname, top_fname, trj_fname = \
            get_struc_top_traj_fnames(system, join_path=cur_path)
    except (ValueError, KeyError) as e:
        sys.stderr.write("Error reading filenames from system dictionary.\n")
        sys.stderr.write(str(type(e)) + " => " + str(e))
        return RCODE_ERROR

    try:
        # Set topology file names and make Gromacs trajectories whole
        if 'gromacs' in software:
            if top_fname is None:
                raise ValueError("TPR is required for OP calculations!")

            xtcwhole = os.path.join(NMLDB_SIMU_PATH, path, 'whole.xtc')
            if not os.path.isfile(xtcwhole):
                exec_str = (
                    f"echo System | {trconv_command} -f {trj_fname} "
                    f"-s {top_fname} -o {xtcwhole} -pbc mol -b {str(eq_time)}"
                )
                print("Make molecules whole in the trajectory")
                if united_atom and system['TRAJECTORY_SIZE'] > 15e9:
                    print("United atom trajectry larger than 15 Gb. "
                          "Using only every third frame to reduce memory usage.")
                    exec_str += " -skip 3"
                rcode = os.system(exec_str)
                if rcode != 0:
                    raise RuntimeError("trjconv exited with error (see above)")
        elif 'openMM' in software or 'NAMD' in software:
            if not os.path.isfile(struc_fname):
                pdb_url = resolve_download_file_url(system.get('DOI'), struc_fname)
                _ = urllib.request.urlretrieve(pdb_url, struc_fname)
        else:
            print("Order parameter calculation for other than gromacs, "
                  "openMM and NAMD are yet to be implemented.")
            return RCODE_ERROR

        # Calculate order parameters
        # -----------------------------------
        # united-atom switch -> engine switch
        if united_atom:
            if 'gromacs' not in software:
                raise ValueError("UNITED_ATOMS is supported only for GROMACS engine!")
            frame0struc = os.path.join(NMLDB_SIMU_PATH, path, 'frame0.gro')
            if g3switch:
                rcode = os.system(
                    f'echo System | editconf -f {top_fname} -o {frame0struc}')
                if rcode != 0:
                    raise RuntimeError("editconf exited with error (see above)")
            else:
                rcode = os.system(f"echo System | {trconv_command} -f {xtcwhole}"
                                  f" -s {top_fname} -dump 0 -o {frame0struc}")
                if rcode != 0:
                    raise RuntimeError(
                        f"trjconv ({trconv_command}) exited with error (see above)")

            for key in system['UNITEDATOM_DICT']:
                # construct order parameter definition file for CH bonds from
                # mapping file
                mapping_dict = system.content[key].mapping_dict

                def_fname = os.path.join(NMLDB_SIMU_PATH, path, key + '.def')
                def_file = open(def_fname, 'w')

                previous_line = ""

                regexp1_H = re.compile(r'M_[A-Z0-9]*C[0-9]*H[0-9]*_M')  # noqa: N806
                regexp2_H = re.compile(r'M_G[0-9]*H[0-9]*_M')  # noqa: N806
                regexp1_C = re.compile(r'M_[A-Z0-9]*C[0-9]*_M')  # noqa: N806
                regexp2_C = re.compile(r'M_G[0-9]_M')  # noqa: N806

                # Note that mapping_dict's keys must be ordered
                # C11 H111 H112 C12 H121 H122 ...
                # otherwize algorithm will fail
                for mapping_key in mapping_dict:
                    if regexp1_C.search(mapping_key) or regexp2_C.search(mapping_key):
                        atom_c = [mapping_key, mapping_dict[mapping_key]['ATOMNAME']]
                        atom_h = []
                    elif regexp1_H.search(mapping_key) or regexp2_H.search(mapping_key):
                        atom_h = [mapping_key, mapping_dict[mapping_key]['ATOMNAME']]
                    else:
                        atom_c = []
                        atom_h = []

                    if atom_h:
                        assert atom_c
                        items = [atom_c[1], atom_h[1], atom_c[0], atom_h[0]]
                        def_line = (items[2] + "&" + items[3] + " " +
                                    key + " " + items[0] + " " + items[1] + "\n")
                        if def_line != previous_line:
                            def_file.write(def_line)
                            previous_line = def_line
                def_file.close()

                # Add hydrogens to trajectory and calculate order parameters with buildH
                op_file = os.path.join(
                    NMLDB_SIMU_PATH, path, key + 'OrderParameters.dat')

                lipid_json_file = [
                    os.path.join(NMLDB_DATA_PATH, 'lipid_json_buildH',
                                 system['UNITEDATOM_DICT'][key] + '.json')
                ]

                if not os.path.isfile(lipid_json_file[0]):
                    lipid_json_file = None

                print(system['UNITEDATOM_DICT'][key])
                buildh.launch(coord_file=frame0struc, def_file=def_fname,
                              lipid_type=system['UNITEDATOM_DICT'][key],
                              lipid_jsons=lipid_json_file, traj_file=xtcwhole,
                              out_file=f"{op_file}.buildH", ignore_CH3s=True)

                outfile = open(op_file, 'w')
                outfile.write("Atom     Average OP     OP stem\n")

                data = {}
                outfile2 = os.path.join(
                    NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')

                with open(op_file + '.buildH') as op_file:
                    lines = op_file.readlines()
                    for line in lines:
                        if "#" in line:
                            continue
                        line2 = (
                            line.split()[0].replace('&', ' ') + "  " +
                            line.split()[4] + "  " + line.split()[5] +
                            " " + line.split()[6] + "\n")
                        outfile.write(line2)

                        op_name = line.split()[0].replace('&', ' ')
                        # -- line.split()[0] + " " + line.split()[1]
                        op_values = [
                            float(line.split()[4]),
                            float(line.split()[5]),
                            float(line.split()[6])
                        ]
                        data[str(op_name)] = []
                        data[str(op_name)].append(op_values)

                with open(outfile2, 'w') as f:
                    json.dump(data, f, cls=CompactJSONEncoder)

                outfile.close()

        # not united-atom cases
        else:
            if 'gromacs' in software:
                gro = os.path.join(NMLDB_SIMU_PATH, path, 'conf.gro')

                print("\n Making gro file")
                if g3switch:
                    rcode = os.system(f'echo System | editconf -f {top_fname} -o {gro}')
                    if rcode != 0:
                        raise RuntimeError("editconf exited with error (see above)")
                else:
                    rcode = os.system(f"echo System | gmx trjconv "
                                      f"-f {trj_fname} -s {top_fname} -dump 0 -o {gro}")
                    if rcode != 0:
                        raise RuntimeError("trjconv exited with error (see above)")

            for key in system['COMPOSITION']:

                if (
                    'WARNINGS' in system.keys() and
                    system['WARNINGS'] is not None and
                    'AMBIGUOUS_ATOMNAMES' in system['WARNINGS'].keys() and
                    key in system['WARNINGS']['AMBIGUOUS_ATOMNAMES']
                ):
                    print(path, key)
                    print("Order parameters cannot be calculated if atom "
                          "names are ambiguous.")
                    continue

                if key in lipids_set:
                    print('Calculating ', key, ' order parameters')
                    resname = system['COMPOSITION'][key]['NAME']
                    outfilename = os.path.join(
                        NMLDB_SIMU_PATH, path, key + 'OrderParameters.dat')
                    outfilename2 = os.path.join(
                        NMLDB_SIMU_PATH, path, key + 'OrderParameters.json')
                    if os.path.isfile(outfilename2):
                        print('Order parameter file already found')
                        continue
                    if 'gromacs' in software:
                        try:
                            op_obj = find_OP(system.content[key].mapping_dict,
                                             top_fname, xtcwhole, resname)
                        except Exception as e:
                            logger.warning(f"We got this exception: \n    {e}")
                            logger.warning("But we will try rebuild the Universe "
                                           "from GROM if using tpr did not work!")
                            op_obj = find_OP(system.content[key].mapping_dict,
                                             gro, xtcwhole, resname)

                    if 'openMM' in software or 'NAMD' in software:
                        op_obj = find_OP(system.content[key].mapping_dict,
                                         struc_fname, trj_fname, resname)

                    data = {}

                    with open(outfilename, 'w') as outfile:
                        outfile.write("Atom     Average OP     OP stem\n")

                        for i, op in enumerate(op_obj):
                            (op.avg, op.std, op.stem) = op.get_avg_std_stem_OP
                            outfile.write(f'{op.name} {str(op.avg)} {str(op.stem)}\n')

                            data[str(op.name)] = []
                            data[str(op.name)].append(op.get_avg_std_stem_OP)

                    with open(outfilename2, 'w') as f:
                        json.dump(data, f, cls=CompactJSONEncoder)

        print("Order parameters calculated and saved to ", path)

    except Exception as e:
        print('Calculation failed for ' + system['path'])
        print(str(e))
        print(traceback.format_exc())
        return RCODE_ERROR

    return RCODE_COMPUTED


def computeFF(  # noqa: N802 (API)
        system: System, logger: Logger, recompute: bool = False) -> int:
    logger.info("System title: " + system['SYSTEM'])
    logger.info("System path: " + system['path'])
    software = system['SOFTWARE']
    # download trajectory and gro files
    system_path = os.path.join(NMLDB_SIMU_PATH, system['path'])
    doi = system.get('DOI')
    skip_downloading: bool = (doi == 'localhost')
    if skip_downloading:
        print("NOTE: The system with 'localhost' DOI should be downloaded by the user.")

    if os.path.isfile(system_path + os.sep + "FormFactor.json") and not recompute:
        return RCODE_SKIPPED

    if system['TYPEOFSYSTEM'] == 'miscellaneous':
        return RCODE_SKIPPED

    print('Analyzing', system_path)

    try:
        if system['UNITEDATOM_DICT']:
            print('United atom simulation')
    except KeyError:
        pass

    try:
        if system['WARNINGS']['ORIENTATION']:
            print('Skipping due to ORIENTATION warning:',
                  system['WARNINGS']['ORIENTATION'])
            return RCODE_SKIPPED
    except (KeyError, TypeError):
        pass

    try:
        if system['WARNINGS']['PBC'] == 'hexagonal-box':
            print('Skipping due to PBC warning:', system['WARNINGS']['PBC'])
            return RCODE_SKIPPED
    except (KeyError, TypeError):
        pass

    try:
        if system['WARNINGS']['NOWATER']:
            print('Skipping because there is not water in the trajectory.')
            return RCODE_SKIPPED
    except (KeyError, TypeError):
        pass

    output_name = ""

    try:
        struc, top, trj = get_struc_top_traj_fnames(system)
        trj_name = os.path.join(system_path, trj)
        if struc is None:
            struc_name = None
        else:
            struc_name = os.path.join(system_path, struc)
        if top is None:
            tpr_name = None
        else:
            top_name = os.path.join(system_path, top)
    except Exception as e:
        logger.error("Error getting structure/topology/trajectory filenames.")
        logger.error(str(e))
        return RCODE_ERROR

    socket.setdefaulttimeout(15)

    try:
        if skip_downloading:
            if not os.path.isfile(trj_name):
                raise FileNotFoundError(
                    f"Trajectory should be downloaded [{trj_name}] by user")
        else:
            trj_url = resolve_download_file_url(system['DOI'], trj)
            if not os.path.isfile(trj_name):
                print('Downloading trajectory with the size of ',
                      system['TRAJECTORY_SIZE'], ' to ', system['path'])
                _ = urllib.request.urlretrieve(trj_url, trj_name)

        # make a function like this
        # TODO TPR should not be obligatory for GROMACS
        if 'gromacs' in software:
            tpr_name = top_name

            if skip_downloading:
                if not os.path.isfile(tpr_name):
                    raise FileNotFoundError(
                        f"TPR should be downloaded [{tpr_name}] by user")
            else:
                tpr_url = resolve_download_file_url(doi, top)
                if not os.path.isfile(tpr_name):
                    _ = urllib.request.urlretrieve(tpr_url, tpr_name)

        if 'openMM' in software or 'NAMD' in software:
            if skip_downloading:
                if (not os.path.isfile(struc_name)):
                    raise FileNotFoundError(
                        f"Structure file should be downloaded [{struc_name}] by user")
            else:
                pdb_url = resolve_download_file_url(doi, struc)
                if not os.path.isfile(struc_name):
                    _ = urllib.request.urlretrieve(pdb_url, struc_name)

        eq_time = float(system['TIMELEFTOUT'])*1000

        # FIND LAST CARBON OF SN-1 TAIL AND G3 CARBON
        for molecule in system['COMPOSITION']:
            if molecule in lipids_set:
                mapping = system.content[molecule].mapping_dict

                # TODO: rewrite via lipid dictionary!
                for nm in ["M_G3_M", "M_G13_M", "M_C32_M"]:
                    try:
                        g3_atom = mapping[nm]['ATOMNAME']
                        continue
                    except (KeyError, TypeError):
                        pass

                # TODO: rewrite via lipid dictionary
                if "M_G1C4_M" in mapping.keys():
                    for c_idx in range(4, 30):
                        atom = 'M_G1C' + str(c_idx) + '_M'
                        try:
                            last_atom = mapping[atom]['ATOMNAME']
                        except (KeyError, TypeError):
                            continue
                elif "M_G11C4_M" in mapping.keys():
                    for c_idx in range(4, 30):
                        atom = 'M_G11C' + str(c_idx) + '_M'
                        try:
                            last_atom = mapping[atom]['ATOMNAME']
                        except (KeyError, TypeError):
                            continue
                elif "M_CA4_M" in mapping.keys():
                    for c_idx in range(4, 30):
                        atom = 'M_CA' + str(c_idx) + '_M'
                        try:
                            last_atom = mapping[atom]['ATOMNAME']
                        except (KeyError, TypeError):
                            continue

        print(last_atom, g3_atom)

        # Center around one lipid tail CH3 to guarantee all lipids in the same box
        if 'gromacs' in system['SOFTWARE']:

            if (
                'WARNINGS' in system and
                system['WARNINGS'] is not None and
                'GROMACS_VERSION' in system['WARNINGS'] and
                system['WARNINGS']['GROMACS_VERSION'] == 'gromacs3'
            ):
                trjconv_cmd = 'trjconv'
                makendx_cmd = 'make_ndx'
            else:
                trjconv_cmd = 'gmx trjconv'
                makendx_cmd = 'gmx make_ndx'

            os.system('rm foo.ndx')
            os.system(f'echo "a {last_atom}\nq" | {makendx_cmd} -f {tpr_name} '
                      f'-o foo.ndx')
            os.system("tail -n1 foo.ndx | awk '{print $NF}'")
            os.system('echo "[ centralAtom ]" >> foo.ndx')
            os.system("tail -n2 foo.ndx | head -n1 |  awk '{print $NF}' >> foo.ndx")

            xtcwhole = os.path.join(system_path, 'whole.xtc')
            xtcfoo = os.path.join(system_path, 'foo2.xtc')
            xtccentered = os.path.join(system_path, 'centered.xtc')
            if (not os.path.isfile(xtccentered)):
                print("Make molecules whole in the trajectory")
                # if unitedAtom and system['TRAJECTORY_SIZE'] > 15000000000:
                #     print("United atom trajectry larger than 15 Gb. "
                #           "Using only every third frame to reduce memory usage.")
                #     os.system('echo System | gmx trjconv -f ' + trj_name + ' -s ' +
                #        tpr_name + ' -o ' + xtcwhole + ' -pbc mol -b ' + str(EQtime) +
                #        ' -skip 3')
                # else:
                if not os.path.isfile(xtcwhole):
                    os.system(f'echo System |  {trjconv_cmd} -f {trj_name} '
                              f'-s {tpr_name} -o {xtcwhole} -pbc mol -b {str(eq_time)}')

                if (not os.path.isfile(xtcfoo)):
                    os.system(f'echo "centralAtom\nSystem" |  {trjconv_cmd} -center'
                              f' -pbc mol -n foo.ndx -f {xtcwhole} -s {tpr_name}'
                              f' -o {xtcfoo}')

                os.system('rm foo.ndx')
                os.system('rm ' + xtcwhole)

            # Center around the center of mass of all the g_3 carbons
            # if (not os.path.isfile(xtccentered)):
                os.system(f'echo "a {g3_atom}\nq" | {makendx_cmd}'
                          f' -f {tpr_name} -o foo.ndx')
                os.system(f'echo "{g3_atom}\nSystem" |  {trjconv_cmd} -center'
                          f' -pbc mol -n foo.ndx -f {xtcfoo} -s {tpr_name}'
                          f' -o {xtccentered}')
                os.system('rm ' + xtcfoo)
        else:
            print("Centering for other than Gromacs may not work if there are"
                  " jumps over periodic boundary conditions in z-direction.")

        if not os.path.isfile(system_path + os.sep + "FormFactor.json"):
            try:
                if 'gromacs' in system['SOFTWARE']:
                    FormFactor(tpr_name, xtccentered, 200,
                               output_name, system)
                if 'openMM' in system['SOFTWARE'] or 'NAMD' in system['SOFTWARE']:
                    FormFactor(struc_name, trj_name, 200,
                               output_name, system)
            except ValueError as e:
                # Here it was expected to have allow_pickle-type errors.
                # but I suppose, we cannot simply ignore them because it means that the
                # code breaks at this place, and the running user should be informed
                # about it!
                raise e

    # finall catch
    except Exception as e:
        print('Calculation failed for ' + system['path'])
        print(str(e))
        print(traceback.format_exc())
        return RCODE_ERROR
    else:
        return RCODE_COMPUTED
