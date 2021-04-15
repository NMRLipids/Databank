import numpy as np
import MDAnalysis as mda
import MDAnalysis.transformations as trans
import argparse

center_trigger = False

allowed_max_distance = 2

u = mda.Universe(topol,traj)

nframes = u.trajectory.n_frames

# select 1 random frame #

random_frame = np.random.randint(0, nframes)

bonds = u.bonds

for i in bonds:

    atom1 = i.atoms.ids[0]
    atom2 = i.atoms.ids[1]

    distance = np.linalg.norm(u.trajectory[random_frame][atom1] - u.trajectory[random_frame][atom2])

    if distance > allowed_max_distance:

        print("unusually long bond between " + str(atom1) + " " + str(atom2) + " with length " + str(distance/10) + " nm")
        center_trigger = True

# the below selections must be adjusted when integrating into the AddData.py
membrane_string = 'resname POPC' #( resname POPC or resname DPPC .... )
not_membrane_string = 'not resname POPC'

if center_trigger:

    u = mda.Universe(topol, traj)

    membrane = u.select_atoms(membrane_string)
    not_membrane = u.select_atoms(not_membrane_string)
    everything = u.select_atoms(everything_string)

    transforms = [trans.unwrap(membrane),
              trans.center_in_box(membrane, wrap=True),
              trans.wrap(not_membrane)]

    u.trajectory.add_transformations(*transforms)


