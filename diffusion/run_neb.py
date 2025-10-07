# script to run NEB optimization after minima1 and minima2 have been computed
import os
import sys
import random
import glob

import ase.collections
import numpy as np
import copy
import yaml

import argparse

import ase.build
import ase.optimize
import ase.io.trajectory
import ase.visualize
import ase.vibrations
import logging

import ase.io
import ase.mep

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt

sys.path.append(os.environ['SURFACE_THERMO_DIR'])
# sys.path.append('/home/moon/surface/surface_thermo')
import util


working_dir = sys.argv[1]


def get_minimum_traj_file(glob_str, z_min=-1e5, z_max=1e5):
    min_traj_file = ''
    min_energy = 1e5

    traj_files = glob.glob(glob_str)
    for i in range(len(traj_files)):
        traj = ase.io.trajectory.Trajectory(traj_files[i])

        # exclude z's outside the range
        z = traj[-1].positions[-1, 2]
        if z < z_min or z > z_max:
            continue

        energy = traj[-1].calc.results['energy']
        if energy < min_energy:
            min_energy = energy
            min_traj_file = traj_files[i]
    return min_traj_file


local_cache = os.environ['FAIRCHEM_LOCAL_CACHE']
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    'GemNet-OC-S2EFS-nsn-OC20+OC22',
    # 'EquiformerV2-31M-S2EF-OC20-All+MD',
    # local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
    local_cache=local_cache
)
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(
    checkpoint_path=checkpoint_path,
    cpu=True,
    seed=400
)

# Read initial and final states:
min_traj_file1 = get_minimum_traj_file(os.path.join(working_dir, f'minima1_*.traj'))
min_traj_file2 = get_minimum_traj_file(os.path.join(working_dir, f'minima2_*.traj'))

initial = ase.io.read(min_traj_file1)
final = ase.io.read(min_traj_file2)

# Make a band consisting of N images:
N = 11
images = [initial]
images += [initial.copy() for i in range(N - 2)]
images += [final]
neb = ase.mep.NEB(images)
# Interpolate linearly the positions of the three middle images:
neb.interpolate()
# Set calculators:
for image in images[1:len(images) - 1]:
    image.calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(
        checkpoint_path=checkpoint_path,
        cpu=True,
        seed=400
    )


# Optimize:
trajectory = os.path.join(working_dir, f'A2B_{N}.traj')
optimizer = ase.optimize.BFGS(neb, trajectory=trajectory)
optimizer.run(fmax=0.05)
