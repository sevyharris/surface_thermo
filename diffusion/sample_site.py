# script to sample individual site using sobol sequence

import os
import sys
import random
import glob

import scipy.stats.qmc

import ase.collections
import numpy as np
import copy
import yaml

import argparse

import ase.build
import ase.optimize
import ase.constraints
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

working_dir = '/scratch/harris.se/guassian_scratch/interstitials/Fe_fcc111_O/sobol_512'
os.makedirs(working_dir, exist_ok=True)

sobol_index = int(sys.argv[1])
print(f'sobol index: {sobol_index}')

slab_name = 'Fe_fcc111'
element = 'O'

results_dir = os.path.join(os.environ['SURFACE_THERMO_DIR'], 'results')
slab_file = os.path.join(results_dir, 'slab', f'{slab_name}_slab.traj')


traj = ase.io.trajectory.Trajectory(slab_file)
slab = traj[-1]

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


np.random.seed(400)
sobol = scipy.stats.qmc.Sobol(3)
N = 2 ** 9
normalized_coords = sobol.random(N)


def normalized2xyz(x_norm, y_norm, z_norm):
    # 13 - 14.8 is layer 2
    # 14.5 - 16 is layer 1
    # 16+ is surface
    minx = 0.0
    miny = 0.0
    minz = 13.0
    maxx = 4.0
    maxy = 4.0
    maxz = 16.0

    x = x_norm * (maxx - minx) + minx
    y = y_norm * (maxy - miny) + miny
    z = z_norm * (maxz - minz) + minz
    return x, y, z


test_system = copy.deepcopy(slab)
ase.build.add_adsorbate(test_system, element, height=0, position=(0, 0))
test_system.positions[-1, :] = normalized2xyz(*normalized_coords[sobol_index, :])
test_system.calc = calc

system_trajectory_file = os.path.join(working_dir, f'sobol_{sobol_index:04}.traj')


fmax = 0.05
MAXSTEP = 100
opt = ase.optimize.BFGS(test_system, trajectory=system_trajectory_file)
opt.run(fmax=fmax, steps=MAXSTEP)
