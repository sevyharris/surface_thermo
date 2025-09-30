# script to relax a Fe2O3 slab using the z direction as the surface normal

import os
import sys
import numpy as np
import yaml
import copy

import ase.build
import ase.optimize

import ase.io.cif
import ase.io.trajectory
import ase.constraints
import logging

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt

local_cache = os.environ['FAIRCHEM_LOCAL_CACHE']
sys.path.append(os.environ['SURFACE_THERMO_DIR'])
import util

slab_name = 'Fe2O3_z'
plotting = True
results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

# Start by reading in .cif file for Fe2O3
cif_file = 'Fe2O3.cif'
slab = ase.io.cif.read_cif('Fe2O3.cif')
logging.info(f"Loaded Fe2O3 from CIF file: {cif_file}")


# assign the layers in the slab, this is only approximate for custom structures
layers = util.enumerate_layers(slab)

# set tags to something so the calculator will work
slab.set_tags(layers)

# get the highest atom in the group and make sure it has a tag of 1 (top layer atom)
indices = np.arange(len(slab))
heights = [atom.position[2] for atom in slab]
sorted_order = [x for _, x in sorted(zip(heights, indices))][::-1]
assert slab.get_tags()[sorted_order[0]] == 1
slab.info['top layer atom index'] = sorted_order[0]


# stack unit cells in x, y, and z directions to make a slab
slab = ase.build.stack(slab, slab, axis=0)
slab = ase.build.stack(slab, slab, axis=1)
height = slab.get_cell()[2][2]
slab = ase.build.stack(slab, slab, axis=2)

# constrain all atoms below height to fix bulk atom positions
fixed_indices = []
for i, pos in enumerate(slab.get_positions()):
    if pos[2] < height:
        fixed_indices.append(i)
fix_bottom_layers = ase.constraints.FixAtoms(indices=fixed_indices)
slab.set_constraint(fix_bottom_layers)

# add vacuum
vacuum = 20
cell = slab.get_cell()
cell[2][2] += vacuum
slab.set_cell(cell)

# TODO see if changing initial magnetic moments changes anything

# keep track of original positions to track displacements
original = copy.deepcopy(slab)


# Run optimization
# fmax = 0.01
fmax = 0.05
max_steps = 10000

trajectory_file = os.path.join(results_dir, 'slab', f'Fe2O3_z_slab.traj')
opt_complete = False
if not os.path.exists(os.path.dirname(trajectory_file)):
    os.makedirs(os.path.dirname(trajectory_file))
# try loading the slab from the trajectory file if it exists
if os.path.exists(trajectory_file):
    logging.info(f"Loading slab from existing trajectory file: {trajectory_file}")
    traj = ase.io.trajectory.Trajectory(trajectory_file)
    slab = traj[-1]  # Get the last frame from the trajectory
    # only run the optimization if the slab is not already relaxed

    if util.atoms_converged(slab, fmax=fmax):
        logging.info(f"Slab for Fe2O3 z facet is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True
else:
    logging.info(f"Creating new slab for Fe2O3 z facet")

# If the slab is not relaxed, run the optimization
if not opt_complete:
    # only load calculator if we need to run the optimization
    logging.info(f"Loading OC20 calculator for slab relaxation")
    checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
        # 'EquiformerV2-31M-S2EF-OC20-All+MD',
        'GemNet-OC-S2EFS-OC20+OC22',
        # local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
        local_cache=local_cache
    )
    calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(checkpoint_path=checkpoint_path, cpu=True, seed=400)
    slab.calc = calc
    opt = ase.optimize.BFGS(slab, trajectory=trajectory_file, append_trajectory=True)
    opt.run(fmax=fmax, steps=max_steps)

if plotting:
    # get energies from the trajectory
    energies = []
    traj = ase.io.trajectory.Trajectory(trajectory_file)
    for frame in traj:
        energies.append(frame.calc.results['energy'])

    energies = np.array(energies)
    # plot the energies
    plt.figure(figsize=(8, 6))
    plt.plot(energies, marker='o', linestyle='-', color='b')
    plt.title(f'{slab_name} - Slab Relaxation Energies')
    plt.xlabel('Step')
    plt.ylabel('Energy (eV)')
    plt.grid()
    plt.savefig(os.path.join(results_dir, 'slab', f'{slab_name}_slab_relaxation_energies.png'))


logging.info(f"Max displacement of atoms during relaxation: {np.max(np.linalg.norm(slab.get_positions() - original.get_positions(), axis=1))} Å")
logging.info(f"Median displacement of atoms during relaxation: {np.median(np.linalg.norm(slab.get_positions() - original.get_positions(), axis=1))} Å")

# fails to draw if occupancy is in slab.info
if 'occupancy' in slab.info:
    logging.info(f"Removing occupancy information from slab info")
    slab.info.pop('occupancy')


# save a picture of the relaxed slab
ase.io.write(os.path.join(results_dir, 'slab', f'{slab_name}_slab_top.png'), slab, rotation="0x,0y,0z")
ase.io.write(os.path.join(results_dir, 'slab', f'{slab_name}_slab_side.png'), slab, rotation="-90x,0y,0z")
