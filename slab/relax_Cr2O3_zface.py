# script to relax a Cr2O3 slab using the z direction as the surface normal

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

sys.path.append('/home/moon/surface/surface_thermo')
import util

slab_name = 'Cr2O3_z'
plotting = True
results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

# Start by reading in .cif file for Cr2O3
cif_file = 'Cr2O3.cif'
slab = ase.io.cif.read_cif('Cr2O3.cif')
logging.info(f"Loaded Cr2O3 from CIF fiele: {cif_file}")

# set tags to something so the calculator will work
tags = slab.get_tags()
tags[1] = 1
slab.set_tags(tags)

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
fmax = 0.01
max_steps = 10000

trajectory_file = os.path.join(results_dir, 'slab', f'Cr2O3_z_slab.traj')
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
        logging.info(f"Slab for Cr2O3 z facet is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True
else:
    logging.info(f"Creating new slab for Cr2O3 z facet")

# If the slab is not relaxed, run the optimization
if not opt_complete:
    # only load calculator if we need to run the optimization
    logging.info(f"Loading OC20 calculator for slab relaxation")
    checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
        # 'EquiformerV2-31M-S2EF-OC20-All+MD',
        'GemNet-OC-S2EFS-OC20+OC22',
        local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
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


# save a picture of the relaxed slab
ase.io.write(os.path.join(results_dir, 'slab', f'{slab_name}_slab_top.png'), slab, rotation="0x,0y,0z")
ase.io.write(os.path.join(results_dir, 'slab', f'{slab_name}_slab_side.png'), slab, rotation="-90x,0y,0z")
