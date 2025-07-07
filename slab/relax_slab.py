import os
import numpy as np
import copy
import yaml
import argparse

import ase.build
import ase.optimize
import ase.io.trajectory
import ase.visualize
import ase.eos
import ase.constraints
import ase.vibrations
import logging

import fairchem.core.models.model_registry
import fairchem.core.models.equiformer_v2.trainers.forces_trainer
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt


# read in the metal, crystal structure, facet, and plotting options from command line
parser = argparse.ArgumentParser(description='Relax a metal slab with a specific facet.')
parser.add_argument('--metal', type=str, default='Pt', help='Metal to use for the slab (default: Pt)')
parser.add_argument('--crystal_structure', type=str, default='fcc', help='Crystal structure (default: fcc)')
parser.add_argument('--facet', type=str, default='111', help='Facet of the metal slab (default: 111)')


args = parser.parse_args()

metal = args.metal
crystal_structure = args.crystal_structure
facet = args.facet
plotting = True
lattice_constant = None  # This will hold the final lattice constant after analysis
results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

plotting = True  # Set to True if you want to plot the results, False otherwise


# start by loading the lattice constant result from the bulk yaml file
bulk_yaml_file = os.path.join(results_dir, 'bulk', f'{metal}_{crystal_structure}_lattice_constant.yaml')
with open(bulk_yaml_file, 'r') as f:
    data = yaml.load(f, Loader=yaml.Loader)
    lattice_constant = data.get('final_lattice_constant', None)
logging.info(f"Loaded lattice constant for {metal} {crystal_structure}: {lattice_constant} Å")

fmax = 0.01
vacuum = 7.5
if crystal_structure == 'fcc':
    if facet == '111':
        slab = ase.build.fcc111(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    elif facet == '100':
        slab = ase.build.fcc100(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    elif facet == '110':
        slab = ase.build.fcc110(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    else:
        raise ValueError(f"Invalid facet: {facet}. Choose from '111', '100', or '110'.")
elif crystal_structure == 'bcc':
    if facet == '110':
        slab = ase.build.bcc110(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    elif facet == '100':
        slab = ase.build.bcc100(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    elif facet == '211':
        slab = ase.build.bcc211(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
    else:
        raise ValueError(f"Invalid facet: {facet}. Choose from '110', '100', or '211'.")


# Fix the bottom two layers
bottom_layer = []
second_layer = []
fixed_indices = []
z_values = list(set([pos[2] for pos in slab.get_positions()]))
z_values.sort()

for i, pos in enumerate(slab.get_positions()):
    if pos[2] == z_values[0]:
       bottom_layer.append(slab[i].index)
    if pos[2] == z_values[1]:
       second_layer.append(slab[i].index)
fixed_indices = bottom_layer + second_layer
fix_bottom_layers = ase.constraints.FixAtoms(indices=fixed_indices)
slab.set_constraint(fix_bottom_layers)

if metal == 'Cr':
    # For Cr, we need to set the magnetic moments as antiferromagnetic
    z_levels = sorted(list(set(slab.positions[:, 2])))
    layer_index_sets = []
    for j, z_level in enumerate(z_levels):
        level_set = []
        for i in range(len(slab)):
            if slab.positions[i, 2] == z_level:
                level_set.append(i)
        layer_index_sets.append(level_set)

    magmoms = np.zeros_like(slab.get_initial_magnetic_moments())

    STARTING_MAG = 2.0
    for i, layer in enumerate(layer_index_sets):
        for atom_index in layer:
            if i % 2 == 0:
                magmoms[atom_index] = STARTING_MAG
            else:
                magmoms[atom_index] = -STARTING_MAG
    slab.set_initial_magnetic_moments(magmoms)

elif metal == 'Fe':
    # For Fe, we need to set the magnetic moment
    slab.set_initial_magnetic_moments([2.0] * len(slab))


checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file('EquiformerV2-31M-S2EF-OC20-All+MD', local_cache='/home/moon/surface/tmp/fairchem_checkpoints/')
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(checkpoint_path=checkpoint_path, cpu=True, seed=400)

logging.info(f'Running slab relaxation for {metal} {facet} with lattice constant {lattice_constant} Å')
slab.calc = calc
trajectory_file = os.path.join(results_dir, 'slab', f'{metal}{facet}_slab.traj')
opt_complete = False
if not os.path.exists(os.path.dirname(trajectory_file)):
    os.makedirs(os.path.dirname(trajectory_file))
# try loading the slab from the trajectory file if it exists
if os.path.exists(trajectory_file):
    logging.info(f"Loading slab from existing trajectory file: {trajectory_file}")
    traj = ase.io.trajectory.Trajectory(trajectory_file)
    slab = traj[-1]  # Get the last frame from the trajectory
    # only run the optimization if the slab is not already relaxed
    try:
        final_force = slab.calc.results['forces']
    except KeyError:
        final_force = np.zeros_like(slab.get_positions()) + np.inf

    # only count the forces on the atoms that are not fixed
    final_force = [final_force[i, :] if i not in fixed_indices else np.zeros(3) for i in range(final_force.shape[0])]

    logging.info(f"Final forces on slab: {np.max(np.linalg.norm(final_force, axis=1))}")
    if np.max(np.linalg.norm(final_force, axis=1)) < fmax:
    
        logging.info(f"Slab for {metal} {facet} is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True
else:
    logging.info(f"Creating new slab for {metal} {facet} with lattice constant {lattice_constant} Å")

# If the slab is not relaxed, run the optimization
if not opt_complete:
    slab.calc = calc
    opt = ase.optimize.BFGS(slab, trajectory=trajectory_file, append_trajectory=True)
    opt.run(fmax=fmax, steps=10000)

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
    plt.title(f'{metal} {facet} Slab Relaxation Energies')
    plt.xlabel('Step')
    plt.ylabel('Energy (eV)')
    plt.grid()
    plt.savefig(os.path.join(results_dir, 'slab', f'{metal}{facet}_slab_relaxation_energies.png'))


# save a picture of the relaxed slab
ase.io.write(os.path.join(results_dir, 'slab', f'{metal}{facet}_slab_top.png'), slab, rotation="0x,0y,0z")
ase.io.write(os.path.join(results_dir, 'slab', f'{metal}{facet}_slab_side.png'), slab, rotation="-90x,0y,0z")
