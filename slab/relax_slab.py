import os
import numpy as np
import copy
import yaml

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


metal = 'Pt'  # Change this to the desired metal
crystal_structure = 'fcc'  # Change this to the desired crystal structure (e.g., 'fcc', 'bcc', 'hcp')
facet = '111'  # Change this to the desired facet (e.g., '111', '100', '110')
plotting = True  # Set to True if you want to plot the results, False otherwise
final_lattice_constant = None  # This will hold the final lattice constant after analysis
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
slab = ase.build.fcc111('Pt', size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)

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
    final_force = slab.calc.results['forces'] if 'forces' in slab.calc.results else np.inf
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
    opt = ase.optimize.BFGS(slab, trajectory=trajectory_file)
    opt.run(fmax=fmax, steps=1000)

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

