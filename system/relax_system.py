import os
import sys

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

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt

sys.path.append(os.environ['SURFACE_THERMO_DIR'])
#sys.path.append('/home/moon/surface/surface_thermo')
import util



# parse arguments:
parser = argparse.ArgumentParser(description='Relax a system with an adsorbate on a metal slab.')
parser.add_argument('--metal', type=str, default='Pt', help='Metal to use for the slab (default: Pt)')
parser.add_argument('--crystal_structure', type=str, default='fcc', help='Crystal structure (default: fcc)')
parser.add_argument('--facet', type=str, default='111', help='Facet of the metal slab (default: 111)')
parser.add_argument('--adsorbate', type=str, default='H2', help='Adsorbate to use (default: H2)')
parser.add_argument('--site', type=str, default='hcp', help='Adsorption site (default: hcp)')
parser.add_argument('--plotting', action='store_true', help='Enable plotting of optimization energies')
parser.add_argument('--rotate', type=str, default='0', help='Rotate the adsorbate degrees (default: 0)')
args = parser.parse_args()


adsorbate_label = args.adsorbate
site = args.site
metal = args.metal
facet = args.facet
crystal_structure = args.crystal_structure
plotting = args.plotting
rotate = float(args.rotate)


# skip monatomic rotations:
if len(adsorbate_label) == 1 and rotate != 0:
    print('Skipping redundant monatomic rotations')
    exit(0)

# Initialize fairchem ocp calculator
local_cache = os.environ['FAIRCHEM_LOCAL_CACHE']
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    'GemNet-OC-S2EFS-nsn-OC20+OC22',
    # 'EquiformerV2-31M-S2EF-OC20-All+MD',
    #local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
    local_cache=local_cache
)
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(
    checkpoint_path=checkpoint_path,
    cpu=True,
    seed=400
)


# assert adsorbate_label in ase.collections.g2.keys()
# assert site in ['fcc', 'hcp', 'bridge', 'ontop'], f"Invalid site: {site}. Choose from 'fcc', 'hcp', 'bridge', or 'ontop'"


results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

plotting = False  # Set to True if you want to plot the results, False otherwise


# Start by loading the slab from the trajectory file if it exists
def load_slab_from_trajectory(metal, facet, crystal_structure, results_dir):
    """
    Load the slab from the trajectory file if it exists.
    """
    trajectory_file = os.path.join(results_dir, 'slab', f'{metal}_{crystal_structure}{facet}_slab.traj')
    if os.path.exists(trajectory_file):
        logging.info(f"Loading slab from existing trajectory file: {trajectory_file}")
        traj = ase.io.trajectory.Trajectory(trajectory_file)
        slab = traj[-1]  # Get the last frame from the trajectory
        return slab
    else:
        logging.info(f"No existing trajectory file found for {metal} {crystal_structure} {facet}.")
        return None


# Load the slab from the trajectory file if it exists
slab = load_slab_from_trajectory(metal, facet, crystal_structure, results_dir)

# fmax = 0.01
fmax = 0.05
if metal == 'Fe' and facet == '110':
    fmax = 0.1

MAXSTEP = 500  # Maximum number of optimization steps -- change this later to be larger
opt_complete = False  # Flag to check if the optimization is complete

# Try loading the system from the trajectory file if it exists
system_trajectory_file = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}_rot{rotate}.traj')
if not os.path.exists(os.path.dirname(system_trajectory_file)):
    os.makedirs(os.path.dirname(system_trajectory_file))
if os.path.exists(system_trajectory_file):
    logging.info(f"Loading system from existing trajectory file: {system_trajectory_file}")
    traj = ase.io.trajectory.Trajectory(system_trajectory_file, mode='r')
    system = traj[-1]  # Get the last frame from the trajectory
    # Check if the system is converged
    if util.atoms_converged(system, fmax=fmax):
        logging.info(f"System for {metal} {crystal_structure} {facet} + {adsorbate_label} {site} is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True
else:
    logging.info(f"No existing trajectory file found for {metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}. Starting from slab.")
    system = copy.deepcopy(slab)

    # Load the adsorbate geometry from trajectory file
    adsorbate_trajectory_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}.traj')
    if os.path.exists(adsorbate_trajectory_file):
        logging.info(f"Loading adsorbate from existing trajectory file: {adsorbate_trajectory_file}")
        adsorbate_traj = ase.io.trajectory.Trajectory(adsorbate_trajectory_file, mode='r')
        try:
            adsorbate = adsorbate_traj[-1]  # Get the last frame from the trajectory
        except IndexError:
            adsorbate = ase.build.molecule(adsorbate_label)
        if rotate != 0:
            logging.info(f"Rotating adsorbate {adsorbate_label} by {rotate} degrees")
            adsorbate.rotate(rotate, v='x', center='COM')
    else:
        adsorbate = ase.build.molecule(adsorbate_label)
        if rotate != 0:
            logging.info(f"Rotating adsorbate {adsorbate_label} by {rotate} degrees")
            adsorbate.rotate(rotate, v='x', center='COM')


    # calculate N energies and set at the lowest energy height
    logging.info(f"Finding optimal height for {adsorbate_label} on {metal} {crystal_structure} {facet} at {site}")
    heights = np.linspace(0.5, 3.0, 7)
    height_energies = np.zeros(len(heights))
    test_system = copy.deepcopy(slab)
    for i, height in enumerate(heights):
        ase.build.add_adsorbate(test_system, adsorbate, height, site)
        test_system.calc = calc
        height_energies[i] = test_system.get_potential_energy()
        # remove the adsorbate for the next iteration
        test_system = test_system[:len(slab)]
    best_height = heights[np.argmin(height_energies)]
    logging.info(f"Best height for {adsorbate_label} on {metal} {crystal_structure} {facet} at {site} is {best_height:.2f} Å")

    ase.build.add_adsorbate(system, adsorbate, best_height, site)


# system_constraint_indices = system.constraints[0].get_indices()
# # clear all constrants
# system.set_constraint()
# # fix everything but the adsorbate
# fixed_indices = [i for i in range(len(slab))]
# system.set_constraint(ase.constraints.FixAtoms(indices=fixed_indices))
# now, fix ALL the slab atoms in the system


system.calc = calc
logfile = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'ase_{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}_rot{rotate}.log')
if not opt_complete:
    logging.info(f'Running optimization for {metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}')
    opt = ase.optimize.BFGS(system, logfile=logfile, trajectory=system_trajectory_file, append_trajectory=True)
    opt.run(fmax=fmax, steps=MAXSTEP)


# make plot of optimization energies
if plotting:
    result_traj = ase.io.trajectory.Trajectory(system_trajectory_file, mode='r')
    energies = [frame.calc.results['energy'] for frame in result_traj]
    plt.figure()
    plt.plot(energies, label='Total Energy')
    plt.xlabel('Optimization Step')
    plt.ylabel('Energy (eV)')
    plt.title(f'Optimization Energy for {metal} {crystal_structure}{facet}_{adsorbate_label}_{site}')
    plt.legend()
    system_opt_plot_file = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'opt_energy_{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}.png')
    plt.savefig(system_opt_plot_file)
    plt.close()


# Run vibrational analysis
# Get the metal number and adsorbate indices
# This assumes the metal is the most abundant element in the system
# and the adsorbate is any other element.
logging.info(f'Running vibrational analysis for {metal} {crystal_structure}{facet}_{adsorbate_label}_{site}')
vals, counts = np.unique(system.get_atomic_numbers(), return_counts=True)
metal_number = vals[np.argmax(counts)]
adsorbate_indices = []
for i, atom in enumerate(system):
    if atom.number != metal_number:
        adsorbate_indices.append(atom.index)


vib = ase.vibrations.Vibrations(system, indices=adsorbate_indices)

# vib = ase.vibrations.Vibrations(system, name=f'{metal}{facet}_{adsorbate_label}_{site}', indices=adsorbate_indices)
vib.clean()  # Clean previous results
vib.run()
vib.summary()
freq = vib.get_frequencies()
logging.info(f'Vibrational frequencies for {metal} {crystal_structure}{facet}_{adsorbate_label}_{site}: {freq} (cm^-1)')
logging.info(f'ZPE for {metal} {crystal_structure}{facet}_{adsorbate_label}_{site}: {vib.get_zero_point_energy()} eV')

# Save result as a yaml file
result = {
    'frequencies': freq.tolist(),
    'zpe': float(vib.get_zero_point_energy()),
}
result_file = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}_rot{rotate}_vib.yaml')
with open(result_file, 'w') as f:
    yaml.dump(result, f, default_flow_style=False)


# Save a picture of the relaxed system
side_pic = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}_rot{rotate}_side.png')
top_pic = os.path.join(results_dir, 'system', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}', f'{metal}_{crystal_structure}{facet}_{adsorbate_label}_{site}_rot{rotate}_top.png')
ase.io.write(side_pic, system, rotation='-90x,0y,0z')
ase.io.write(top_pic, system, rotation='0x,0y,0z')
