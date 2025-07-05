import os
import numpy as np
import copy
import yaml

import ase.build
import ase.optimize
import ase.io.trajectory
import ase.visualize
import ase.vibrations
import logging

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt


# Initialize fairchem ocp calculator
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    'EquiformerV2-31M-S2EF-OC20-All+MD',
    local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
)
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(
    checkpoint_path=checkpoint_path,
    cpu=True,
    seed=400
)

adsorbate_label = 'H2'  # Change this to the desired adsorbate (e.g., 'H2', 'O2', 'CO', etc.)
site = 'hcp'  # Change this to the desired adsorption site ('hcp', 'ontop', 'fcc', 'bridge')
metal = 'Pt'  # Change this to the desired metal
crystal_structure = 'fcc'  # Change this to the desired crystal structure (e.g., 'fcc', 'bcc', 'hcp')
facet = '111'  # Change this to the desired facet (e.g., '111', '100', '110')
plotting = True  # Set to True if you want to plot the results, False otherwise

results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

plotting = True  # Set to True if you want to plot the results, False otherwise


# Start by loading the slab from the trajectory file if it exists
def load_slab_from_trajectory(metal, facet, results_dir):
    """
    Load the slab from the trajectory file if it exists.
    """
    trajectory_file = os.path.join(results_dir, 'slab', f'{metal}{facet}_slab.traj')
    if os.path.exists(trajectory_file):
        logging.info(f"Loading slab from existing trajectory file: {trajectory_file}")
        traj = ase.io.trajectory.Trajectory(trajectory_file)
        slab = traj[-1]  # Get the last frame from the trajectory
        return slab
    else:
        logging.info(f"No existing trajectory file found for {metal} {facet}.")
        return None

# Load the slab from the trajectory file if it exists
slab = load_slab_from_trajectory(metal, facet, results_dir)


fmax = 0.01
MAXSTEP = 500  # Maximum number of optimization steps -- change this later to be larger
opt_complete = False  # Flag to check if the optimization is complete

# Try loading the system from the trajectory file if it exists
trajectory_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'{metal}{facet}_{adsorbate_label}_{site}.traj')
if not os.path.exists(os.path.dirname(trajectory_file)):
    os.makedirs(os.path.dirname(trajectory_file))
if os.path.exists(trajectory_file):
    logging.info(f"Loading system from existing trajectory file: {trajectory_file}")
    traj = ase.io.trajectory.Trajectory(trajectory_file)
    system = traj[-1]  # Get the last frame from the trajectory

    # check if the slab is already relaxed
    fixed_indices = system.constraints[0].get_indices()
    final_force = slab.calc.results['forces'] if 'forces' in slab.calc.results else np.inf
    # only count the forces on the atoms that are not fixed
    final_force = [final_force[i, :] if i not in fixed_indices else np.zeros(3) for i in range(final_force.shape[0])]

    logging.info(f"Final forces on slab: {np.max(np.linalg.norm(final_force, axis=1))}")
    if np.max(np.linalg.norm(final_force, axis=1)) < fmax:
        logging.info(f"System for {metal} {facet} + {adsorbate_label} {site} is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True

else:
    logging.info(f"No existing trajectory file found for {metal}{facet}_{adsorbate_label}_{site}. Starting from slab.")
    system = copy.deepcopy(slab)

    # Load the adsorbate geometry from trajectory file
    adsorbate_trajectory_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}.traj')
    logging.info(f"Loading adsorbate from existing trajectory file: {adsorbate_trajectory_file}")
    adsorbate_traj = ase.io.trajectory.Trajectory(adsorbate_trajectory_file)
    adsorbate = adsorbate_traj[-1]  # Get the last frame from the trajectory
    ase.build.add_adsorbate(system, adsorbate, 2.0, site)

    system.calc = calc
    logfile = 'ase.log'
    trajectory_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'{metal}{facet}_{adsorbate_label}_{site}.traj')
    logfile = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'ase_{metal}{facet}_{adsorbate_label}_{site}.log')


if not opt_complete:
    logging.info(f'Running optimization for {metal}{facet}_{adsorbate_label}_{site}')
    opt = ase.optimize.BFGS(system, logfile=logfile, trajectory=trajectory_file)
    opt.run(fmax=fmax, steps=MAXSTEP)


# make plot of optimization energies
if plotting:
    result_traj = ase.io.trajectory.Trajectory(trajectory_file)
    energies = [frame.calc.results['energy'] for frame in result_traj]
    plt.figure()
    plt.plot(energies, label='Total Energy')
    plt.xlabel('Optimization Step')
    plt.ylabel('Energy (eV)')
    plt.title(f'Optimization Energy for {metal}{facet}_{adsorbate_label}_{site}')
    plt.legend()
    plt.savefig(os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'opt_energy_{metal}{facet}_{adsorbate_label}_{site}.png'))
    plt.close()


# Run vibrational analysis
# Get the metal number and adsorbate indices
# This assumes the metal is the most abundant element in the system
# and the adsorbate is any other element.
logging.info(f'Running vibrational analysis for {metal}{facet}_{adsorbate_label}_{site}')
vals, counts = np.unique(system.get_atomic_numbers(), return_counts=True)
metal_number = vals[np.argmax(counts)]
adsorbate_indices = []
for i, atom in enumerate(system):
    if atom.number != metal_number:
        adsorbate_indices.append(atom.index)


vib = ase.vibrations.Vibrations(system, name=f'{metal}{facet}_{adsorbate_label}_{site}', indices=adsorbate_indices)
vib.run()
vib.summary()
freq = vib.get_frequencies()
logging.info(f'Vibrational frequencies for {metal}{facet}_{adsorbate_label}_{site}: {freq} (cm^-1)')
logging.info(f'ZPE for {metal}{facet}_{adsorbate_label}_{site}: {vib.get_zero_point_energy()} eV')

# Save result as a yaml file
result = {
    'frequencies': freq.tolist(),
    'zpe': float(vib.get_zero_point_energy()),
}
result_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'{metal}{facet}_{adsorbate_label}_{site}_vib.yaml')
with open(result_file, 'w') as f:
    yaml.dump(result, f, default_flow_style=False)


# Save a picture of the relaxed system
side_pic = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'{metal}{facet}_{adsorbate_label}_{site}_side.png')
top_pic = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', f'{metal}{facet}_{adsorbate_label}_{site}_top.png')
ase.io.write(side_pic, system, rotation='-90x,0y,0z')
ase.io.write(top_pic, system, rotation='0x,0y,0z')
