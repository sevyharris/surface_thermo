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

local_cache = os.environ['FAIRCHEM_LOCAL_CACHE']
sys.path.append(os.environ['SURFACE_THERMO_DIR'])
import util


#  python relax_custom_system.py --slabname=Cr2O3_z --adsorbate=H2 --rotate=0 --site=0


# Some Cr2O3 sites
sites_dict = {
    'Cr2O3_z': [
        (0.0, 0.0),
        (2.5018702100000003, 1.4444554392210054),
        (9.067513484506406e-16, 2.8889108784420108),
        (0.7547939021521936, 4.19625226621278),
        (0.7547939021521931, 1.5815694906712414)
    ],
    'Fe2O3_z': [
        (0, 0),
        (0.777535734, 1.3467314),
        (1.55507147, 0.0),
        (-0.777535734, 1.59307962),
        (0.388767867, 0.6733657),
        (1.166303602, 0.6733657),
        (-0.388767867, 0.79653981)
    ]
}

# parse arguments:
parser = argparse.ArgumentParser(description='Relax a system with an adsorbate on a metal slab.')
parser.add_argument('--slabname', type=str, default='Cr2O3_z', help='Name of slab to use for optimization (default: Cr2O3_z)')
parser.add_argument('--adsorbate', type=str, default='H2', help='Adsorbate to use (default: H2)')
parser.add_argument('--rotate', type=str, default='0', help='Rotate the adsorbate degrees (default: 0)')
parser.add_argument('--site', type=str, default=0, help='Site index to place the adsorbate on the slab (default: 0)')
args = parser.parse_args()


adsorbate_label = args.adsorbate
site = int(args.site)
slab_name = args.slabname
rotate = float(args.rotate)

results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)

# skip monatomic rotations:
if len(adsorbate_label) == 1 and rotate != 0:
    print('Skipping redundant monatomic rotations')
    exit(0)


# Start by loading the slab from the trajectory file if it exists
def load_slab_from_trajectory(slab_name, results_dir):
    """
    Load the slab from the trajectory file if it exists.
    """
    trajectory_file = os.path.join(results_dir, 'slab', f'{slab_name}_slab.traj')
    if os.path.exists(trajectory_file):
        logging.info(f"Loading slab from existing trajectory file: {trajectory_file}")
        traj = ase.io.trajectory.Trajectory(trajectory_file)
        slab = traj[-1]  # Get the last frame from the trajectory
        return slab
    else:
        logging.info(f"No existing trajectory file found for {slab_name}.")
        return None

# Load the slab from the trajectory file if it exists
slab = load_slab_from_trajectory(slab_name, results_dir)


# fmax = 0.01
fmax = 0.05
MAXSTEP = 1000  # Maximum number of optimization steps -- change this later to be larger
opt_complete = False  # Flag to check if the optimization is complete


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

adsorbate = None
# Try loading the system from the trajectory file if it exists
system_trajectory_file = os.path.join(results_dir, 'system', f'{slab_name}_{adsorbate_label}', f'{slab_name}_{adsorbate_label}_{site}_rot{rotate}.traj')
if not os.path.exists(os.path.dirname(system_trajectory_file)):
    os.makedirs(os.path.dirname(system_trajectory_file))
if os.path.exists(system_trajectory_file):
    logging.info(f"Loading system from existing trajectory file: {system_trajectory_file}")
    traj = ase.io.trajectory.Trajectory(system_trajectory_file, mode='r')
    system = traj[-1]  # Get the last frame from the trajectory
    # Check if the system is converged
    if util.atoms_converged(system, fmax=fmax):
        logging.info(f"System for {slab_name} + {adsorbate_label} {site} is already relaxed with forces below {fmax} eV/Å")
        opt_complete = True
else:
    logging.info(f"No existing trajectory file found for {slab_name}_{adsorbate_label}_{site}. Starting from slab.")
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

    # get the highest atom in the group
    indices = np.arange(len(system))
    heights = [atom.position[2] for atom in system]
    sorted_order = [x for _, x in sorted(zip(heights, indices))][::-1]

    # set tags to something so the calculator will work
    tags = system.get_tags()
    tags[sorted_order[0]] = 1
    system.set_tags(tags)
    system.info['top layer atom index'] = system[sorted_order[0]]

    # calculate N energies and set at the lowest energy height
    logging.info(f"Finding optimal height for {adsorbate_label} on {slab_name} at {site}")
    heights = np.linspace(0.1, 3.0, 11)
    height_energies = np.zeros(len(heights))
    test_system = copy.deepcopy(slab)
    for i, height in enumerate(heights):
        ase.build.add_adsorbate(test_system, adsorbate, height, position=sites_dict[slab_name][site])
        test_system.calc = calc
        height_energies[i] = test_system.get_potential_energy()
        # remove the adsorbate for the next iteration
        test_system = test_system[:len(slab)]
    best_height = heights[np.argmin(height_energies)]
    logging.info(f"Best height for {adsorbate_label} on {slab_name} at {site} is {best_height:.2f} Å")

    ase.build.add_adsorbate(system, adsorbate, best_height, position=sites_dict[slab_name][site])

# assume adsorbate gets added onto the end of the system
if adsorbate is None:
    adsorbate = ase.build.molecule(adsorbate_label)

adsorbate_indices = [i for i in np.arange(len(slab), len(slab) + len(adsorbate))]

logfile = os.path.join(results_dir, 'system', f'{slab_name}_{adsorbate_label}', f'ase_{slab_name}_{adsorbate_label}_{site}_rot{rotate}.log')
# Initialize fairchem ocp calculator

system.calc = calc
if not opt_complete:

    logging.info(f'Running optimization for {slab_name}_{adsorbate_label}_{site}')
    opt = ase.optimize.BFGS(system, logfile=logfile, trajectory=system_trajectory_file, append_trajectory=True)
    opt.run(fmax=fmax, steps=MAXSTEP)


# Run vibrational analysis
logging.info(f'Running vibrational analysis for {slab_name}_{adsorbate_label}_{site}')
vib = ase.vibrations.Vibrations(system, indices=adsorbate_indices)
vib.clean()  # Clean previous results
vib.run()
vib.summary()
freq = vib.get_frequencies()
logging.info(f'Vibrational frequencies for {slab_name}_{adsorbate_label}_{site}: {freq} (cm^-1)')
logging.info(f'ZPE for {slab_name}_{adsorbate_label}_{site}: {vib.get_zero_point_energy()} eV')

# Save result as a yaml file
result = {
    'frequencies': freq.tolist(),
    'zpe': float(vib.get_zero_point_energy()),
}
result_file = os.path.join(results_dir, 'system', f'{slab_name}_{adsorbate_label}', f'{slab_name}_{adsorbate_label}_{site}_rot{rotate}_vib.yaml')
with open(result_file, 'w') as f:
    yaml.dump(result, f, default_flow_style=False)


# fails to draw if occupancy is in system.info
if 'occupancy' in system.info:
    logging.info(f"Removing occupancy information from system info")
    system.info.pop('occupancy')

# Save a picture of the relaxed system
side_pic = os.path.join(results_dir, 'system', f'{slab_name}_{adsorbate_label}', f'{slab_name}_{adsorbate_label}_{site}_rot{rotate}_side.png')
top_pic = os.path.join(results_dir, 'system', f'{slab_name}_{adsorbate_label}', f'{slab_name}_{adsorbate_label}_{site}_rot{rotate}_top.png')
ase.io.write(side_pic, system, rotation='-90x,0y,0z')
ase.io.write(top_pic, system, rotation='0x,0y,0z')
