import os
import numpy as np
import yaml
import argparse

import ase.build
import ase.eos
import logging

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt




# read in the metal from command line
parser = argparse.ArgumentParser(description='Calculate the lattice constant of a metal.')
parser.add_argument('--metal', type=str, default='Pt', help='Metal to use for the calculation (default: Pt)')
parser.add_argument('--crystal_structure', type=str, default='fcc', help='Crystal structure (default: fcc)')
parser.add_argument('--initial_guess', type=float, default=3.5, help='Initial guess for the lattice constant (default: 3.5 Å)')
parser.add_argument('--facet', type=str, default='111', help='Facet if using the slab instead of bulk (default: 111)')

args = parser.parse_args()

metal = args.metal
crystal_structure = args.crystal_structure
initial_guess = args.initial_guess
facet = args.facet
plotting = True
results_dir = '../results'  # Directory to save results

use_slab = True  # Use slab for the calculation, otherwise use bulk
#use_slab = False  # Use slab for the calculation, otherwise use bulk

output_file = os.path.join(results_dir, 'bulk', f'{metal}_{crystal_structure}', f'{metal}_{crystal_structure}_lattice_constant.yaml')
if not os.path.exists(os.path.dirname(output_file)):
    os.makedirs(os.path.dirname(output_file))

# set up logging
logging.basicConfig(level=logging.INFO)


# Define the OCP calculator
# 'EquiformerV2-31M-S2EF-OC20-All+MD',
local_cache = os.environ['FAIRCHEM_LOCAL_CACHE']
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    'GemNet-OC-S2EFS-nsn-OC20+OC22',
    # 'EquiformerV2-31M-S2EF-OC20-All+MD',
    # local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
    # local_cache='/projects/westgroup/harris.se/surface_thermo/tmp/fairchem_checkpoints/'
    local_cache=local_cache
)

# checkpoint_path = '/home/moon/surface/tmp/fairchem_checkpoints/gnoc_oc22_oc20_all_s2ef.pt'
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(checkpoint_path=checkpoint_path, cpu=True, seed=400)


# Fit a polynomial to the energies
def too_much_noise(lattice_constants, energies, threshold=0.1):
    if np.max(lattice_constants) - np.min(lattice_constants) > 2.0:
        # The first round is too coarse, so don't trust the eos fit
        return True

    coeffs = np.polyfit(lattice_constants, energies, 2)
    # Calculate the fitted energies
    fitted_energies = np.polyval(coeffs, lattice_constants)

    # check if the fitted energies deviate too much from the original energies
    mean = np.mean(np.abs(energies - fitted_energies))

    range_energies = np.max(energies) - np.min(energies)

    # return variance > threshold * range_energies
    return mean > threshold * range_energies


def make_plot(lattice_constants, energies, title, filename, a0=None):
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants, energies, 'o-', label='Energies')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(title)
    if a0 is not None:
        plt.axvline(a0, color='g', linestyle='--', label='Final Lattice Constant')
        # annotate the final lattice constant
        plt.annotate(f'a0 = {a0:.4f} Å', xy=(a0, np.min(energies)), xytext=(a0 + np.std(lattice_constants), np.min(energies)),
                     arrowprops=dict(facecolor='black', shrink=0.05))
    plt.legend()
    plt.grid()
    plt.savefig(filename)
    plt.close()


def run_eos_analysis(description, a0_init, half_range):
    lattice_constants = np.linspace(a0_init - half_range, a0_init + half_range, 21)
    energies = np.zeros_like(lattice_constants)
    volumes = np.zeros_like(lattice_constants)
    for i in range(len(energies)):

        if use_slab:
            # Build a slab for the metal
            if crystal_structure == 'fcc':
                slab = ase.build.fcc111(metal, size=(3, 3, 4), vacuum=10.0, a=lattice_constants[i])
            elif crystal_structure == 'bcc':
                if facet == '110':
                    slab = ase.build.bcc110(metal, size=(3, 3, 4), vacuum=10.0, a=lattice_constants[i])
                elif facet == '100':
                    slab = ase.build.bcc100(metal, size=(3, 3, 4), vacuum=10.0, a=lattice_constants[i])
                else:
                    print('must specify a facet')       

                #if metal == 'Fe':
                #   # Chromium has a magnetic ground state, so we need to set the spin polarization
                #    slab.set_initial_magnetic_moments([2.0] * len(slab))
                #elif metal == 'Cr':
                #    slab.set_initial_magnetic_moments([2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                #                                       -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0,
                #                                       2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,])  # Set magnetic moments for bcc Cr

            else:
                raise ValueError(f"Unsupported crystal structure: {crystal_structure}")
            slab.calc = calc
            energies[i] = slab.get_potential_energy()
            equivalent_bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants[i], cubic=True)
            volumes[i] = equivalent_bulk.cell.volume

        else:
            bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants[i], cubic=True)

            if metal == 'Cr':
                # Chromium has a magnetic ground state, so we need to set the spin polarization
                pass
                # bulk.set_initial_magnetic_moments([2.0, -2.0])
            elif metal == 'Fe':
                # Iron also has a magnetic ground state, so we need to set the spin polarization
                pass
                bulk.set_initial_magnetic_moments([2.0, 2.0])

            bulk.calc = calc

            energies[i] = bulk.get_potential_energy()
            volumes[i] = bulk.cell.volume

    # Fit the equation of state
    eos = ase.eos.EquationOfState(volumes, energies, eos='sj')
    try:
        v0, e0, B = eos.fit()
        if use_slab:
            raise ValueError("Using slab, so EOS fit is not reliable, using minimum energy point as v0")
        if too_much_noise(lattice_constants, energies):
            raise ValueError("Too much noise in EOS, using minimum energy point as v0")
        
    except ValueError as e:
        logging.error(f"Error fitting EOS: {e}")
        # just take the minimum energy point as v0, or the previous guess if the half range is too small
        if half_range < 0.05:
            v0 = a0_init ** 3
        else:
            min_energy_index = np.argmin(energies)
            v0 = volumes[min_energy_index]
    a0 = np.float_power(v0, 1.0 / 3.0)
    logging.info(f"Lattice constant for {metal} in {crystal_structure} structure: {a0:.4f} Å")

    if plotting:
        make_plot(lattice_constants, energies,
                  title=f'Energy of State Analysis for {metal} ({crystal_structure}) - {description}',
                  filename=os.path.join(results_dir, 'bulk', f'{metal}_{crystal_structure}', f'{metal}_{crystal_structure}_{description}_analysis.png'),
                  a0=a0)
    return a0


levels_of_analysis = {
    # 'very_coarse': 1.500,
    'coarse': 0.500,
    'medium': 0.250,
    'medium_fine': 0.125,
    'fine': 0.050,
    'very_fine': 0.015,
    'ultra_fine': 0.005
}

level_results = {}

final_a0 = initial_guess  # Initial guess for the very coarse level
for level, half_range in levels_of_analysis.items():
    logging.info(f"Running {level} energy of state analysis...")
    final_a0 = run_eos_analysis(level, final_a0, half_range)
    level_results[level] = final_a0

logging.info(f'Final lattice constant for {metal} in {crystal_structure} structure: {final_a0:.12f} Å')


# Save final results to a yaml file
results = {
    'metal': metal,
    'crystal_structure': crystal_structure,
    'final_lattice_constant': float(final_a0),
}

with open(output_file, 'w') as f:
    yaml.dump(results, f, default_flow_style=False)
logging.info(f"Results saved to {output_file}")
