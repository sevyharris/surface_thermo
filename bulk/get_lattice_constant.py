import os
import numpy as np
import copy

import ase.build
import ase.optimize
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
plotting = True  # Set to True if you want to plot the results, False otherwise
final_lattice_constant = None  # This will hold the final lattice constant after analysis
# set up logging
logging.basicConfig(level=logging.INFO)


# Define the OCP calculator
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file('EquiformerV2-31M-S2EF-OC20-All+MD', local_cache='./tmp/fairchem_checkpoints/')
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(checkpoint_path=checkpoint_path, cpu=True, seed=400)


# Fit a polynomial to the energies
def too_much_deviation(lattice_constants, energies, threshold=0.1):
    coeffs = np.polyfit(lattice_constants, energies, 2)
    # Calculate the fitted energies
    fitted_energies = np.polyval(coeffs, lattice_constants)

    # check if the fitted energies deviate too much from the original energies
    mean = np.mean(np.abs(energies - fitted_energies))

    range_energies = np.max(energies) - np.min(energies)

    # return variance > threshold * range_energies
    return mean > threshold * range_energies



# Run very coarse energy of state analysis:
logging.info("Running very coarse energy of state analysis...")
lattice_constants_very_coarse = np.linspace(2.0, 5.0, 15)
energies_very_coarse = np.zeros_like(lattice_constants_very_coarse)
volumes_very_coarse = np.zeros_like(lattice_constants_very_coarse)
for i in range(len(energies_very_coarse)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_very_coarse[i], cubic=True)
    bulk.calc = calc
    energies_very_coarse[i] = bulk.get_potential_energy()
    volumes_very_coarse[i] = bulk.cell.volume

# don't trust eos fitting for the first round, just take the minimum energy point
min_energy_index = np.argmin(energies_very_coarse)
v0 = volumes_very_coarse[min_energy_index]
very_coarse_a0 = np.float_power(v0, 1.0 / 3.0)
logging.info(f"Very coarse lattice constant for {metal} in {crystal_structure} structure: {very_coarse_a0:.4f} Å")

if plotting:
    
    # Plot each step of the analysis separately
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_very_coarse, energies_very_coarse, 'o-', label='Very Coarse')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Very Coarse')
    plt.axvline(very_coarse_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {very_coarse_a0:.4f} Å', xy=(very_coarse_a0, np.min(energies_very_coarse)),
                xytext=(very_coarse_a0 + 0.05, np.min(energies_very_coarse) + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_very_coarse_analysis.png')



# Run coarse energy of state analysis:
logging.info("Running coarse energy of state analysis...")
lattice_constants_coarse = np.linspace(very_coarse_a0 - 0.5, very_coarse_a0 + 0.5, 15)
energies_coarse = np.zeros_like(lattice_constants_coarse)
volumes_coarse = np.zeros_like(lattice_constants_coarse)
for i in range(len(energies_coarse)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_coarse[i], cubic=True)
    bulk.calc = calc
    energies_coarse[i] = bulk.get_potential_energy()
    volumes_coarse[i] = bulk.cell.volume

coarse_eos = ase.eos.EquationOfState(volumes_coarse, energies_coarse, eos='sj')
try:
    v0, e0, B = coarse_eos.fit()
    if too_much_deviation(lattice_constants_coarse, energies_coarse):
        raise ValueError("Too much deviation in coarse EOS, using minimum energy point as v0")
except ValueError as e:
    logging.error(f"Error fitting coarse EOS: {e}")
    # just take the minimum energy point as v0
    min_energy_index = np.argmin(energies_coarse)
    v0 = volumes_coarse[min_energy_index]
coarse_a0 = np.float_power(v0, 1.0 / 3.0)
logging.info(f"Coarse lattice constant for {metal} in {crystal_structure} structure: {coarse_a0:.4f} Å")

if plotting:
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_coarse, energies_coarse, 'o-', label='Coarse')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Coarse')
    plt.axvline(coarse_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {coarse_a0:.4f} Å', xy=(coarse_a0, np.min(energies_coarse)),
                xytext=(coarse_a0 + 0.05, np.min(energies_coarse) + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_coarse_analysis.png')


# Run medium energy of state analysis:
logging.info("Running medium energy of state analysis...")
lattice_constants_medium = np.linspace(coarse_a0 - 0.2, coarse_a0 + 0.2, 21)
energies_medium = np.zeros_like(lattice_constants_medium)
volumes_medium = np.zeros_like(lattice_constants_medium)
for i in range(len(energies_medium)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_medium[i], cubic=True)
    bulk.calc = calc
    energies_medium[i] = bulk.get_potential_energy()
    volumes_medium[i] = bulk.cell.volume

# Fit the equation of state
medium_eos = ase.eos.EquationOfState(volumes_medium, energies_medium, eos='sj')
try:
    v0, e0, B = medium_eos.fit()
    if too_much_deviation(lattice_constants_medium, energies_medium):
        raise ValueError("Too much deviation in medium EOS, using minimum energy point as v0")
except ValueError as e:
    logging.error(f"Error fitting medium EOS: {e}")
    # just take the minimum energy point as v0
    min_energy_index = np.argmin(energies_medium)
    v0 = volumes_medium[min_energy_index]
medium_a0 = np.float_power(v0, 1.0 / 3.0)
logging.info(f"Medium lattice constant for {metal} in {crystal_structure} structure: {medium_a0:.4f} Å")
if plotting:
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_medium, energies_medium, 'o-', label='Medium')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Medium')
    plt.axvline(medium_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {medium_a0:.4f} Å', xy=(medium_a0, np.min(energies_medium)),
                xytext=(medium_a0 + 0.05, np.min(energies_medium) + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_medium_analysis.png')



# Run fine energy of state analysis:  # here the danger is being too zoomed in, so stop analysis if the deviation is too high
logging.info("Running fine energy of state analysis...")
lattice_constants_fine = np.linspace(medium_a0 - 0.05, medium_a0 + 0.05, 21)
energies_fine = np.zeros_like(lattice_constants_fine)
volumes_fine = np.zeros_like(lattice_constants_fine)
for i in range(len(energies_fine)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_fine[i], cubic=True)
    bulk.calc = calc
    energies_fine[i] = bulk.get_potential_energy()
    volumes_fine[i] = bulk.cell.volume

fine_eos = ase.eos.EquationOfState(volumes_fine, energies_fine, eos='sj')
try:
    v0, e0, B = fine_eos.fit()
    if too_much_deviation(lattice_constants_fine, energies_fine):
        # we already have a good estimate, so just take the minimum energy point as v0
        logging.info("Too much deviation in fine EOS, using medium EOS minimum")
        logging.info(f"Using medium lattice constant {medium_a0:.4f} Å as final estimate")
        final_lattice_constant = medium_a0
except ValueError as e:
    logging.error(f"Error fitting fine EOS: {e}")
    # just take the minimum energy point as v0
    min_energy_index = np.argmin(energies_fine)
    v0 = volumes_fine[min_energy_index]
fine_a0 = np.float_power(v0, 1.0 / 3.0)
logging.info(f"Fine lattice constant for {metal} in {crystal_structure} structure: {fine_a0:.4f} Å")

if plotting:
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_fine, energies_fine, 'o-', label='Fine')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Fine')
    plt.axvline(fine_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {fine_a0:.4f} Å', xy=(fine_a0, np.min(energies_fine)),
                xytext=(fine_a0 + 0.05, np.min(energies_fine)),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_fine_analysis.png')



# Run very fine energy of state analysis:
logging.info("Running very fine energy of state analysis...")
lattice_constants_very_fine = np.linspace(fine_a0 - 0.02, fine_a0 + 0.02, 21)
energies_very_fine = np.zeros_like(lattice_constants_very_fine)
volumes_very_fine = np.zeros_like(lattice_constants_very_fine)
for i in range(len(energies_very_fine)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_very_fine[i], cubic=True)
    bulk.calc = calc
    energies_very_fine[i] = bulk.get_potential_energy()
    volumes_very_fine[i] = bulk.cell.volume

very_fine_eos = ase.eos.EquationOfState(volumes_very_fine, energies_very_fine, eos='sj')
try:
    v0, e0, B = very_fine_eos.fit()
    if too_much_deviation(lattice_constants_very_fine, energies_very_fine):
        logging.info("Too much deviation in very fine EOS, using fine EOS minimum")
        logging.info(f"Using fine lattice constant {fine_a0:.4f} Å as final estimate")
        if final_lattice_constant is None:
            final_lattice_constant = fine_a0
except ValueError as e:
    logging.error(f"Error fitting very fine EOS: {e}")
    # just take the minimum energy point as v0
    min_energy_index = np.argmin(energies_very_fine)
    v0 = volumes_very_fine[min_energy_index]
very_fine_a0 = np.float_power(v0, 1.0 / 3.0)

if plotting:
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_very_fine, energies_very_fine, 'o-', label='Very Fine')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Very Fine')
    plt.axvline(very_fine_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {very_fine_a0:.4f} Å', xy=(very_fine_a0, np.min(energies_very_fine)),
                xytext=(very_fine_a0 + 0.05, np.min(energies_very_fine) + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_very_fine_analysis.png')




# Run ultra fine energy of state analysis:
logging.info("Running ultra fine energy of state analysis...")
lattice_constants_ultra_fine = np.linspace(fine_a0 - 0.005, fine_a0 + 0.005, 21)
energies_ultra_fine = np.zeros_like(lattice_constants_ultra_fine)
volumes_ultra_fine = np.zeros_like(lattice_constants_ultra_fine)
for i in range(len(energies_ultra_fine)):
    bulk = ase.build.bulk(metal, crystalstructure=crystal_structure, a=lattice_constants_ultra_fine[i], cubic=True)
    bulk.calc = calc
    energies_ultra_fine[i] = bulk.get_potential_energy()
    volumes_ultra_fine[i] = bulk.cell.volume

ultra_fine_eos = ase.eos.EquationOfState(volumes_ultra_fine, energies_ultra_fine, eos='sj')
try:
    v0, e0, B = ultra_fine_eos.fit()
    if too_much_deviation(lattice_constants_ultra_fine, energies_ultra_fine):
        logging.info("Too much deviation in ultra fine EOS, using very fine EOS minimum")
        logging.info(f"Using very fine lattice constant {very_fine_a0:.4f} Å as final estimate")
        if final_lattice_constant is None:
            final_lattice_constant = very_fine_a0
except ValueError as e:
    logging.error(f"Error fitting ultra fine EOS: {e}")
    # just take the minimum energy point as v0
    min_energy_index = np.argmin(energies_ultra_fine)
    v0 = volumes_ultra_fine[min_energy_index]
ultra_fine_a0 = np.float_power(v0, 1.0 / 3.0)

if plotting:
    plt.figure(figsize=(10, 6))
    plt.plot(lattice_constants_ultra_fine, energies_ultra_fine, 'o-', label='Ultra Fine')
    plt.xlabel('Lattice Constant (Å)')
    plt.ylabel('Energy (eV)')
    plt.title(f'Energy of State Analysis for {metal} ({crystal_structure}) - Ultra Fine')
    plt.axvline(ultra_fine_a0, color='r', linestyle='--', label='Final Lattice Constant')
    # add annotation of lattice constant
    plt.annotate(f'Lattice Constant: {ultra_fine_a0:.4f} Å', xy=(ultra_fine_a0, np.min(energies_ultra_fine)),
                xytext=(ultra_fine_a0 + 0.05, np.min(energies_ultra_fine) + 0.1),
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='black')
    plt.legend()
    plt.grid()
    plt.savefig(f'{metal}_{crystal_structure}_ultra_fine_analysis.png')



if final_lattice_constant is None:
    final_lattice_constant = ultra_fine_a0
logging.info(f'Final lattice constant for {metal} in {crystal_structure} structure: {final_lattice_constant:.12f} Å')