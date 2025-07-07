import os
import matplotlib.pyplot as plt
import numpy as np
import rmgpy.constants
import adsorbate_thermo
# import rmgpy.species
import sys
sys.path.append('/home/moon/surface/surface_thermo')
import util
import yaml
import ase.io.trajectory
import logging


logging.basicConfig(level=logging.INFO)

# script to convert binding energies to NASA polynomials
metal = 'Pt'
facet = '111'

if facet == '111':
    sites = ['ontop', 'hcp', 'bridge']  # TODO, finish converging bridge for Pt111
    # sites = ['ontop', 'hcp']
    # sites = ['hcp']
elif facet == '110':
    sites = ['ontop', 'hollow']


site = 'ontop'
adsorbate_label = 'H'
results_dir = f'../results'

eV_to_kJ = 96.485332

# -------------- Gather the binding energies ----------------

# Pick the site with the minimum energy
site_energies = np.zeros(len(sites))
for i, site in enumerate(sites):
    system_traj_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', 
                                    f'{metal}{facet}_{adsorbate_label}_{site}.traj')
    system_trajectory = ase.io.trajectory.Trajectory(system_traj_file)
    system = system_trajectory[-1]
    assert util.atoms_converged(system), f'System {metal}{facet}_{adsorbate_label}_{site} is not converged'
    system_energy = system.calc.results['energy']  # in eV
    site_energies[i] = system_energy

# # pick the hcp site if the bridge site is close (bridge likely fell into hcp)
# if facet == '111' and np.abs(site_energies[1] - site_energies[2]) < 0.1:
#     print(f'Bridge site is close to HCP site for {metal}{facet}, using HCP site instead.')
#     site_energies[1] = site_energies[2]

# Print the minimum energy site
min_energy_site = sites[np.argmin(site_energies)]
print(f'Minimum energy site for {adsorbate_label} on {metal}{facet}: {min_energy_site}')

# reload the system with the minimum energy site
system_traj_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', 
                                f'{metal}{facet}_{adsorbate_label}_{min_energy_site}.traj')
system_trajectory = ase.io.trajectory.Trajectory(system_traj_file)
system = system_trajectory[-1]
site = min_energy_site  # update the site to the minimum energy site
system_energy = system.calc.results['energy']  # in eV


# get the slab energy
slab_traj_file = os.path.join(results_dir, 'slab', f'{metal}{facet}_slab.traj')
slab_trajectory = ase.io.trajectory.Trajectory(slab_traj_file)
slab = slab_trajectory[-1]
slab_energy = slab.calc.results['energy']  # in eV

# Gather the gas energy
gas_traj_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}.traj')
gas_trajectory = ase.io.trajectory.Trajectory(gas_traj_file)
gas = gas_trajectory[-1]
gas_energy = gas.calc.results['energy']  # in eV



# Calculate the binding energy
binding_energy = system_energy - (slab_energy + gas_energy)


print(f'Gas energy for {adsorbate_label}: {gas_energy:.4f} eV')
print(f'Slab energy for {metal}{facet}: {slab_energy:.4f} eV')
print(f'System energy for {metal}{facet}_{adsorbate_label}_{site}: {system_energy:.4f} eV')
print()
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy:.4f} eV')
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy * eV_to_kJ:.4f} kJ/mol')

print((binding_energy * eV_to_kJ + 239.2) / -239.2)



# -------------- Convert to NASA polynomials ----------------
# get molecular weight for gas
molecular_weight = np.sum(gas.get_masses())  # in amu
print(f'Molecular weight of {adsorbate_label}: {molecular_weight:.4f} amu')

# read in the frequencies from the system vibrational yaml file
vibrational_yaml_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', 
                                     f'{metal}{facet}_{adsorbate_label}_{site}_vib.yaml')
# read yaml
with open(vibrational_yaml_file, 'r') as f:
    vib_data = yaml.load(f, Loader=yaml.FullLoader)

# get magnitude to avoid imaginary components
frequencies = np.abs(vib_data['frequencies'])
print(f'Vibrational frequencies for {metal}{facet}_{adsorbate_label}_{site}: {frequencies} (cm^-1)')

# build the composition dictionary
composition = {}
for atom in gas:
    if atom.symbol not in composition:
        composition[atom.symbol] = 0
    composition[atom.symbol] += 1
print(f'Composition of {adsorbate_label}: {composition}')

#------------------- Grab the reference gases -------------------
# Reference gas for H2
H2_traj_file = os.path.join(results_dir, 'gas', f'H2.traj')
H2_trajectory = ase.io.trajectory.Trajectory(H2_traj_file)
H2_atoms = H2_trajectory[-1]
H2_energy = H2_atoms.calc.results['energy']  # in eV

# Reference gas for H2O
H2O_traj_file = os.path.join(results_dir, 'gas', f'H2O.traj')
H2O_trajectory = ase.io.trajectory.Trajectory(H2O_traj_file)
H2O_atoms = H2O_trajectory[-1]
H2O_energy = H2O_atoms.calc.results['energy']  # in eV

# Reference gas for CH4
CH4_traj_file = os.path.join(results_dir, 'gas', f'CH4.traj')
CH4_trajectory = ase.io.trajectory.Trajectory(CH4_traj_file)
CH4_atoms = CH4_trajectory[-1]
CH4_energy = CH4_atoms.calc.results['energy']  # in eV

# Reference gas for NH3
NH3_traj_file = os.path.join(results_dir, 'gas', f'NH3.traj')
NH3_trajectory = ase.io.trajectory.Trajectory(NH3_traj_file)
NH3_atoms = NH3_trajectory[-1]
NH3_energy = NH3_atoms.calc.results['energy']  # in eV


# determine a, b, c, d for the adsorbate
a = composition['C'] if 'C' in composition else 0
b = composition['O'] if 'O' in composition else 0
c = composition['N'] if 'N' in composition else 0
d = composition['H'] if 'H' in composition else 0

# heat of formation at 
if adsorbate_label == 'H':
    # H2 reference gas
    working_rxn = 2.239  # eV

    heat_of_formation_0K = -0.5 * working_rxn + system_energy - slab_energy # in kJ/mol
# heat_of_formation_0K = binding_energy * eV_to_kJ + 239.2  # in kJ/mol


# heat_of_formation_0K = -259.05  # kJ/mol
# my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K, twoD_gas=True)
# a_low2, a_high2 = my_calc.get_thermo2()