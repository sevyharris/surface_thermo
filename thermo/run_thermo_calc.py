import os
import matplotlib.pyplot as plt
import numpy as np
import rmgpy.constants
import argparse
import adsorbate_thermo
# import rmgpy.species
import sys
sys.path.append('/home/moon/surface/surface_thermo')
import util
import yaml
import ase.io.trajectory
import logging

import rmgpy.thermo
import rmgpy.chemkin


logging.basicConfig(level=logging.INFO)

# script to convert binding energies to NASA polynomials
metal = 'Pt'
facet = '111'

if facet == '111':
    # sites = ['ontop']  # TODO, finish converging bridge for Pt111
    sites = ['ontop', 'hcp', 'bridge']  # TODO, finish converging bridge for Pt111
    # sites = ['ontop', 'hcp']
    # sites = ['hcp']
elif facet == '110':
    sites = ['ontop', 'hollow']


# get adsorbate label from input
parser = argparse.ArgumentParser(description='Run thermo calculation for adsorbate on metal facet.')
parser.add_argument('--adsorbate', type=str, default='H', help='Adsorbate label (default: H)')
args = parser.parse_args()
adsorbate_label = args.adsorbate



species_dictionary = rmgpy.chemkin.load_species_dictionary(
    '../my_dictionary.txt'
)

translator = {
    'H': 'HX',
    'H2': 'H2X',
    'O': 'OX',
    'H2O': 'H2OX',
    'CH3': 'CH3X',
    # 'CH4': 'CH4X',
    'CH': 'CHX',
    'CH2': 'CH2X',
    'C': 'CX',
    'CO': 'OCX',
    'CO2': 'CO2X',
    'NH3': 'NH3X',
    'NH2': 'H2NX',
    'NH': 'HNX',
    'N': 'NX',
    'N2': 'N2X',
}


# adsorbate_label = 'H2'
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

    # include the system ZPE
    system_vib_yaml_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', 
                                         f'{metal}{facet}_{adsorbate_label}_{site}_vib.yaml')
    with open(system_vib_yaml_file, 'r') as f:
        system_vib_data = yaml.load(f, Loader=yaml.FullLoader)  # Load the YAML file
    system_zpe = system_vib_data['zpe']  # in eV
    total_system_energy = system_energy + system_zpe  # in eV

    site_energies[i] = total_system_energy

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

# include the system ZPE
system_vib_yaml_file = os.path.join(results_dir, 'system', f'{metal}{facet}_{adsorbate_label}', 
                                     f'{metal}{facet}_{adsorbate_label}_{site}_vib.yaml')
with open(system_vib_yaml_file, 'r') as f:
    system_vib_data = yaml.load(f, Loader=yaml.FullLoader)
system_zpe = system_vib_data['zpe']  # in eV
total_system_energy = system_energy + system_zpe  # in eV


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
# include the gas ZPE
gas_vib_yaml_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}_vib.yaml')
with open(gas_vib_yaml_file, 'r') as f:
    gas_vib_data = yaml.load(f, Loader=yaml.FullLoader)
gas_zpe = gas_vib_data['zpe']  # in eV
total_gas_energy = gas_energy + gas_zpe  # in eV

# Calculate the binding energy
binding_energy = total_system_energy - (slab_energy + total_gas_energy)
binding_energy_no_zpe = system_energy - (slab_energy + gas_energy)

print('-------------- Binding Energies (with ZPE)----------------')
print(f'Gas energy for {adsorbate_label}: {total_gas_energy:.4f} eV')
print(f'Slab energy for {metal}{facet}: {slab_energy:.4f} eV')
print(f'System energy for {metal}{facet}_{adsorbate_label}_{site}: {total_system_energy:.4f} eV')
print()
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy:.4f} eV')
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy * eV_to_kJ:.4f} kJ/mol')
print()
print()
# print((binding_energy * eV_to_kJ + 239.2) / -239.2)
print('-------------- Binding Energies (no ZPE)----------------')
print(f'Gas energy for {adsorbate_label}: {gas_energy:.4f} eV')
print(f'Slab energy for {metal}{facet}: {slab_energy:.4f} eV')
print(f'System energy for {metal}{facet}_{adsorbate_label}_{site}: {system_energy:.4f} eV')
print()
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy_no_zpe:.4f} eV')
print(f'Binding energy for {adsorbate_label} on {metal}{facet} ({site}): {binding_energy_no_zpe * eV_to_kJ:.4f} kJ/mol')
print()
print()


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
composition = {'H': 0, 'C': 0, 'O': 0, 'N': 0}
for atom in gas:
    if atom.symbol not in composition:
        composition[atom.symbol] = 0
    composition[atom.symbol] += 1
print(f'Composition of {adsorbate_label}: {composition}')

#------------------- Grab the reference gases from DFT -------------------
# Reference gas for H2
H2_traj_file = os.path.join(results_dir, 'gas', f'H2.traj')
H2_trajectory = ase.io.trajectory.Trajectory(H2_traj_file, mode='r')
H2_atoms = H2_trajectory[-1]
H2_energy = H2_atoms.calc.results['energy']  # in eV
# get ZPE from yaml file
H2_vib_yaml_file = os.path.join(results_dir, 'gas', f'H2_vib.yaml')
with open(H2_vib_yaml_file, 'r') as f:
    H2_vib_data = yaml.load(f, Loader=yaml.FullLoader)
H2_frequencies = np.abs(H2_vib_data['frequencies'])
H2_zpe = H2_vib_data['zpe']  # in e
H2_total_energy = H2_energy + H2_zpe  # in eV


# Reference gas for H2O
H2O_traj_file = os.path.join(results_dir, 'gas', f'H2O.traj')
H2O_trajectory = ase.io.trajectory.Trajectory(H2O_traj_file, mode='r')
H2O_atoms = H2O_trajectory[-1]
H2O_energy = H2O_atoms.calc.results['energy']  # in eV
# get ZPE from yaml file
H2O_vib_yaml_file = os.path.join(results_dir, 'gas', f'H2O_vib.yaml')
with open(H2O_vib_yaml_file, 'r') as f:
    H2O_vib_data = yaml.load(f, Loader=yaml.FullLoader)
H2O_frequencies = np.abs(H2O_vib_data['frequencies'])
H2O_zpe = H2O_vib_data['zpe']  # in eV
H2O_total_energy = H2O_energy + H2O_zpe  # in eV


# Reference gas for CH4
CH4_traj_file = os.path.join(results_dir, 'gas', f'CH4.traj')
CH4_trajectory = ase.io.trajectory.Trajectory(CH4_traj_file, mode='r')
CH4_atoms = CH4_trajectory[-1]
CH4_energy = CH4_atoms.calc.results['energy']  # in eV
# get ZPE from yaml file
CH4_vib_yaml_file = os.path.join(results_dir, 'gas', f'CH4_vib.yaml')
with open(CH4_vib_yaml_file, 'r') as f:
    CH4_vib_data = yaml.load(f, Loader=yaml.FullLoader)
CH4_frequencies = np.abs(CH4_vib_data['frequencies'])
CH4_zpe = CH4_vib_data['zpe']  # in eV
CH4_total_energy = CH4_energy + CH4_zpe  # in eV


# Reference gas for NH3
NH3_traj_file = os.path.join(results_dir, 'gas', f'NH3.traj')
NH3_trajectory = ase.io.trajectory.Trajectory(NH3_traj_file, mode='r')
NH3_atoms = NH3_trajectory[-1]
NH3_energy = NH3_atoms.calc.results['energy']  # in eV
# get ZPE from yaml file
NH3_vib_yaml_file = os.path.join(results_dir, 'gas', f'NH3_vib.yaml')
with open(NH3_vib_yaml_file, 'r') as f:
    NH3_vib_data = yaml.load(f, Loader=yaml.FullLoader)
NH3_frequencies = np.abs(NH3_vib_data['frequencies'])
NH3_zpe = NH3_vib_data['zpe']  # in eV
NH3_total_energy = NH3_energy + NH3_zpe  # in eV


# determine a, b, c, d for the adsorbate
a = composition['C'] if 'C' in composition else 0
b = composition['O'] if 'O' in composition else 0
c = composition['N'] if 'N' in composition else 0
d = composition['H'] if 'H' in composition else 0

print(f'Adsorbate composition: C={a}, O={b}, N={c}, H={d}')
print()

# Get the ATcT heats of formation at 0K
Hf_CH4_ATcT = -66.556  # kJ/mol
Hf_H2O_ATcT = -238.938  # kJ/mol
Hf_NH3_ATcT = -38.565  # kJ/mol
Hf_H2_ATcT = 0.0  # kJ/mol


# Get the heat of reaction for the fictitious working reaction [1]
heat_of_working_rxn_dft = total_gas_energy - a * CH4_total_energy - b * H2O_total_energy - c * NH3_total_energy - \
     (d / 2.0 - 2.0 * a - b - 1.5 * c) * H2_total_energy  # in eV

print(f'Heat of working reaction (DFT) [2]: {heat_of_working_rxn_dft:.4f} eV')
print(f'Heat of working reaction (DFT) [2]: {heat_of_working_rxn_dft * eV_to_kJ:.4f} kJ/mol')
print()

Hf_gas_kJ_mol = a * Hf_CH4_ATcT + b * Hf_H2O_ATcT + c * Hf_NH3_ATcT + (d / 2.0 - 2.0 * a - b - 1.5 * c) * Hf_H2_ATcT + \
        heat_of_working_rxn_dft * eV_to_kJ  # in kJ/mol
print(f'Heat of formation (DFT) of {adsorbate_label} at 0K [3]: {Hf_gas_kJ_mol / eV_to_kJ:.4f} eV')
print(f'Heat of formation (DFT) of {adsorbate_label} at 0K [3]: {Hf_gas_kJ_mol:.4f} kJ/mol')
print()


# next up is binding energy
print(f'Heat of adsorption at 0K [4]: {binding_energy:.4f} eV')
print(f'Heat of adsorption at 0K [4]: {binding_energy * eV_to_kJ:.4f} kJ/mol')
print()

# Final heat of formation of system at 0K
heat_of_formation_0K_kJ_mol = Hf_gas_kJ_mol + binding_energy * eV_to_kJ
print(f'Final heat of formation of {adsorbate_label} on {metal}{facet} ({site}) at 0K [5]: {heat_of_formation_0K_kJ_mol / eV_to_kJ:.4f} eV')
print(f'Final heat of formation of {adsorbate_label} on {metal}{facet} ({site}) at 0K [5]: {heat_of_formation_0K_kJ_mol:.4f} kJ/mol')
print()

print(*frequencies, sep=", ")
# put the result through the adsorbate thermo calculator

# Using my refactor of Katrin's code
my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K_kJ_mol, twoD_gas=True)
a_low, a_high = my_calc.get_thermo2()
thermo_data = rmgpy.thermo.NASA(
    polynomials = [
        rmgpy.thermo.NASAPolynomial(coeffs=a_low, Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
        rmgpy.thermo.NASAPolynomial(coeffs=a_high, Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
    ],
    Tmin = (298.0, 'K'),
    Tmax = (3000.0, 'K'),
)
# print(f'Adsorbate thermo coefficients for {adsorbate_label} on {metal}{facet} ({site}):')
# print(thermo_data)
# print(f'a_low: {a_low}')
# print(f'a_high: {a_high}')
# print()


# Katrin's original format
# frequencies_str = [freq for freq in frequencies]
# frequencies_str.append('cm-1')  # Append the unit to the frequencies list
# frequencies_str = ', '.join(map(str, frequencies_str))  # Convert to a string for output
# frequencies_str = f'[{frequencies_str}]'  # Format as a list string

# # Save results to a file;
# my_result_file = os.path.join(results_dir, 'thermo', f'{adsorbate_label}-ads.dat')
# with open(my_result_file, 'w') as f:
#     f.write(f'name = {adsorbate_label}_ads\n')
#     f.write(f'DFT_binding_energy = [{binding_energy:.4f}, eV]\n')
#     f.write(f'heat_of_formation_0K = [{heat_of_formation_0K_kJ_mol:.4f}, kJ/mol]\n')
#     f.write(f'composition = {composition}\n')
#     f.write(f'sites = 1\n')
#     f.write(f'adsorbate_mass = [{molecular_weight:.4f}, amu]\n')
#     f.write(f'linear_scaling_binding_atom = [-2.479, eV]\n')
#     f.write(f'linear_scaling_gamma(X) = [1.0] \n')
#     f.write(f'linear_scaling_psi = [0, eV]\n')
#     f.write(f'frequencies = {frequencies_str}\n')

# also save as yaml file
my_result_yaml = os.path.join(results_dir, 'thermo', f'{adsorbate_label}-ads.yaml')
results_dict = {
    'name': f'{adsorbate_label}_ads',
    'DFT_binding_energy': [binding_energy, 'eV'],
    'heat_of_formation_0K': [heat_of_formation_0K_kJ_mol, 'kJ/mol'],
    'composition': composition,
    'sites': 1,
    'adsorbate_mass': [float(molecular_weight), 'amu'],
    'linear_scaling_binding_atom': [-2.479, 'eV'],
    'linear_scaling_gamma(X)': [1.0],
    'linear_scaling_psi': [0, 'eV'],
    'frequencies': [float(f) for f in frequencies], # Convert frequencies to float and append unit
    'adjacency_list': species_dictionary[translator[adsorbate_label]].to_adjacency_list() if adsorbate_label in translator else None,
}
with open(my_result_yaml, 'w') as f:
    yaml.dump(results_dict, f, default_flow_style=False)

print(f'Results saved to {os.path.abspath(my_result_yaml)}')