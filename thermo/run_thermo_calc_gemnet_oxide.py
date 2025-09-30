# script to compile important information from adsorbate quantum chemistry calculation

import os
import matplotlib.pyplot as plt
import numpy as np
import rmgpy.constants
import argparse
import adsorbate_thermo
import shutil
# import rmgpy.species
import sys

sys.path.append(os.environ['SURFACE_THERMO_DIR'])
import util
import yaml
import ase.io.trajectory
import ase.build
import logging

import rmgpy.thermo
import rmgpy.chemkin

import matplotlib


logging.basicConfig(level=logging.INFO)


# get adsorbate label from input
parser = argparse.ArgumentParser(description='Run thermo calculation for adsorbate on oxide facet.')
parser.add_argument('--adsorbate', type=str, default='H', help='Adsorbate label (default: H)')
parser.add_argument('--slabname', type=str, default='Cr2O3_z', help='Oxide label (default: Cr2O3_z)')


args = parser.parse_args()
adsorbate_label = args.adsorbate
slabname = args.slabname

fmax = 0.05

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
sites = [i for i in range(len(sites_dict[slabname]))]


species_dictionary = rmgpy.chemkin.load_species_dictionary(
    '../my_dictionary.txt'
)

translator = {  # going from adsorbate name in my workflow to the species dictionary label of the adsorbate
    'H': 'HX',
    'H2': 'H2X',
    'O': 'OX',
    'OH': 'HOX',
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
    'O2': 'O2X',
    'NO': 'NOX',
}


results_dir = f'../results'
my_result_yaml = os.path.join(results_dir, 'thermo', slabname, f'{slabname}_{adsorbate_label}-ads.yaml')
if not os.path.exists(os.path.dirname(my_result_yaml)):
    os.makedirs(os.path.dirname(my_result_yaml), exist_ok=True)


eV_to_kJ = 96.485332

# -------------- Gather the binding energies ----------------
possible_rotations = ['0.0', '90.0', '180.0']  # in degrees

# Pick the site with the minimum energy
site_energies = np.zeros(len(sites))
site_rotations = np.zeros(len(sites), dtype=int)  # to store the rotation for each site
for i, site in enumerate(sites):
    rot_energies = 1e5 + np.zeros(3)  # for 0, 90, 180 degree rotations
    for j in range(3):
        # load the system from the trajectory file
        system_traj_file = os.path.join(results_dir, 'system', f'{slabname}_{adsorbate_label}',
                                        f'{slabname}_{adsorbate_label}_{site}_rot{float(j * 90.0):.1f}.traj')
        if not os.path.exists(system_traj_file):
            logging.warning(f'System trajectory file {system_traj_file} does not exist. Skipping site {site}.')
            continue

        system_trajectory = ase.io.trajectory.Trajectory(system_traj_file)
        system = system_trajectory[-1]
        if not util.atoms_converged(system, fmax=fmax):
            logging.warning(f'System {slabname}_{adsorbate_label}_{site}_rot{j * 90.0} is not converged')
            continue

        # reject systems where the relaxation took too many steps. This is a sign of restructuring
        if len(system_trajectory) > 100:
            logging.warning(f'System {slabname}_{adsorbate_label}_{site}_rot{j * 90.0} took too many steps: {len(system_trajectory)}')
            continue

        # reject systems where the adsorbate is below the surface
        slab_traj_file = os.path.join(results_dir, 'slab', f'{slabname}_slab.traj')
        slab_traj = ase.io.trajectory.Trajectory(slab_traj_file)
        slab = slab_traj[-1]

        metal_atoms = [atom for atom in system[:len(slab)] if atom.symbol not in ['C', 'H', 'O', 'N']]
        highest_metal_z = np.max([atom.position[2] for atom in metal_atoms])  # highest z-coordinate of metal atoms
        # adsorbate_atoms = [atom for atom in system if atom.symbol != metal]
        adsorbate_z = np.min([atom.position[2] for atom in system[len(slab):]])  # lowest z-coordinate of adsorbate atoms
        if adsorbate_z < highest_metal_z:
            logging.warning(f'Adsorbate {adsorbate_label} is below the surface for {slabname}_{adsorbate_label}_{site}_rot{j * 90}. Skipping.')
            continue

        # reject systems where the adsorbate is too far away from the surface
        THRESHOLD = 3.0  # in Angstroms
        # if vdW, use a larger threshold
        if util.is_vdW_species(species_dictionary[translator[adsorbate_label]]):
            THRESHOLD = 5.0

        if adsorbate_z - highest_metal_z > THRESHOLD:
            logging.warning(f'Adsorbate {adsorbate_label} is too far away from the surface for {slabname}_{adsorbate_label}_{site}_rot{j * 90}. Skipping.')
            continue

        # reject vdW species that are too close to the surface
        if util.is_vdW_species(species_dictionary[translator[adsorbate_label]]):
            if adsorbate_z - highest_metal_z < 2.0:
                logging.warning(f'Adsorbate {adsorbate_label} is too close to the surface for {slabname}_{adsorbate_label}_{site}_rot{j * 90}. Skipping.')
                continue

        # Reject systems where the adsorbate has fallen apart
        if not util.adsorbate_intact(system, slab, adsorbate_label):
            logging.warning(f'Adsorbate {adsorbate_label} has fallen apart for {slabname}_{adsorbate_label}_{site}_rot{j * 90}. Skipping.')
            continue

        system_energy = system.calc.results['energy']  # in eV

        # include the system ZPE
        system_vib_yaml_file = os.path.join(results_dir, 'system', f'{slabname}_{adsorbate_label}',
                                            f'{slabname}_{adsorbate_label}_{site}_rot{j * 90.0:.1f}_vib.yaml')
        with open(system_vib_yaml_file, 'r') as f:
            system_vib_data = yaml.load(f, Loader=yaml.FullLoader)  # Load the YAML file

        system_zpe = system_vib_data['zpe']  # in eV
        total_system_energy = system_energy + system_zpe  # in eV
        rot_energies[j] = total_system_energy

    site_energies[i] = np.min(rot_energies)
    site_rotations[i] = np.argmin(rot_energies)  # store the rotation with the minimum energy


# Print the minimum energy site
min_energy_site = sites[np.argmin(site_energies)]
print(f'Minimum energy site for {adsorbate_label} on {slabname}: {min_energy_site} {np.min(site_energies)} (includes ZPE)')

rotation = site_rotations[np.argmin(site_energies)] * 90.0
print(f'adsorbate was rotated by {rotation} degrees')


# reload the system with the minimum energy site
system_traj_file = os.path.join(results_dir, 'system', f'{slabname}_{adsorbate_label}',
                                f'{slabname}_{adsorbate_label}_{min_energy_site}_rot{rotation}.traj')
system_trajectory = ase.io.trajectory.Trajectory(system_traj_file)
system = system_trajectory[-1]
site = min_energy_site  # update the site to the minimum energy site
system_energy = system.calc.results['energy']  # in eV


# Save a picture of the relaxed system
side_pic = os.path.join(results_dir, 'thermo', slabname, 'images', f'{slabname}_{adsorbate_label}_{site}_rot{rotation}_side.png')
top_pic = os.path.join(results_dir, 'thermo', slabname, 'images', f'{slabname}_{adsorbate_label}_{site}_rot{rotation}_top.png')
if not os.path.exists(os.path.dirname(side_pic)):
    os.makedirs(os.path.dirname(side_pic), exist_ok=True)

if str(matplotlib.__version__) != '3.9.4':
    ase.io.write(side_pic, system, rotation='-90x,0y,0z')
    ase.io.write(top_pic, system, rotation='0x,0y,0z')
else:
    existing_side_pic = os.path.join(
        results_dir, 'system', f'{slabname}_{adsorbate_label}', f'{slabname}_{adsorbate_label}_{site}_rot{rotation}_side.png'
    )
    existing_top_pic = os.path.join(
        results_dir, 'system', f'{slabname}_{adsorbate_label}', f'{slabname}_{adsorbate_label}_{site}_rot{rotation}_top.png'
    )
    shutil.copyfile(existing_side_pic, side_pic)
    shutil.copyfile(existing_top_pic, top_pic)


# include the system ZPE
system_vib_yaml_file = os.path.join(results_dir, 'system', f'{slabname}_{adsorbate_label}',
                                    f'{slabname}_{adsorbate_label}_{site}_rot{rotation}_vib.yaml')
with open(system_vib_yaml_file, 'r') as f:
    system_vib_data = yaml.load(f, Loader=yaml.FullLoader)
system_zpe = system_vib_data['zpe']  # in eV
total_system_energy = system_energy + system_zpe  # in eV


# get the slab energy
slab_traj_file = os.path.join(results_dir, 'slab', f'{slabname}_slab.traj')
slab_trajectory = ase.io.trajectory.Trajectory(slab_traj_file)
slab = slab_trajectory[-1]
slab_energy = slab.calc.results['energy']  # in eV
print(f'Slab energy for {slabname}: {slab_energy:.4f} eV')
# # Gather the gas energy
gas = ase.build.molecule(adsorbate_label)


# -------------- Convert to NASA polynomials ----------------
# get molecular weight for gas
molecular_weight = np.sum(gas.get_masses())  # in amu
print(f'Molecular weight of {adsorbate_label}: {molecular_weight:.4f} amu')

# read yaml
with open(system_vib_yaml_file, 'r') as f:
    vib_data = yaml.load(f, Loader=yaml.FullLoader)

# get magnitude to avoid imaginary components
frequencies = np.abs(vib_data['frequencies'])
print(f'Vibrational frequencies for {slabname}_{adsorbate_label}_{site}: {frequencies} (cm^-1)')

# build the composition dictionary
composition = {'H': 0, 'C': 0, 'O': 0, 'N': 0}
for atom in gas:
    if atom.symbol not in composition:
        composition[atom.symbol] = 0
    composition[atom.symbol] += 1
print(f'Composition of {adsorbate_label}: {composition}')

# Reference gas energies:
e_C = -7.282  # eV
e_O = -7.204  # eV
e_N = -8.083  # eV
e_H = -3.477  # eV

# Reference ZPEs:
zpe_CO = 0.147  # eV
zpe_H2O = 0.609
zpe_N2 = 0.1458  # not real yet
zpe_H2 = 0.270  # chao wrote 0.272  # eV

# ATcT heats of formation at 0K
Hf_CO_ATcT = -113.799  # kJ/mol
Hf_H2O_ATcT = -238.938  # kJ/mol
Hf_N2_ATcT = 0.0  # kJ/mol
Hf_H2_ATcT = 0.0  # kJ/mol

# determine a, b, c, d for the adsorbate
a = composition['C'] if 'C' in composition else 0
b = composition['O'] if 'O' in composition else 0
c = composition['N'] if 'N' in composition else 0
d = composition['H'] if 'H' in composition else 0

print(f'Adsorbate composition: C={a}, O={b}, N={c}, H={d}')
print()

DFT_CO_total_energy = e_C + e_O + zpe_CO        # in eV
DFT_H2O_total_energy = e_H * 2 + e_O + zpe_H2O  # in eV
DFT_N2_total_energy = e_N * 2 + zpe_N2          # in eV
DFT_H2_total_energy = e_H * 2 + zpe_H2          # in eV

print()
print(DFT_H2O_total_energy)

# Get the DFT heat of reaction for the fictitious working reaction
heat_of_working_rxn_dft = -a * DFT_CO_total_energy - (b - a) * DFT_H2O_total_energy - \
    (d / 2.0 - b + a) * DFT_H2_total_energy - (c / 2.0) * DFT_N2_total_energy

print(f'Heat of working reaction (DFT) [2]: {heat_of_working_rxn_dft:.4f} eV')
print(f'Heat of working reaction (DFT) [2]: {heat_of_working_rxn_dft * eV_to_kJ:.4f} kJ/mol')
print()

Hf_gas_kJ_mol = a * Hf_CO_ATcT + (b - a) * Hf_H2O_ATcT + (d / 2.0 - b + a) * Hf_H2_ATcT + \
        (c / 2.0) * Hf_N2_ATcT + \
        heat_of_working_rxn_dft * eV_to_kJ  # in kJ


print(f'Heat of formation (DFT) of {adsorbate_label} at 0K [3]: {Hf_gas_kJ_mol / eV_to_kJ:.4f} eV')
print(f'Heat of formation (DFT) of {adsorbate_label} at 0K [3]: {Hf_gas_kJ_mol:.4f} kJ/mol')
print()


# Final heat of formation of system at 0K
heat_of_formation_0K_kJ_mol = Hf_gas_kJ_mol + (total_system_energy - slab_energy) * eV_to_kJ
print(f'Final heat of formation of {adsorbate_label} on {slabname} ({site}) at 0K [5]: {heat_of_formation_0K_kJ_mol / eV_to_kJ:.4f} eV')
print(f'Final heat of formation of {adsorbate_label} on {slabname} ({site}) at 0K [5]: {heat_of_formation_0K_kJ_mol:.4f} kJ/mol')
print()

print(*frequencies, sep=", ")
# put the result through the adsorbate thermo calculator

# Using my refactor of Katrin's code
my_calc = adsorbate_thermo.AdsorbateThermoCalc(molecular_weight, frequencies, composition, heat_of_formation_0K_kJ_mol, twoD_gas=False)
a_low, a_high = my_calc.get_thermo2()
thermo_data = rmgpy.thermo.NASA(
    polynomials = [
        rmgpy.thermo.NASAPolynomial(coeffs=a_low, Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
        rmgpy.thermo.NASAPolynomial(coeffs=a_high, Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
    ],
    Tmin = (298.0, 'K'),
    Tmax = (3000.0, 'K'),
)


print()
print(thermo_data)

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
my_result_yaml = os.path.join(results_dir, 'thermo', slabname, f'{slabname}_{adsorbate_label}-ads.yaml')
results_dict = {
    'name': f'{adsorbate_label}_ads',
    'DFT_binding_energy': [0, 'eV'],
    'heat_of_formation_0K': [heat_of_formation_0K_kJ_mol, 'kJ/mol'],
    'composition': composition,
    'sites': 1,
    'adsorbate_mass': [float(molecular_weight), 'amu'],
    'linear_scaling_binding_atom': [-2.479, 'eV'],
    'linear_scaling_gamma(X)': [1.0],
    'linear_scaling_psi': [0, 'eV'],
    'frequencies': [float(f) for f in frequencies],  # Convert frequencies to float and append unit
    'adjacency_list': species_dictionary[translator[adsorbate_label]].to_adjacency_list() if adsorbate_label in translator else None,
}
with open(my_result_yaml, 'w') as f:
    yaml.dump(results_dict, f, default_flow_style=False)

print(f'Results saved to {os.path.abspath(my_result_yaml)}')
