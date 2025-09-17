# script to print the binding energy block
import argparse
import yaml
import os

# take in the metal and facet as arguments
parser = argparse.ArgumentParser(description='Print binding energy block for a system.')
parser.add_argument('--metal', type=str, required=True, help='Metal to use for the calculation')
parser.add_argument('--facet', type=str, required=True, help='Facet of the metal slab')


# parse the argumentsargs = parser.parse_args()
args = parser.parse_args()
metal = args.metal
facet = args.facet


# binding energy is the heat of formation of the system minus the ATcT heat of formation of the gas
ATcT_energies = {
    'H': 2.239042338,  # eV
    'C': 7.373000985,
    'O': 2.558366585,
    'N': 4.877566461,
}

# read in the H and metal heat of formation from the yaml file
print('bindingEnergies = {')
for atom in ['H', 'C', 'O', 'N']:
    thermo_file = os.path.join(os.path.dirname(__file__), f'../results/thermo/{metal}{facet}/{metal}{facet}_{atom}-ads.yaml')
    with open(thermo_file, 'r') as f:
        data = yaml.safe_load(f)

        heat_of_formation_0K = data['heat_of_formation_0K'][0]
        units = data['heat_of_formation_0K'][1]
        assert units == 'kJ/mol', f"Expected units 'kJ/mol', got {units}"

        binding_energy = heat_of_formation_0K / 96.485 - ATcT_energies[atom]
        print(f"\t'{atom}': ({binding_energy}, 'eV/molecule'),")
print('}')


# output format like this 
bindingEnergies = {
    'C': (-6.102800985, 'eV/molecule'),
    'H': (-2.750942338, 'eV/molecule'),
    'O': (-6.130366585, 'eV/molecule'),
    'N': (-5.987266461, 'eV/molecule'),
}
