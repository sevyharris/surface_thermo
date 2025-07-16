import os
import numpy as np
import sys
import yaml

import ase.optimize
import ase.collections
import ase.vibrations
import logging

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils


# set up logging
logging.basicConfig(level=logging.INFO)

if len(sys.argv) < 2:
    adsorbate_label = 'H2'  # Default adsorbate if not provided
else:
    adsorbate_label = sys.argv[1]


# initialize fairchem ocp calculator
checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    # 'EquiformerV2-31M-S2EF-OC20-All+MD',
    # 'gnoc_oc22_oc20_all_s2ef.pt',
    'GemNet-OC-S2EFS-OC20+OC22',
    local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
)
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(
    checkpoint_path=checkpoint_path,
    cpu=True,
    seed=400
)


results_dir = '../results'  # Directory to save results

# optimize the adsorbate by itself
adsorbate = ase.collections.g2[adsorbate_label]
adsorbate.set_cell([10, 10, 10])  # Set a large cell to avoid interactions
adsorbate.set_pbc(True)  # Set periodic boundary conditions
adsorbate.calc = calc

# Optimize the adsorbate by itself
logging.info(f'Running optimization for {adsorbate_label}')
trajectory_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}.traj')
if not os.path.exists(os.path.dirname(trajectory_file)):
    os.makedirs(os.path.dirname(trajectory_file))
opt = ase.optimize.BFGS(adsorbate, trajectory=trajectory_file, append_trajectory=True)
opt.run(fmax=0.01)

# Run vibrational analysis
vib = ase.vibrations.Vibrations(adsorbate)
# vib = ase.vibrations.Vibrations(adsorbate, name=adsorbate_label)
vib.run()
vib.summary()
freq = vib.get_frequencies()
logging.info(f'Vibrational frequencies for {adsorbate_label}: {freq} (cm^-1)')
logging.info(f'ZPE for {adsorbate_label}: {vib.get_zero_point_energy()} eV')

# save result as a yaml file
result = {
    'frequencies': freq.tolist(),
    'zpe': float(vib.get_zero_point_energy()),
}
result_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}_vib.yaml')
with open(result_file, 'w') as f:
    yaml.dump(result, f, default_flow_style=False)
