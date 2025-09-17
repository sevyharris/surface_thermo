import os
import numpy as np
import copy
import yaml
import argparse

import ase.build
import ase.optimize

import ase.io.trajectory
import ase.constraints
import logging

import fairchem.core.models.model_registry
import fairchem.core.common.relaxation.ase_utils

import matplotlib.pyplot as plt



metal = 'Fe'
crystal_structure = 'bcc'
facet = '110'
plotting = True
lattice_constant = None  # This will hold the final lattice constant after analysis
results_dir = '../results'  # Directory to save results
# set up logging
logging.basicConfig(level=logging.INFO)



checkpoint_path = fairchem.core.models.model_registry.model_name_to_local_file(
    # 'EquiformerV2-31M-S2EF-OC20-All+MD',
    'GemNet-OC-S2EFS-OC20+OC22',
    local_cache='/home/moon/surface/tmp/fairchem_checkpoints/'
)
calc = fairchem.core.common.relaxation.ase_utils.OCPCalculator(checkpoint_path=checkpoint_path, cpu=True, seed=400)


for lattice_constant in np.linspace(2.6, 2.8, 11):
    print(f"Running relaxation for {metal} {crystal_structure} with lattice constant {lattice_constant:.3f} Å")
    fmax = 0.01
    vacuum = 7.5

    slab = ase.build.bcc110(metal, size=(3, 3, 4), vacuum=vacuum, a=lattice_constant)
 
    # Fix the bottom two layers
    bottom_layer = []
    second_layer = []
    fixed_indices = []
    z_values = list(set([pos[2] for pos in slab.get_positions()]))
    z_values.sort()

    for i, pos in enumerate(slab.get_positions()):
        if pos[2] == z_values[0]:
            bottom_layer.append(slab[i].index)
        if pos[2] == z_values[1]:
            second_layer.append(slab[i].index)
    fixed_indices = bottom_layer + second_layer
    fix_bottom_layers = ase.constraints.FixAtoms(indices=fixed_indices)
    slab.set_constraint(fix_bottom_layers)

    slab.set_initial_magnetic_moments([2.0] * len(slab))

    slab.calc = calc
    opt = ase.optimize.BFGS(slab, logfile=f'Fe_{lattice_constant}.log')
    opt.run(fmax=fmax, steps=20)

# # save a picture of the relaxed slab
# ase.io.write(os.path.join(results_dir, 'slab', f'{metal}{facet}_slab_top.png'), slab, rotation="0x,0y,0z")
# ase.io.write(os.path.join(results_dir, 'slab', f'{metal}{facet}_slab_side.png'), slab, rotation="-90x,0y,0z")
