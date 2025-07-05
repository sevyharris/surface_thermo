# surface_thermo
Using ML energies to compute adsorbate thermo

# Installation:
1. Make a virtual environment with Python 3.12
2. Install fairchem-core1.1 into that environment (I had to manually install pytorch-scatter and pytorch-sparse)

# Instructions:
1. Run bulk optimizer to get lattice constant of relaxed metal
2. Run slab optimizer to get relaxed slab object
3. Run system optimizer to adsorb species onto the slab
4. Run thermo converter to get NASA polynomial
5. Do a sanity check with something from literature
