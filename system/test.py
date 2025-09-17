import ase.io.trajectory
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import util

# trajectory_file = '/home/moon/surface/surface_thermo/results/slab/Cr110_slab.traj'
# trajectory_file = '/home/moon/surface/surface_thermo/results/slab/Fe110_slab.traj'
trajectory_file = '/home/moon/surface/surface_thermo/results/slab/Pt111_slab.traj'

traj = ase.io.trajectory.Trajectory(trajectory_file)
slab = traj[-1]  # Get the last frame from the trajectory

if util.atoms_converged(slab):
    print("Slab is converged")
else:
    print("Slab is not converged")

# print(len(slab))
# # ase.io.write('/home/moon/surface/surface_thermo/results/slab/Cr110_slab_side.png', slab, rotation="-90x,0y,0z")
# # ase.io.write('/home/moon/surface/surface_thermo/results/slab/Pt111_slab_side.png', slab, rotation="-90x,0y,0z")

# ase.io.write('/home/moon/surface/surface_thermo/results/slab/Fe110_slab.png', slab, rotation="0x,0y,0z")
# ase.io.write('/home/moon/surface/surface_thermo/results/slab/Fe110_slab_side.png', slab, rotation="-90x,0y,0z")


# # print(slab.calc.results['forces'])
# # print(slab.constraints[0].get_indices())



# Gather the gas energy
# adsorbate_label = 'H'
# results_dir = '../results'
# gas_traj_file = os.path.join(results_dir, 'gas', f'{adsorbate_label}.traj')
# gas_trajectory = ase.io.trajectory.Trajectory(gas_traj_file)
# for i in range(len(gas_trajectory)):
#     print(i, gas_trajectory[i].calc.results['energy'])
# gas = gas_trajectory[-1]
# gas_energy = gas.calc.results['energy']  # in 
# print(gas_energy)
