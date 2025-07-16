import numpy as np
import ase.io.trajectory

# check converged
def atoms_converged(atoms, fmax=0.01):
    """
    Check if the atoms are converged based on the maximum force. Only considers the force on uncontrained atoms.
    
    Parameters:
    atoms (ase.Atoms): The atoms object to check.
    fmax (float): The maximum force threshold for convergence.
    
    Returns:
    bool: True if converged, False otherwise.
    """
    # check if the slab is already relaxed
    fixed_indices = atoms.constraints[0].get_indices()
    final_force = atoms.calc.results['forces'] if 'forces' in atoms.calc.results else np.inf
    # only count the forces on the atoms that are not fixed
    final_force = [final_force[i, :] if i not in fixed_indices else np.zeros(3) for i in range(final_force.shape[0])]
    max_force = np.max(np.linalg.norm(final_force, axis=1))
    return max_force < fmax

def get_final_energy(trajectory_file):
    """
    Get the final energy from the trajectory file.
    
    Parameters:
    trajectory_file (str): The path to the trajectory file.
    
    Returns:
    float: The final energy in eV.
    """
    if type(trajectory_file) == str:
        traj = ase.io.trajectory.Trajectory(trajectory_file, mode='r')
    elif type(trajectory_file) == ase.io.trajectory.Trajectory:
        # If it's already a Trajectory object, use it directly
        traj = trajectory_file
    
    return traj[-1].calc.results['energy'] if traj else None
