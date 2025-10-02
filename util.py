import numpy as np
import scipy.cluster.vq
import ase.build
import ase.io.trajectory
import ase.geometry.analysis


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
    if isinstance(trajectory_file, str):
        traj = ase.io.trajectory.Trajectory(trajectory_file, mode='r')
    elif isinstance(trajectory_file, ase.io.trajectory.Trajectory):
        # If it's already a Trajectory object, use it directly
        traj = trajectory_file

    return traj[-1].calc.results['energy'] if traj else None


def adsorbate_intact(system, slab, adsorbate_label):
    """
    Check if the adsorbate is intact in the system.

    Parameters:
    system (ase.Atoms): The system containing the adsorbate.
    adsorbate_label (str): The label of the adsorbate to check.

    Returns:
    bool: True if the adsorbate is intact, False otherwise.
    """

    original_adsorbate = ase.build.molecule(adsorbate_label)
    original_analysis = ase.geometry.analysis.Analysis(original_adsorbate)
    original_count = np.sum([len(x) for x in original_analysis.all_bonds[0]])

    # # assume metal is most popular element in the system
    # metal = max(set(atom.symbol for atom in system), key=lambda x: sum(atom.symbol == x for atom in system))
    # adsorbate_indices = [i for i in range(len(system)) if system[i].symbol != metal]

    # adsorbate = ase.Atoms()
    # for i in adsorbate_indices:
    #     adsorbate += system[i]

    analysis = ase.geometry.analysis.Analysis(system[len(slab):])
    count = np.sum([len(x) for x in analysis.all_bonds[0]])
    return count == original_count


def is_vdW_species(species):
    """
    Check if the species is a van der Waals species.

    Parameters:
    species (str): The species to check.

    Returns:
    bool: True if the species is a van der Waals species, False otherwise.
    """
    if not species.contains_surface_site():
        return False
    bonded_to_surface = False
    for atom in species.molecule[0].atoms:
        if atom.is_bonded_to_surface():
            bonded_to_surface = True
    return not bonded_to_surface


def enumerate_layers(atoms):
    """
    Enumerate the layers in the slab based on the z-coordinates of the atoms.

    Parameters:
    atoms (ase.Atoms): The slab atoms.

    Returns:
    list: A list of layer indices.
    """

    def compute_error(atoms, centroids, labels):
        a = 1.0  # characteristic length between layers
        error = 0
        assert len(atoms) == len(labels)
        max_mass = np.max(atoms.get_masses())
        for i in range(len(atoms)):
            error += np.float_power(atoms[i].position[2] - centroids[labels[i]], 2.0) * atoms.get_masses()[i] / max_mass

        # get centroid's nearest neighbor
        centroid_min_distances = np.zeros_like(centroids) + np.inf
        for i in range(len(centroids)):
            for j in range(len(centroids)):
                if i == j:
                    continue
                if np.abs(centroids[i] - centroids[j]) < centroid_min_distances[i]:
                    centroid_min_distances[i] = np.abs(centroids[i] - centroids[j])

        # add penalty for non-uniformity and deviation from characteristic length
        nonuniformity_error = np.sum(np.float_power(centroid_min_distances - np.median(centroid_min_distances), 2.0))
        characteristic_length_error = 10 * np.sum(np.float_power(centroid_min_distances - a, 2.0))
        error += nonuniformity_error + characteristic_length_error
        return error

    Zs = np.array([atom.position[2] for atom in atoms])

    nlayers = np.arange(2, 32)
    errors = np.zeros_like(nlayers) + np.inf
    for i, nlayer in enumerate(nlayers):
        centroid, label = scipy.cluster.vq.kmeans2(Zs, nlayer)

        # if a bin only contains 1, break here
        bin_counts = np.zeros(nlayer)
        for j in range(len(Zs)):
            bin_counts[label[j]] += 1

        if 1 in bin_counts or 0 in bin_counts:
            break

        errors[i] = compute_error(atoms, centroid, label)

    best_nlayer = nlayers[np.argmin(errors)]
    centroid, label = scipy.cluster.vq.kmeans2(Zs, best_nlayer)
    print(centroid)

    # sort the labels by the centroid positions
    indices = np.arange(len(centroid))
    sorted_indices = [x for _, x in sorted(zip(centroid, indices))][::-1]

    new_centroid = [centroid[i] for i in sorted_indices]

    new_label = [new_centroid.index(centroid[i]) for i in label]

    # print(label)
    # label = np.array([sorted_indices[i] for i in label])
    # print(label)
    return list(np.array(new_label) + 1)  # make it 1-indexed


def get_site_names(crystal_structure=None, facet=None, slabname=None):
    sites_dict = {
        'Cr2O3_z': [
            (0.0, 0.0),
            (2.5018702100000003, 1.4444554392210054),
            (9.067513484506406e-16, 2.8889108784420108),
            (0.7547939021521936, 4.19625226621278),
            (0.7547939021521931, 1.5815694906712414),
            (1.74707631e+00, 2.75179683e+00)
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
    if slabname is not None:
        return [i for i in range(len(sites_dict[slabname]))]

    assert crystal_structure is not None
    assert facet is not None

    sites = []
    tmp_slab = None
    if crystal_structure == 'fcc':
        if facet == '111':
            tmp_slab = ase.build.fcc111('Pt', (3, 3, 4))
        elif facet == '110':
            tmp_slab = ase.build.fcc110('Pt', (3, 3, 4))
        elif facet == '100':
            tmp_slab = ase.build.fcc100('Pt', (3, 3, 4))
    elif crystal_structure == 'bcc':
        if facet == '111':
            tmp_slab = ase.build.bcc111('Fe', (3, 3, 4), a=3.0)
        elif facet == '110':
            tmp_slab = ase.build.bcc110('Fe', (3, 3, 4), a=3.0)
        elif facet == '100':
            tmp_slab = ase.build.bcc100('Fe', (3, 3, 4), a=3.0)
    assert tmp_slab is not None, 'unrecognized crystal structure/facet combo'
    sites = list(tmp_slab.info['adsorbate_info']['sites'].keys())
    return sites


def get_max_surface_displacement(slab_init, system_end):
    # check whether the surface has rearranged itself after optimizing the system
    N_atoms = len(slab_init)
    return np.max(np.linalg.norm(system_end[:N_atoms].positions - slab_init.positions, axis=1))