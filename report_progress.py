import os
import ase.io.trajectory
import util

# sets to check
my_systems = [
    {'slabname': 'Cr2O3_z'},
    {'slabname': 'Fe2O3_z'},
    {'metal': 'Pt', 'crystal_structure': 'fcc', 'facet': '111'},
    {'metal': 'Pt', 'crystal_structure': 'fcc', 'facet': '100'},
    {'metal': 'Cr', 'crystal_structure': 'bcc', 'facet': '100'},
    {'metal': 'Cr', 'crystal_structure': 'bcc', 'facet': '110'},
    {'metal': 'Fe', 'crystal_structure': 'bcc', 'facet': '100'},
    {'metal': 'Fe', 'crystal_structure': 'fcc', 'facet': '111'},
]

important_adsorbates = ['C', 'H', 'O', 'N']

# slab completed?
for i in range(len(my_systems)):
    if 'slabname' in my_systems[i]:
        slabname = my_systems[i]['slabname']
        site_names = util.get_site_names(slabname=slabname)
    else:
        metal = my_systems[i]['metal']
        crystal_structure = my_systems[i]['crystal_structure']
        facet = my_systems[i]['facet']
        slabname = f'{metal}_{crystal_structure}{facet}'
        site_names = util.get_site_names(crystal_structure=crystal_structure, facet=facet)

    slab_traj_file = os.path.join('results', 'slab', f'{slabname}_slab.traj')
    if not os.path.exists(slab_traj_file):
        print(f'Slab {slabname} incomplete')

    traj = ase.io.trajectory.Trajectory(slab_traj_file)
    if not util.atoms_converged(traj[-1], fmax=0.05):
        print(f'Slab {slabname} not converged')

    # check on the adsorbates
    for j in range(len(important_adsorbates)):
        ads = important_adsorbates[j]
        rotations = [0.0, 90.0, 180.0]
        if len(ads) == 1:
            rotations = [0.0]  # skip other rotations for monatomic species

        for k in range(len(site_names)):
            site_name = site_names[k]
            for r in range(len(rotations)):
                rot = rotations[r]
                system_traj_file = os.path.join(
                    'results', 'system', f'{slabname}_{ads}', f'{slabname}_{ads}_{site_name}_rot{rot}.traj'
                )
                if not os.path.exists(system_traj_file):
                    print(f'Missing system traj file {system_traj_file}')
                    continue
                traj = ase.io.trajectory.Trajectory(system_traj_file)
                if not util.atoms_converged(traj[-1], fmax=0.05):
                    print(f'System {system_traj_file} not converged')
                    # probably fine if it doesn't converge. Indicates it's not stable

                system_vib_file = system_traj_file.replace('.traj', '_vib.yaml')
                if not os.path.exists(system_vib_file):
                    print(f'Missing system vib file {system_vib_file}')
