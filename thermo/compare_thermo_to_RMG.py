# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# +
# Script to compare my calculated results to what's in RMG-database

# +
import os
import rmgpy.kinetics
import rmgpy.data.thermo
import numpy as np

import matplotlib.pyplot as plt
# %matplotlib inline

# -

def plot_thermos(thermos, labels=None):
    if type(thermos) != list:
        thermos = [thermos]
    if labels is None:
        labels = ['' for t in thermos]
    linestyles = ['solid', 'solid', 'dashed', 'dashed']
    fig, ax = plt.subplots(1, 3)
    fig.set_size_inches(12, 3)
    fig.tight_layout()
    ax[0].set_xlabel('Temperature (K)')
    ax[0].set_ylabel('H (kJ / mol)')
    ax[0].set_title('Enthalpy vs. Temperature')
    ax[1].set_xlabel('Temperature (K)')
    ax[1].set_ylabel('S (kJ / mol K)')
    ax[1].set_title('Entropy vs. Temperature')
    ax[2].set_xlabel('Temperature (K)')
    ax[2].set_ylabel('Cp (kJ / mol K)')
    ax[2].set_title('Heat Capacity vs. Temperature')
    T = np.linspace(300, 2000, 1001)
#     T = np.linspace(300, 3000, 1001)
    for n, thermo in enumerate(thermos):
        H = np.zeros(len(T))
        S = np.zeros(len(T))
        Cp = np.zeros(len(T))
        for i in range(0, len(T)):
            H[i] = thermo.get_enthalpy(T[i]) / 1000.0
            S[i] = thermo.get_entropy(T[i]) / 1000.0
            Cp[i] = thermo.get_heat_capacity(T[i]) / 1000.0
        ax[0].plot(T, H, linestyle=linestyles[n % len(linestyles)])
        ax[1].plot(T, S, linestyle=linestyles[n % len(linestyles)])
        ax[2].plot(T, Cp, linestyle=linestyles[n % len(linestyles)])
    ax[0].legend(labels)
    ax[1].legend(labels)
    ax[2].legend(labels)
    ax[2].yaxis.get_major_formatter().set_useOffset(False)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


def plot_gibbs(thermos, labels=None):
    if type(thermos) != list:
        thermos = [thermos]
    if labels is None:
        labels = ['' for t in thermos]
    linestyles = ['solid', 'solid', 'dashed', 'dashed']

    T = np.linspace(300, 2000, 1001)
#     T = np.linspace(300, 3000, 1001)
    for n, thermo in enumerate(thermos):
        G = np.zeros(len(T))
        for i in range(0, len(T)):
            G[i] = thermo.get_free_energy(T[i]) / 4184
        plt.plot(T, G, linestyle=linestyles[n % len(linestyles)])
       
    ax = plt.gca()
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('G (kcal / mol)')
    ax.set_title('Gibbs Energy vs. Temperature')
    ax.legend(labels)
   
    ax.yaxis.get_major_formatter().set_useOffset(False)
    plt.subplots_adjust(wspace=0.25)
    plt.show()


def get_i_thing(thing, thing_list):
    for i in range(len(thing_list)):
        if thing.is_isomorphic(thing_list[i]):
            return i
    return -1


# +
# Load my thermo calculations and a copy of the RMG surfaceThermoPt111 library
thermo_db = rmgpy.data.thermo.ThermoDatabase()
my_thermo_lib = '/home/moon/surface/surface_thermo/thermo/Pt_thermodata_adsorbates'
thermo_db.load_libraries(my_thermo_lib)

mylibname = 'Pt_thermodata_adsorbates'
assert mylibname in thermo_db.libraries
# -

rmg_items = [thermo_db.libraries['surfaceThermoPt111'].entries[x].item for x in thermo_db.libraries['surfaceThermoPt111'].entries]

# # Fe results

entry.item

# Plot thermo
mylibname = 'Fe_thermodata_adsorbates_gemnet'
# for entry in thermo_db.libraries[mylibname].entries:
for entry_name in thermo_db.libraries[mylibname].entries:
    
    entry = thermo_db.libraries[mylibname].entries[entry_name]
    print(entry.label)
    display(entry.item)
    # also grab the RMG version of the molecule
    rmg_idx = get_i_thing(entry.item, rmg_items)
    rmg_data = [thermo_db.libraries['surfaceThermoFe110'].entries[x].data for x in thermo_db.libraries['surfaceThermoFe110'].entries][rmg_idx]
    
    
    plot_thermos([entry.data, rmg_data], ['Calculation', 'RMG'])
    plot_gibbs([entry.data, rmg_data], ['Calculation', 'RMG'])





thermo_db.libraries

# # DFT results

# Plot thermo
mylibname = 'Pt111_gemnet'
# for entry in thermo_db.libraries[mylibname].entries:
for entry_name in thermo_db.libraries[mylibname].entries:
    
    entry = thermo_db.libraries[mylibname].entries[entry_name]
    print(entry.label)
    display(entry.item)
    # also grab the RMG version of the molecule
    rmg_idx = get_i_thing(entry.item, rmg_items)
    rmg_data = [thermo_db.libraries['surfaceThermoPt111'].entries[x].data for x in thermo_db.libraries['surfaceThermoPt111'].entries][rmg_idx]
    
    
    plot_thermos([entry.data, rmg_data], ['Calculation', 'RMG'])
    plot_gibbs([entry.data, rmg_data], ['Calculation', 'RMG'])


# # Equiformer Results

# Plot thermo
mylibname = 'Pt_thermodata_adsorbates'
# for entry in thermo_db.libraries[mylibname].entries:
for entry_name in thermo_db.libraries[mylibname].entries:
    
    entry = thermo_db.libraries[mylibname].entries[entry_name]
    print(entry.label)
    display(entry.item)
    # also grab the RMG version of the molecule
    rmg_idx = get_i_thing(entry.item, rmg_items)
    rmg_data = [thermo_db.libraries['surfaceThermoPt111'].entries[x].data for x in thermo_db.libraries['surfaceThermoPt111'].entries][rmg_idx]
    
    
    plot_thermos([entry.data, rmg_data], ['Calculation', 'RMG'])
    plot_gibbs([entry.data, rmg_data], ['Calculation', 'RMG'])


thermo_db.libraries



# read the results again
kinetics_lib = os.path.join('harris_kinetics')
ark_kinetics_database = rmgpy.data.kinetics.KineticsDatabase()
ark_kinetics_database.load_libraries(kinetics_lib)
# print(ark_kinetics_database.libraries)
print(f'{len(ark_kinetics_database.libraries["kinetics"].entries)} entries loaded')




