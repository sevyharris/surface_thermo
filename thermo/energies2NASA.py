import numpy as np
import scipy
import glob
import pylab
import yaml
import sys
import os
import collections
import matplotlib
import matplotlib.pyplot  as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import NullFormatter, MaxNLocator, LogLocator
import rmgpy.data.thermo
import rmgpy.species
import argparse


# get adsorbate label from input
parser = argparse.ArgumentParser(description='Run thermo calculation for all adsorbates on a slab')
parser.add_argument('--system_name', type=str, required=True, help='System name, like Pt_fcc111')

args = parser.parse_args()
system_name = args.system_name



# define constants
R = 8.3144621E-3  # ideal Gas constant in kJ/mol-K
kB = 1.38065e-23  # Boltzmann constant in J/K
h = 6.62607e-34 #Planck constant in J*s
c = 2.99792458e8 #speed of light in m/s
amu = 1.6605e-27 #atomic mass unit in kg
Avogadro = 6.0221E23 #mole^-1
GHz_to_Hz = 1.0E9 #convert rotational constants from GHz to Hz
invcm_to_invm = 1.0E2 #convert cm^-1 to m^-1, for frequencies
P_ref = 1.0E5 #reference pressure, 1 bar = 1E5 Pascal
hartree_to_kcalpermole = 627.5095 #convert hartree/molecule to kcal/mol
hartree_to_kJpermole = 2627.25677 #convert hartree/molecule to kJ/mol
eV_to_kJpermole = 96.485 #convert eV/molecule to kJ/mol
# T_switch = 1000.0 #K, switching temperature in NASA polynomial. Default. can overwrite.
# site_occupation_number = 1 #number of sites occupied by adsorbate
# unit_cell_area = 62.10456e-20/9.0 #m2 - using surface area per binding site (nine binding sites per cell)
cutoff_frequency = 50.0 #cm^-1
# twoD_gas = False

# # declare a class for molecules
class Molecule:
    def __init__(self):
        self.T_switch = 1000.0 #K, switching temperature in NASA polynomial. Default. can overwrite.
        self.site_occupation_number = 1 #number of sites occupied by adsorbate
        self.unit_cell_area = 62.10456e-20/9.0 #m2 - using surface area per binding site (nine binding sites per cell)
        self.twoD_gas = False
        self.frequencies_units = 'cm-1' #default units for frequencies
                

# create the array of temperatures in 10 degree increments
temperature = [298.15] #NOTE 298.15 must be first for the NASA polynomial routine to work!
T_low = 300.0
T_high = 2000.0
dT = 10.0 #temperature increment
temperature = np.append(temperature, np.arange(T_low, T_high+dT, dT) )


# HERE BEGINS THE LONG LIST OF SUBROUTINES
#-------------------------------------------------------------------------
# subroutine for the translational mode
def get_translation_thermo(molecule, temperature):
    # unpack the constants (not essential, but makes it easier to read)

    area = molecule.unit_cell_area
    sites = molecule.site_occupation_number
    m = molecule.adsorbate_mass

    #initialize the arrays for the partition function, entropy, enthalpy,
    #and heat capacity.
    Q_trans  = np.ones(len(temperature)) 
    S_trans  = np.zeros(len(temperature))
    dH_trans  = np.zeros(len(temperature))
    Cp_trans  = np.zeros(len(temperature))

    if molecule.twoD_gas:
        print("switching to 2D-gas for 2 lowest modes for %s"%molecule.name)
        # cycle through each temperature
        for (i,T) in enumerate(temperature):
            # partition function is: (2*pi*mass*kB*T/h**2)^(2/2) * area
            if (1==0): #3D gas, really here just for inspiration
                V = kB*T/P_ref
                Q_trans[i] = (2 * np.pi * m * amu*kB*T/h**2)**(1.5) * V
                S_trans[i] = R * (2.5 + np.log( Q_trans[i] )) #
                Cp_trans[i] = R * 2.5 #NOTE: Cp = Cv + R
                dH_trans[i] = R * 2.5 * T      
            else: #surface
                if (1==0): #Campbell + Arnadottir
                    V = kB*T/P_ref
                    Q_trans[i] = (2 * np.pi * m * amu*kB*T/h**2)**(1.0) *V**0.66667
                    S_trans[i] = R * (2.0 + np.log( Q_trans[i] ))
                    Cp_trans[i] = R * 1.66667 #NOTE: Cp = Cv + 2/3R
                    dH_trans[i] = R * 1.66667 * T            
        
                else: #area is not a function of temperature
                    if sites > 0:
                        Q_trans[i] = (2*np.pi*m*amu*kB*T/h**2) * area * sites
                    else:
                        Q_trans[i] = (2*np.pi*m*amu*kB*T/h**2) * area
                    S_trans[i] = R * (2.0 + np.log( Q_trans[i] ))
                    Cp_trans[i] = R * 1.0 #NOTE: Cp = Cv 
                    dH_trans[i] = R * 1.0 * T            

    # add the results to the thermo object
    molecule.Q_trans = Q_trans
    molecule.S_trans = S_trans
    molecule.dH_trans = dH_trans
    molecule.Cp_trans = Cp_trans 
    

    return


# subroutine for the vibrational mode
def get_vibrational_thermo(molecule,temperature):
    units = 1.0
    units *= h * c / kB * invcm_to_invm # K * cm
    
    #initialize the arrays for the partition function, entropy, enthalpy,
    #and heat capacity.
    Q_vib  = np.ones(len(temperature))
    S_vib  = np.zeros(len(temperature))
    dH_vib  = np.zeros(len(temperature))
    Cv_vib  = np.zeros(len(temperature))
    
    for (t,temp) in enumerate(temperature):
        for (n,nu) in enumerate(molecule.frequencies):
            if molecule.twoD_gas==True and n <= 1: #skip the first two if we do 2D gas
                #do nothing!
                Q_vib[t] *= 1.0
                S_vib[t] += 0.0
                dH_vib[t] += 0.0
                Cv_vib[t] += 0.0
            else:
                # if nu < 50:  # mine
                #     nu = 50.0 #avoid numerical issues with very low frequencies
                x = nu * units / temp #cm^-1 * K cm / K = dimensionless
                Q_vib[t]  *= 1.0 / (1.0 - np.exp( - x) )
                S_vib[t]  += -np.log( 1.0 - np.exp( - x ) ) + x * np.exp( - x) / (1.0 - np.exp( - x) ) 
                dH_vib[t] += x * np.exp( - x) / (1.0 - np.exp( - x) ) 
                Cv_vib[t] += x**2.0 * np.exp( - x) / (1.0 - np.exp( - x) )**2.0
        S_vib[t]  *= R
        dH_vib[t] *= R * temp
        Cv_vib[t] *= R

    # add the results to the thermo object
    molecule.Q_vib = Q_vib
    molecule.S_vib = S_vib
    molecule.dH_vib = dH_vib
    molecule.Cv_vib = Cv_vib #NOTE: the correction from Cv to Cp is handled in the translation partition function.
                             #if the molecule is tightly bound and thus the 2D-gas is not used, 
                             #then we assume that Cp=Cv for the adsorbate.

    return

#-------------------------------------------------------------------------
#create the main thermo function that calls the individual modes
def thermo(molecule, temperature):

    # call the subroutine for the vibrational partition function
    get_translation_thermo(molecule,temperature)
    get_vibrational_thermo(molecule,temperature)

    
    #now compute the correction to the heat of formation as you go from 0 to 298 K
    h_correction = 4.234 #kJ/mol. enthalpy_H(298) - enthalpy_H(0)
    c_correction = 1.051 #kJ/mol. enthalpy_C(298) - enthalpy_C(0)
    n_correction = 4.335 #kJ/mol. enthalpy_N(298) - enthalpy_N(0)
    o_correction = 4.340 #kJ/mol. enthalpy_O(298) - enthalpy_O(0)
    
    molecule.heat_of_formation_correction = 0.0
    molecule.heat_of_formation_correction += molecule.composition['H'] * h_correction
    molecule.heat_of_formation_correction += molecule.composition['C'] * c_correction    
    molecule.heat_of_formation_correction += molecule.composition['N'] * n_correction
    molecule.heat_of_formation_correction += molecule.composition['O'] * o_correction        
    
    # note that the partition function is the production of the individual terms,
    # whereas the thermodynamic properties are additive
    molecule.Q = molecule.Q_trans * molecule.Q_vib 
    molecule.S = molecule.S_trans + molecule.S_vib 
    molecule.dH = molecule.dH_trans + molecule.dH_vib 
    molecule.Cp = molecule.Cp_trans + molecule.Cv_vib # see comments in each section regarding Cp vs Cv
    molecule.heat_of_formation_298K = molecule.heat_of_formation_0K + molecule.dH[0] - molecule.heat_of_formation_correction
    molecule.H = molecule.heat_of_formation_298K + molecule.dH - molecule.dH[0]
    
    print(molecule.heat_of_formation_298K)
    print(molecule.H[0])
    #This writes H_298, S_298 and appropriate indices of Cp to file (preparation for computing adsorption corrections)
    g = open("Pt_thermodata_adsorbates.py",'a+')
    g.write('[' + str(molecule.name) + ', Cpdata:, ' +  str(molecule.Cp[np.where(temperature==300)]*239.0057)[1:-1] + ', ' + str(molecule.Cp[np.where(temperature==400)]*239.0057)[1:-1] + ', '+ str(molecule.Cp[np.where(temperature==500)]*239.0057)[1:-1] + ', ' + str(molecule.Cp[np.where(temperature==600)]*239.0057)[1:-1] + ', ' + str(molecule.Cp[np.where(temperature==800)]*239.0057)[1:-1] + ', ' + str(molecule.Cp[np.where(temperature==1000)]*239.0057)[1:-1] + ', ' + str(molecule.Cp[np.where(temperature==1500)]*239.0057)[1:-1] + ', ' + ",'cal/(mol*K)', H298, " + str(molecule.H[0]*0.2390057) + ", 'kcal/mol', S298, " + str(molecule.S[0]*239.0057) + ", 'cal/(mol*K)']")
    g.write('\n')
    g.close()
    
    # now that we've computed the thermo properties, go ahead and fit them to a NASA polynomial
    fit_NASA(temperature, molecule)
    format_output(molecule)
    return

#-------------------------------------------------------------------------
#compute thermo properties from nasa polynomials
def get_thermo_from_NASA(temperature, molecule):
    
    a_low = molecule.a_low
    a_high = molecule.a_high
    T_switch = molecule.T_switch
    
    i_switch = -1
    for i in range(len(temperature)):
        if temperature[i]==T_switch:
            i_switch = i
    
    cp_fit = np.zeros(len(temperature))
    h_fit = np.zeros(len(temperature))
    s_fit = np.zeros(len(temperature))
    for (i,temp) in enumerate(temperature):
        if temp <= T_switch:
            cp_fit[i] = a_low[0] + a_low[1]*temp + a_low[2]*temp**2.0  + a_low[3]*temp**3.0  + a_low[4]*temp**4.0
            h_fit[i] = a_low[0]*temp + a_low[1]/2.0*temp**2.0 + a_low[2]/3.0*temp**3.0  + a_low[3]/4.0*temp**4.0  + a_low[4]/5.0*temp**5.0 + a_low[5]
            s_fit[i] = a_low[0]*np.log(temp) + a_low[1]*temp + a_low[2]/2.0*temp**2.0  + a_low[3]/3.0*temp**3.0  + a_low[4]/4.0*temp**4.0 + a_low[6]
        else:
            cp_fit[i] = a_high[0] + a_high[1]*temp + a_high[2]*temp**2.0  + a_high[3]*temp**3.0  + a_high[4]*temp**4.0
            h_fit[i] = a_high[0]*temp + a_high[1]/2.0*temp**2.0 + a_high[2]/3.0*temp**3.0  + a_high[3]/4.0*temp**4.0  + a_high[4]/5.0*temp**5.0 + a_high[5]
            s_fit[i] = a_high[0]*np.log(temp) + a_high[1]*temp + a_high[2]/2.0*temp**2.0  + a_high[3]/3.0*temp**3.0  + a_high[4]/4.0*temp**4.0 + a_high[6]

    cp_fit *= R        
    h_fit *= R  
    s_fit *= R  
    
    molecule.Cp_fit = cp_fit
    molecule.H_fit = h_fit
    molecule.S_fit = s_fit
    return 


#-------------------------------------------------------------------------
#fit nasa coefficients
def fit_NASA(temperature, molecule):
    
    heat_capacity = molecule.Cp
    reference_enthalpy = molecule.H[0]
    reference_entropy = molecule.S[0]
    T_switch = molecule.T_switch
    
    i_switch = -1
    for i in range(len(temperature)):
        if temperature[i]==T_switch:
            i_switch = i
    if i_switch==-1:
        print("We have a problem! Cannot find switching temperature")
        
    
    #start by creating the independent variable matrix for the low-temperature fit
    YT = np.array( [ np.ones(len(temperature[:i_switch+1])), temperature[:i_switch+1], temperature[:i_switch+1]**2.0, temperature[:i_switch+1]**3.0, temperature[:i_switch+1]**4.0 ],dtype=np.float64 ) #this is transpose of our Y
    Y = YT.transpose() #this is the desired Y

    b = heat_capacity[:i_switch+1] / R  
    a_low = np.linalg.lstsq(Y, b, rcond=None)[0]

    T_ref = 298.15
    #now determine the enthalpy coefficient for the low-T region
    subtract = a_low[0] + a_low[1]/2.0*T_ref + a_low[2]/3.0*T_ref**2.0 + a_low[3]/4.0*T_ref**3.0  + a_low[4]/5.0*T_ref**4.0
    a_low = np.append(a_low, reference_enthalpy / R - subtract * T_ref)
    #now determine the entropy coefficient for the low-T region
    subtract = a_low[0] * np.log(T_ref) + a_low[1]*T_ref     + a_low[2]/2.0*T_ref**2.0  + a_low[3]/3.0*T_ref**3.0  + a_low[4]/4.0*T_ref**4.0
    a_low = np.append(a_low, reference_entropy / R - subtract )

    #
    # NOW SWITCH TO HIGH-TEMPERATURE REGIME!
    #
    T_ref = T_switch
    #compute the heat capacity, enthalpy, and entropy at the switching point
    Cp_switch = a_low[0] + a_low[1]*T_ref + a_low[2]*T_ref**2.0  + a_low[3]*T_ref**3.0  + a_low[4]*T_ref**4.0
    H_switch = a_low[0]*T_ref + a_low[1]/2.0*T_ref**2.0 + a_low[2]/3.0*T_ref**3.0  + a_low[3]/4.0*T_ref**4.0  + a_low[4]/5.0*T_ref**5.0 + a_low[5]
    S_switch = a_low[0]*np.log(T_ref) + a_low[1]*T_ref + a_low[2]/2.0*T_ref**2.0  + a_low[3]/3.0*T_ref**3.0  + a_low[4]/4.0*T_ref**4.0 + a_low[6]
    
    #now repeat the process for the high-temperature regime
    a_high = [0.0]
    YT = np.array( [ temperature[i_switch:], temperature[i_switch:]**2.0, temperature[i_switch:]**3.0, temperature[i_switch:]**4.0 ],dtype=np.float64 ) #this is transpose of our Y
    Y = YT.transpose() #this is the desired Y

    b = heat_capacity[i_switch:] / R - Cp_switch
    a_high = np.append(a_high, np.linalg.lstsq(Y, b, rcond=None)[0])
    a_high[0] = Cp_switch - (a_high[0] + a_high[1]*T_switch + a_high[2]*T_switch**2.0  + a_high[3]*T_switch**3.0  + a_high[4]*T_switch**4.0)
    
    a_high = np.append(a_high, H_switch - (a_high[0] + a_high[1]/2.0*T_ref + a_high[2]/3.0*T_ref**2.0  + a_high[3]/4.0*T_ref**3.0  + a_high[4]/5.0*T_ref**4.0)*T_ref )
    a_high = np.append(a_high, S_switch - (a_high[0]*np.log(T_ref) + a_high[1]*T_ref + a_high[2]/2.0*T_ref**2.0  + a_high[3]/3.0*T_ref**3.0  + a_high[4]/4.0*T_ref**4.0) )

    #Check to see if there is a discontinuity
    if (1==0):
        print("\ncheck for discontinuities:")
        cp_low_Tswitch = a_low[0] + a_low[1]*T_switch + a_low[2]*T_switch**2.0  + a_low[3]*T_switch**3.0  + a_low[4]*T_switch**4.0
        cp_high_Tswitch = a_high[0] + a_high[1]*T_switch + a_high[2]*T_switch**2.0  + a_high[3]*T_switch**3.0  + a_high[4]*T_switch**4.0
        H_low_Tswitch = a_low[0]*T_switch + a_low[1]/2.0*T_switch**2.0 + a_low[2]/3.0*T_switch**3.0  + a_low[3]/4.0*T_switch**4.0  + a_low[4]/5.0*T_switch**5.0 + a_low[5]
        H_high_Tswitch = a_high[0]*T_switch + a_high[1]/2.0*T_switch**2.0 + a_high[2]/3.0*T_switch**3.0  + a_high[3]/4.0*T_switch**4.0  + a_high[4]/5.0*T_switch**5.0 + a_high[5]
        S_low_Tswitch = a_low[0]*np.log(T_switch) + a_low[1]*T_switch + a_low[2]/2.0*T_switch**2.0  + a_low[3]/3.0*T_switch**3.0  + a_low[4]/4.0*T_switch**4.0 + a_low[6]
        S_high_Tswitch = a_high[0]*np.log(T_switch) + a_high[1]*T_switch + a_high[2]/2.0*T_switch**2.0  + a_high[3]/3.0*T_switch**3.0  + a_high[4]/4.0*T_switch**4.0 + a_high[6]    

        print("discontinuity at T_switch for Cp/R is %.4F"%(cp_low_Tswitch - cp_high_Tswitch))
        print("discontinuity at T_switch for H/R is %.4F"%(H_low_Tswitch - H_high_Tswitch))   
        print("discontinuity at T_switch for S/R is %.4F"%(S_low_Tswitch - S_high_Tswitch))        
    
    #line = '\n\t !cut and paste this value into the cti file!\n'
    line = '\tthermo = (\n'
    line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n"%(300.0, 1000.0, a_low[0], a_low[1], a_low[2], a_low[3], a_low[4], a_low[5], a_low[6])
    line += "\t\tNASA( [%.1F, %.1F], [%.8E, %.8E,\n \t\t %.8E, %.8E, %.8E,\n \t\t %.8E, %.8E]), \n"%(1000.0, max(temperature), a_high[0], a_high[1], a_high[2], a_high[3], a_high[4], a_high[5], a_high[6])
    line += "\t\t ),\n"

    molecule.thermo_lines = line

    molecule.a_low = a_low
    molecule.a_high = a_high
    
    return 


#-------------------------------------------------------------------------
#compare NASA fits to computed fits
def compare_NASA_to_thermo(temperature, molecule):
    
    fig = pylab.figure(dpi=300,figsize=(12,4))
    gs = gridspec.GridSpec(1, 3)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])

    if (1==1): #use this to plot the absolute curves
        ax0.plot(temperature, molecule.Cp, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
        ax0.plot(temperature, molecule.Cp_fit, 'b', linewidth=2)
        ax1.semilogy(temperature, molecule.H - molecule.heat_of_formation_298K, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
        ax1.semilogy(temperature, molecule.H_fit - molecule.heat_of_formation_298K, 'b', linewidth=2)
        ax2.semilogy(temperature, molecule.S, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
        ax2.semilogy(temperature, molecule.S_fit, 'b', linewidth=2)
        ax0.set_ylim(min(molecule.Cp_fit)*0.9, max(molecule.Cp_fit)*1.025)
        ax1.set_ylim(min(molecule.H - molecule.heat_of_formation_298K)*0.9, max(molecule.H - molecule.heat_of_formation_298K)*1.025)
        ax2.set_ylim(min(molecule.S_fit)*0.9, max(molecule.S_fit)*1.025)
        ax1.yaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
        ax2.yaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
        
    else: #use this one to plot the percent change    
        ax0.plot(temperature, 1.0 - molecule.Cp/molecule.Cp_fit, 'b', linewidth=2)
        ax1.plot(temperature, 1.0 - molecule.H/molecule.H_fit, 'b', linewidth=2)
        ax2.plot(temperature, 1.0 - molecule.S/molecule.S_fit, 'b', linewidth=2)
        ax0.set_ylim(-5E-3, 5E-3)
        ax1.set_ylim(-5E-3, 5E-3)
        ax2.set_ylim(-5E-3, 5E-3)
        ax1.yaxis.set_major_locator(MaxNLocator(4))
        ax2.yaxis.set_major_locator(MaxNLocator(4))
        
    # now make it look better
    ax0.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax0.xaxis.set_major_locator(MaxNLocator(4))
    ax0.yaxis.set_major_locator(MaxNLocator(4))
    ax0.tick_params(axis='both', which='major', labelsize=12)
    ax0.set_title("heat capacity")
    ax0.set_xlabel("temperature [K]", fontsize=12)

    ax1.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.set_title("change in enthalpy")
    ax1.set_xlabel("temperature [K]", fontsize=12)

    ax2.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax2.xaxis.set_major_locator(MaxNLocator(4))
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.set_title("entropy")
    ax2.set_xlabel("temperature [K]", fontsize=12)
    
    return


#-------------------------------------------------------------------------
#compare NASA fits to computed fits
def compare_Cantera_to_thermo(temperature, molecule):
    
    fig = pylab.figure(dpi=300,figsize=(12,4))
    gs = gridspec.GridSpec(1, 3)
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1])
    ax2 = plt.subplot(gs[2])

    ax0.plot(temperature, molecule.Cp, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
    ax0.plot(temperature, molecule.cantera_cp, 'b', linewidth=2)

    ax1.semilogy(temperature, molecule.H- molecule.heat_of_formation_298K, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
    ax1.semilogy(temperature, molecule.cantera_h- molecule.heat_of_formation_298K, 'b', linewidth=2)

    ax2.semilogy(temperature, molecule.S, marker='o', markeredgecolor='r',color='w',alpha=0.5,linestyle='None')
    ax2.semilogy(temperature, molecule.cantera_s, 'b', linewidth=2)

    # now make it look better
    ax0.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax0.set_ylim(min(molecule.Cp_fit)*0.9, max(molecule.Cp_fit)*1.025)
    ax0.xaxis.set_major_locator(MaxNLocator(4))
    ax0.yaxis.set_major_locator(MaxNLocator(4))
    ax0.tick_params(axis='both', which='major', labelsize=12)
    ax0.set_title("heat capacity")
    ax0.set_xlabel("temperature [K]", fontsize=12)

    ax1.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax1.set_ylim(min(molecule.H - molecule.heat_of_formation_298K)*0.9, max(molecule.H - molecule.heat_of_formation_298K)*1.025)
    ax1.xaxis.set_major_locator(MaxNLocator(4))
    ax1.yaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
    ax1.tick_params(axis='both', which='major', labelsize=12)
    ax1.set_title("change in enthalpy")
    ax1.set_xlabel("temperature [K]", fontsize=12)

    ax2.set_xlim(min(temperature)*0.95, max(temperature)*1.025)
    ax2.set_ylim(min(molecule.S_fit)*0.9, max(molecule.S_fit)*1.025)
    ax2.xaxis.set_major_locator(MaxNLocator(4))
    ax2.yaxis.set_major_locator(LogLocator(base=10.0, numticks=4))
    ax2.tick_params(axis='both', which='major', labelsize=12)
    ax2.set_title("entropy")
    ax2.set_xlabel("temperature [K]", fontsize=12)
    
    return



#-------------------------------------------------------------------------
#print(output in cti format)
def format_output(molecule):
    
    line = '\n'
    line += 'species(name = "%s",\n'%(molecule.name)
    line += '\tatoms = "'
    for element in molecule.composition:
        if molecule.composition[element]>0:
            line += " %s:%d"%(element, molecule.composition[element])
    line += '",\n'
    line += "\tsize = %d,\n"%(molecule.site_occupation_number)
    line += molecule.thermo_lines
    line += '    longDesc = u"""Calculated by x at x University using statistical mechanics (file: compute_NASA_for_Pt-adsorbates.ipynb). \n' 
    line += "                   Based on DFT calculations by x at x.\n"
    line += "            DFT binding energy: %.3F %s.\n" %(molecule.DFT_binding_energy, molecule.DFT_binding_energy_units.replace("'",""))
    # if molecule.site_occupation_number == 1:
    #     line += "            Linear scaling parameters: ref_adatom_%s = %.3F %s, psi = %.5F %s, gamma_%s(X) = %.3F."%(molecule.binding_atom1, molecule.ref_adatom_Eb1, molecule.ref_adatom_Eb1_units.replace("'",""), molecule.linear_scaling_psi, molecule.linear_scaling_psi_units.replace("'",""), molecule.binding_atom1, molecule.linear_scaling_gamma)
    # else:
    #     line += "            Linear scaling parameters: ref_adatom_%s1 = %.3F %s, ref_adatom_%s2 = %.3F %s, psi = %.5F %s, gamma_%s1(X) = %.3F, gamma_%s2(X) = %.3F."%(molecule.binding_atom1, molecule.ref_adatom_Eb1, molecule.ref_adatom_Eb1_units.replace("'",""), molecule.binding_atom2, molecule.ref_adatom_Eb2, molecule.ref_adatom_Eb2_units.replace("'",""), molecule.linear_scaling_psi, molecule.linear_scaling_psi_units.replace("'",""), molecule.binding_atom1, molecule.linear_scaling_gamma, molecule.binding_atom2, molecule.linear_scaling_gamma_B)
    if molecule.twoD_gas:
        line += '\n            The two lowest frequencies, %.1F and %.1F %s, where replaced by the 2D gas model.' %(molecule.frequencies[0], molecule.frequencies[1], molecule.frequencies_units.replace("'",""))
    line += '""",\n\t)\n'
    
    molecule.species_lines = line
    
    return


#-------------------------------------------------------------------------
#Define the input parser
def parse_input_file(input_file, molecule):

    # read yaml file
    with open(input_file, 'r') as f:
        input_data = yaml.load(f, Loader=yaml.FullLoader)

    
    if input_data['adjacency_list'] is None:
        print(f'Error: adjacency_list is None for {input_file}. Skipping this file.')
        return False
    molecule.adjacency_list = input_data['adjacency_list'].replace('Pt', 'X')


    vdw_translater = {
        'O': 'O',
        '[H][H].Pt': 'H',
    }

    # get the binding atom
    sp = rmgpy.species.Species().from_adjacency_list(molecule.adjacency_list)
    print(f'Parsing species: {sp.smiles}')
    for atom in sp.molecule[0].atoms:
        if atom.is_surface_site():
            bonds = sp.molecule[0].get_bonds(atom)

            # if it's not vdW
            if len(bonds) > 0:
                element1 = list(bonds.keys())[0].symbol
                print(f'binding atom is {element1}')
            else:
                # get the element with the biggest atomic number
                possible_atoms = [atom for atom in sp.molecule[0].atoms if not atom.is_surface_site()]
                element1 = possible_atoms[np.argmax([possible_atoms[i].number for i in range(len(possible_atoms))])].symbol

                # element1 = vdw_translater[sp.smiles]
                print(f'binding atom is {element1}')


    molecule.binding_atom1 = str(element1)
    molecule.name = input_data['name']

    molecule.composition = input_data['composition']
    N_adsorbate_atoms = 0
    for element in molecule.composition:
        if element != 'Pt' and element != 'X':
            N_adsorbate_atoms += molecule.composition[element]  

    molecule.site_occupation_number = input_data['sites']
    molecule.DFT_binding_energy = input_data['DFT_binding_energy'][0]
    molecule.DFT_binding_energy_units = input_data['DFT_binding_energy'][1]
    molecule.heat_of_formation_0K = input_data['heat_of_formation_0K'][0]
    molecule.heat_of_formation_0K_units = input_data['heat_of_formation_0K'][1]

    # if provided in eV, convert to kJ/mol
    if molecule.heat_of_formation_0K_units == 'eV':
        molecule.heat_of_formation_0K = molecule.heat_of_formation_0K * eV_to_kJpermole
        molecule.heat_of_formation_0K_units = 'kJ/mol'

    molecule.adsorbate_mass = input_data['adsorbate_mass'][0]
    molecule.adsorbate_mass_units = input_data['adsorbate_mass'][1]

    molecule.ref_adatom_Eb1 = input_data['linear_scaling_binding_atom'][0]
    molecule.ref_adatom_Eb1_units = input_data['linear_scaling_binding_atom'][1]

    molecule.linear_scaling_gamma = input_data['linear_scaling_gamma(X)']
    molecule.linear_scaling_gamma_B = input_data.get('linear_scaling_gamma(X)_B', None)
    molecule.linear_scaling_psi = input_data['linear_scaling_psi'][0]
    molecule.linear_scaling_psi_units = input_data['linear_scaling_psi'][1]
    molecule.frequencies = input_data['frequencies']
    # N_freq_computed = 3 * N_adsorbate_atoms
    if molecule.frequencies[1] < cutoff_frequency:
        #print("switching to 2D-gas for 2 lowest modes for %s"%name
        print(f'Switching to 2D-gas for 2 lowest modes for {molecule.name}')
        molecule.twoD_gas = True


    molecule.ref_adatom_Eb2 = input_data['linear_scaling_binding_atom_B'][0] if 'linear_scaling_binding_atom_B' in input_data else None
    molecule.ref_adatom_Eb2_units = input_data['linear_scaling_binding_atom_B'][1] if 'linear_scaling_binding_atom_B' in input_data else None
    molecule.binding_atom2 = input_data.get('binding_atom_B', None)
    
    return True

# -------------------------------------------------------------------------
# Main script to read species list and compute thermodynamic properties

# metal = 'Pt'
# crystal_structure = 'fcc'
# facet = '111'

# system_name = 'Pt_fcc111'
# system_name = 'Pt_fcc100'
# system_name = 'Fe_fcc111'
# system_name = 'Fe_bcc100'
# system_name = 'Cr_bcc100'
# system_name = 'Cr_bcc110'
# system_name = 'Cr2O3_z'
# system_name = 'Fe2O3_z'

metal = system_name.split('_')[0]
facet = system_name.split('_')[-1]

filenames = glob.glob(os.path.join(os.environ['SURFACE_THERMO_DIR'], f'results/thermo/{system_name}/{system_name}_*-ads.yaml'))
new_output = open('my_new_cti.txt', 'w')


name_line = '\n'
species_line = '\n'
counter = -1

# compile it all into a single database and a single library which I'll call harris_butane
# my_library_name = 'Pt_thermodata_adsorbates'
#my_library_name = f'{metal}_{crystal_structure}{facet}_gemnet'
my_library_name = f'{system_name}_gemnet'

output_file = f'{my_library_name}.py'
thermo_database = rmgpy.data.thermo.ThermoDatabase()
thermo_database.libraries[my_library_name] = rmgpy.data.thermo.ThermoLibrary()
thermo_database.libraries[my_library_name].name = my_library_name
thermo_database.libraries[my_library_name].label = my_library_name
thermo_database.libraries[my_library_name].entries = collections.OrderedDict()


# make X
sp = rmgpy.molecule.Molecule().from_adjacency_list("""
1 X u0 p0 c0
""")
index = 0
thermo_data = rmgpy.thermo.NASA(
    polynomials = [
        rmgpy.thermo.NASAPolynomial(coeffs=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
        rmgpy.thermo.NASAPolynomial(coeffs=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
    ],
    Tmin = (298.0, 'K'),
    Tmax = (2000.0, 'K'),
)
# make a thermo entry for the database
entry = rmgpy.data.base.Entry(
    index=index,
    label='X',
    item=sp,
    data=thermo_data,
)
thermo_database.libraries[my_library_name].entries[entry.label] = entry

HX = rmgpy.species.Species().from_adjacency_list("""
    1 X  u0 p0 c0 {2,S}
    2 H  u0 p0 c0 {1,S}
""")

for index, filename in enumerate(filenames):
    print(filename)

    counter += 1
    # filename = species.strip()

    test = Molecule()
    success = parse_input_file(filename, test)
    if not success:
        continue
    thermo(test, temperature)
    
    name_line += ' %s'%(test.name)
    if counter == 4:
        name_line +='\n'
        counter == -1
    species_line += test.species_lines

    get_thermo_from_NASA(temperature, test)
    # compare_NASA_to_thermo(temperature, test)

    # sp = rmgpy.species.Species().from_adjacency_list(test.adjacency_list)
    sp = rmgpy.molecule.Molecule().from_adjacency_list(test.adjacency_list)

    thermo_data = rmgpy.thermo.NASA(
        polynomials = [
            rmgpy.thermo.NASAPolynomial(coeffs=test.a_low, Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
            rmgpy.thermo.NASAPolynomial(coeffs=test.a_high, Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
        ],
        Tmin = (298.0, 'K'),
        Tmax = (2000.0, 'K'),
    )

    # make a thermo entry for the database
    entry = rmgpy.data.base.Entry(
        index=index + 1,
        label=test.name,
        item=sp,
        data=thermo_data,
        metal=metal,
        facet=facet,
    )

    thermo_database.libraries[my_library_name].entries[entry.label] = entry

    if sp.is_isomorphic(HX.molecule[0]):
        HX_thermo = rmgpy.thermo.NASA(
            polynomials = [
                rmgpy.thermo.NASAPolynomial(coeffs=[-1.96702988E+00, 1.67920714E-02,  -2.50314139E-05, 1.80485455E-08, -5.11491197E-12,  -3.21277026E+03, 7.68211257E+00], Tmin=(298.0,'K'), Tmax=(1000.0, 'K')),
                rmgpy.thermo.NASAPolynomial(coeffs=[2.71968546E+00, -1.07696656E-03,  2.00193294E-06, -1.12865983E-09, 2.11269165E-13,  -4.24701712E+03, -1.52793490E+01], Tmin=(1000.0,'K'), Tmax=(2000.0, 'K')),
            ],
            Tmin = (298.0, 'K'),
            Tmax = (2000.0, 'K'),
        )
        print(f'\tKatrin HX: {HX_thermo.get_free_energy(1000) / 4184.0} kcal/mol')
        print(f'\tYour HX: {thermo_data.get_free_energy(1000) / 4184.0} kcal/mol')
        print(thermo_data)


# the real writing
thermo_database.save_libraries(my_library_name)


name_line += '\n\n' 
new_output.write(name_line)
new_output.write(species_line)

new_output.close()



# def convert_energies_to_nasa_format(input_yaml):
#     """
#     Convert DFT energies to NASA format.
#     """

#     nasa_energies = []
#     for energy in energies:
#         # Convert each energy to NASA format
#         nasa_energy = {
#             "temperature": energy["temperature"],
#             "enthalpy": energy["enthalpy"] * 1.0,  # Convert to kJ/mol
#             "entropy": energy["entropy"] * 1.0,    # Convert to J/(mol*K)
#         }
#         nasa_energies.append(nasa_energy)
#     return nasa_energies
