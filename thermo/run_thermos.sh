#!/bin/bash
# adsorbates=("O" "H" "N" "C")
adsorbates=("O" "H" "N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "NH2")


for adsorbate in "${adsorbates[@]}"; do

    python run_thermo_calc_gemnet_oxide.py --slabname=Cr2O3_z --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet_oxide.py --slabname=Fe2O3_z --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Pt --crystal_structure=fcc --facet=111 --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Pt --crystal_structure=fcc --facet=100 --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Cr --crystal_structure=bcc --facet=100 --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Cr --crystal_structure=bcc --facet=110 --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Fe --crystal_structure=fcc --facet=111 --adsorbate="$adsorbate"
    python run_thermo_calc_gemnet.py --metal=Fe --crystal_structure=bcc --facet=100 --adsorbate="$adsorbate"

done


python energies2NASA.py --system_name=Cr2O3_z
python energies2NASA.py --system_name=Fe2O3_z
python energies2NASA.py --system_name=Pt_fcc111
python energies2NASA.py --system_name=Pt_fcc100
python energies2NASA.py --system_name=Cr_bcc100
python energies2NASA.py --system_name=Cr_bcc110
python energies2NASA.py --system_name=Fe_bcc100
python energies2NASA.py --system_name=Fe_fcc111
