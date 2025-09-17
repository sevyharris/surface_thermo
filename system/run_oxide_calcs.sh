#!/bin/bash

# define adsorbates
# adsorbates=("H" "O" "N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "NH2" "CH3" "CH2" "CH" "CO" "CO2")
adsorbates=("N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "NH2" "CH3" "CH2" "CH" "CO" "CO2")
# adsorbates=("O2" "NH3" "NH" "OH" "NO" "ON")
# adsorbates=("NH")
# adsorbates=("NH2")


for slabname in 'Cr2O3_z'; do
    for adsorbate in "${adsorbates[@]}"; do
        for site in '0' '1' '2' '3' '4'; do
            for rotation in 0 90 180; do
                echo "Running relaxation for $slabname + $adsorbate at site $site with rotation $rotation"

                python relax_custom_system.py --slabname "$slabname" --adsorbate "$adsorbate" --rotate "$rotation" --site "$site"

                # python relax_system.py --metal "$metal" --facet "$facet" --adsorbate "$adsorbate" --site "$site" --rotate "$rotation"
                echo "Finished relaxation for $slabname + $adsorbate at site $site with rotation $rotation"
            done
        done
    done

done
