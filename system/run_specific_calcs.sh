#!/bin/bash

# define adsorbates
# adsorbates=("O" "H" "N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "ON")
adsorbates=("O" "H" "N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "NH2" "CH3" "CH2" "CH" "CO" "CO2")
# adsorbates=("O2" "NH3" "NH" "OH" "NO" "ON")
# adsorbates=("NH")
# adsorbates=("NH2")

# for metal in 'Pt' 'Rh' 'Ni'; do
# # for metal in 'Pt'; do
# for metal in 'Fe'; do
#    for facet in '111' ; do
#        for adsorbate in "${adsorbates[@]}"; do
#            for site in 'ontop' 'bridge' 'fcc' 'hcp'; do
#                for rotation in 0 90 180; do
#                    echo "Running relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
#                    python relax_system.py --metal "$metal" --facet "$facet" --adsorbate "$adsorbate" --site "$site" --rotate "$rotation"
#                    echo "Finished relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
#                done
#            done
#        done
#    done
# done

#  for metal in 'Cr' 'Fe'; do
 for metal in 'Fe' ; do
     for facet in '110' ; do
         for adsorbate in "${adsorbates[@]}"; do
             for site in 'ontop' 'hollow' 'shortbridge' 'longbridge'; do
                for rotation in 0 90 180; do
                 echo "Running relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
                 python relax_system.py --metal "$metal" --facet "$facet" --adsorbate "$adsorbate" --site "$site" --rotate "$rotation"
                 echo "Finished relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
                 done
             done
         done
     done
 done

#python relax_system.py --metal "Fe" --facet "111" --adsorbate "NH" --site ontop

