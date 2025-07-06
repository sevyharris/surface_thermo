#!/bin/bash

# define adsorbates
adsorbates=("O" "H" "N" "C" "H2" "N2" "H2O" "O2" "NH3")

for metal in 'Pt' 'Rh' 'Ni'; do
    for facet in '111' ; do
        for adsorbate in "${adsorbates[@]}"; do
            for site in 'ontop' 'bridge' 'fcc' 'hcp'; do
                echo "Running relaxation for $metal $facet + $adsorbate at site $site"
                python relax_system.py --metal "$metal" --facet "$facet" --adsorbate "$adsorbate" --site "$site"
                echo "Finished relaxation for $metal $facet + $adsorbate at site $site"
            done
        done
    done
done

for metal in 'Cr' 'Fe'; do
    for facet in '110' ; do
        for adsorbate in "${adsorbates[@]}"; do
            for site in 'ontop' 'hollow' 'shortbridge' 'longbridge'; do
                echo "Running relaxation for $metal $facet + $adsorbate at site $site"
                python relax_system.py --metal "$metal" --facet "$facet" --adsorbate "$adsorbate" --site "$site"
                echo "Finished relaxation for $metal $facet + $adsorbate at site $site"
            done
        done
    done
done