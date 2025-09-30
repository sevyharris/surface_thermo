#!/bin/bash
#SBATCH --job-name=run_gemnet_Pt_fcc111
#SBATCH --mem=20Gb
#SBATCH --time=12:00:00
#SBATCH --partition=short,west
#SBATCH --export=ALL


# define adsorbates
adsorbates=("O" "H" "N" "C" "H2" "N2" "H2O" "O2" "NH3" "NH" "OH" "NO" "NH2" "CH3" "CH2" "CH" "CO" "CO2")

for metal in 'Pt' ; do
   for facet in '111' ; do
       for adsorbate in "${adsorbates[@]}"; do
           for site in 'ontop' 'bridge' 'fcc' 'hcp'; do
               for rotation in 0 90 180; do
                   echo "Running relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
                   python relax_system.py --metal "$metal" --facet "$facet" --crystal_structure=fcc --adsorbate "$adsorbate" --site "$site" --rotate "$rotation"
                   echo "Finished relaxation for $metal $facet + $adsorbate at site $site with rotation $rotation"
               done
           done
       done
   done
done
