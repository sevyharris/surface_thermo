#!/bin/bash
#SBATCH --job-name=run_neb
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --partition=short,west
#SBATCH --export=ALL

python /projects/westgroup/harris.se/surface_thermo/diffusion/run_neb.py $1