#!/bin/bash
#SBATCH --job-name=sobol_sample
#SBATCH --partition=sharing,west,short
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=4
#SBATCH --array=0-512%20


python /projects/westgroup/harris.se/surface_thermo/diffusion/sample_site.py $SLURM_ARRAY_TASK_ID