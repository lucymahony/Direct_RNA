#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 1
#SBATCH --mem 50G				# memory pool for all cores
#SBATCH --time=0-01:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

# Plot the distibution of polyA tails as reported from the simplex model. 
#Simplex
source ~/.bashrc
mamba activate /hpc-home/mahony/miniforge3
conda run -n miniconda_dna python measure_polyA_counts.py 
conda run -n miniconda_dna python plot_polyA_counts.py 

