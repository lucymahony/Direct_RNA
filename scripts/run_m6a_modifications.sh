#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 1
#SBATCH --mem 70G				# memory pool for all cores
#SBATCH --time=0-06:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

#Simplex
source ~/.bashrc
mamba activate /hpc-home/mahony/miniforge3
conda run -n miniconda_dna python m6a_modifications.py

