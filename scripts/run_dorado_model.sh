#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 1	
#SBATCH --mem 5G				# memory pool for all cores
#SBATCH --time=0-05:00:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

# Note downloading models requires run on software node as requires internet access

# Set file paths
simplex_model=../dorado_models/simplex/rna004_130bps_sup@v3.0.1
modification_model=../dorado_models/modification/rna004_130bps_sup@v3.0.1_m6A_DRACH@v1
pod5_folder=/ei/data/reads/PDEV-88/2024_10_24_PSEQ2766/R9_Wheat_root/20241024_1559_MN32565_FBA19032_d1057abb/pod5/
reference=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
output_bam_simplex=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/simplex/20241024_1559_MN32565_FBA19032_d1057abb.bam
output_summary_simplex=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/simplex/20241024_1559_MN32565_FBA19032_d1057abb.summary.txt
output_bam_modification=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.bam
output_summary_modification=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.summary.txt

# Source dorado
source package b1490e9c-4180-47c9-85c0-bb300845b7be # dorado - 0.7.2 

# List models available for download 
dorado download --list

# In the end I had issues with the software node so downloaded the model locally and then moved it to the HPC 
# wget https://cdn.oxfordnanoportal.com/software/analysis/dorado/rna004_130bps_sup@v3.0.1.zip
# wget https://cdn.oxfordnanoportal.com/software/analysis/dorado/rna004_130bps_sup@v3.0.1_m6A_DRACH@v1.zip 
# Then move on to the HPC and unzip. 

# Rebase call and align - Simplex Model 
dorado basecaller \
      "$simplex_model" \
      "$pod5_folder" \
      --reference "$reference" \
      --estimate-poly-a \
      > "$output_bam_simplex"

# Summary 
dorado summary $output_bam_simplex > $output_summary_simplex 

# Rebase call and align - Modification model 
dorado basecaller \
      "$simplex_model" \
      --modified-bases-models "$modification_model" \
      "$pod5_folder" \
      --reference "$reference" \
      > "$output_bam_modification"

# Summary 
dorado summary $output_bam_modification > $output_summary_modification 

