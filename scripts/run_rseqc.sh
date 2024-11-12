#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 16	
#SBATCH --mem 150G				# memory pool for all cores
#SBATCH --time=0-008:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address


# File Paths

iso_seq_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.bam
iso_seq_bam_sorted=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/CAD_R.PacBio_sorted.bam
direct_rna_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.bam
direct_rna_bam_sorted=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_sorted.bam

reference_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3.bed


# Sort and Index Bam files 
#source package c92263ec-95e5-43eb-a527-8f1496d56f1a # Samtools 1.18
#samtools sort -@ 16 -o "$iso_seq_bam_sorted" "$iso_seq_bam"
#samtools index -c -@ 16 "$iso_seq_bam_sorted"


#samtools sort -@ 16 -o "$direct_rna_bam_sorted" "$direct_rna_bam"
#samtools index -c -@ 16 "$direct_rna_bam_sorted"

source ~/.bashrc
conda activate /hpc-home/mahony/miniconda3
geneBody_coverage.py -h
geneBody_coverage.py --version
geneBody_coverage.py -r "$reference_bed" -i  "$iso_seq_bam" -o out


