#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 16	
#SBATCH --mem 150G				# memory pool for all cores
#SBATCH --time=0-008:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address


# For comparison to Direct RNA we are mapping some ISO seq data from Cadenza root too.

# File Paths
iso_seq_transcripts=/ei/references/internal/projects/10wheatgenomes/Transcripts/PacBio/CAD_R.PacBio.fasta.gz
reference=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
reference_index=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.mmi
output_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.bam
output_sam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.sam
# Source minimap2 and samtools
source package 222eac79-310f-4d4b-8e1c-0cece4150333 #Minimap2 2.24-41122
source package c92263ec-95e5-43eb-a527-8f1496d56f1a # Samtools 1.18

#minimap2 -x -d "$reference_index" "$reference"

minimap2 -t 16 -x splice:hq --split-prefix tmp -G 50000 -uf -t 30 -L --eqx -2 --secondary=no "$reference" "$iso_seq_transcripts" -a > "$output_sam"

#samtools view "$output_sam" -@ 16 -b | samtools sort -m 4G -@ 16 > "$output_bam" 
