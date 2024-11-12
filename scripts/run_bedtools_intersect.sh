#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 16	
#SBATCH --mem 5G				# memory pool for all cores
#SBATCH --time=0-01:01:00				# time limit
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
output_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.bed
reference_gff=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3
reference_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3.bed
output_bedcoverage=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/bed_coverage.txt
# Source bedtools
source package bedops-2.4.28
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478 # bedtools 2.31.0 
source package c92263ec-95e5-43eb-a527-8f1496d56f1a # samtools Samtools 1.18 

#gff2bed < "$reference_gff" > "$reference_bed"
#bedtools bamtobed -i "$output_bam" > "$output_bed"
# sam to bam 
#samtools view -@ 16 -Sb -o "$output_bam" "$output_sam"
#bedtools coverage -a "$output_bed" -b "$reference_bed" > "$output_bedcoverage"


direct_rna_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb.bed
direct_rna_bed_renamed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_renamed.bed
intersected_read=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/intersected_read.bed
complete_read=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/complete_reads.bed
longest_complete_reads=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/longest_complete_reads.bed
sorted_complete_read=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/longest_complete_reads_sorted.bed
#cp "$direct_rna_bed" "$direct_rna_bed_renamed"
#sed -i 's/^NC_057794.1/Chr1A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057795.1/Chr1B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057796.1/Chr1D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057797.1/Chr2A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057798.1/Chr2B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057799.1/Chr2D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057800.1/Chr3A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057801.1/Chr3B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057802.1/Chr3D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057803.1/Chr4A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057804.1/Chr4B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057805.1/Chr4D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057806.1/Chr5A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057807.1/Chr5B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057808.1/Chr5D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057809.1/Chr6A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057810.1/Chr6B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057811.1/Chr6D/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057812.1/Chr7A/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057813.1/Chr7B/' "$direct_rna_bed_renamed"
#sed -i 's/^NC_057814.1/Chr7D/' "$direct_rna_bed_renamed"
# Remove NW_025226156.1 and NW_025226156.1 ?

bedtools intersect -a "$direct_rna_bed_renamed" -b "$reference_bed" -wa -wb -bed > "$intersected_read"
awk '$14 == "mRNA" { overlap = $3 - $2; annotation_length = $9 - $8; if (annotation_length > 0) { print $0, overlap / annotation_length } }' "$intersected_read" > "$complete_read"
num_cols=$(awk '{print NF; exit}' "$complete_read") # Get the last column 
sort -k${num_cols},${num_cols}nr "$complete_read" > "$sorted_complete_read"
