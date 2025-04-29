#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 1
#SBATCH --mem 5G				# memory pool for all cores
#SBATCH --time=0-001:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

# The ISOseq chromosome names are NC_ e.c.t 
# The Direct chromosomes are NC_ e.ct. 

# Input File Paths

iso_seq_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.bam
reference_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3.bed

direct_rna_sam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.minimap2.sam


# Resultant file paths 
iso_seq_bam_sorted=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/CAD_R.PacBio_sorted.bam
direct_rna_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb.minimap2.bam
direct_rna_bam_sorted=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb.minimap2.sorted.bam

# STEP 1 : Sort and Index Bam files 
source package c92263ec-95e5-43eb-a527-8f1496d56f1a # Samtools 1.18
##samtools sort -@ 16 -o "$iso_seq_bam_sorted" "$iso_seq_bam"
##samtools index -c -@ 16 "$iso_seq_bam_sorted"
##samtools sort -@ 16 -o "$direct_rna_bam_sorted" "$direct_rna_bam"
##samtools index -b -@ 16 "$direct_rna_bam_sorted"


# Direct RNA sam convert to bam , sort and index
# Need a header to the sam file 

reference_fa=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_assembly.fa
reference_fa_renamed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_assembly.fa
#echo "$reference_fa"
#echo "becomes"
#echo "$reference_fa_renamed"
#cp "$reference_fa" "$reference_fa_renamed"
#sed -i -e 's/^>Chr1A/>NC_057794.1/' \
# -e 's/^>Chr1B/>NC_057795.1/' \
# -e 's/^>Chr1D/>NC_057796.1/' \
# -e 's/^>Chr2A/>NC_057797.1/' \
# -e 's/^>Chr2B/>NC_057798.1/' \
# -e 's/^>Chr2D/>NC_057799.1/' \
# -e 's/^>Chr3A/>NC_057800.1/' \
# -e 's/^>Chr3B/>NC_057801.1/' \
# -e 's/^>Chr3D/>NC_057802.1/' \
# -e 's/^>Chr4A/>NC_057803.1/' \
# -e 's/^>Chr4B/>NC_057804.1/' \
# -e 's/^>Chr4D/>NC_057805.1/' \
# -e 's/^>Chr5A/>NC_057806.1/' \
# -e 's/^>Chr5B/>NC_057807.1/' \
# -e 's/^>Chr5D/>NC_057808.1/' \
# -e 's/^>Chr6A/>NC_057809.1/' \
# -e 's/^>Chr6B/>NC_057810.1/' \
# -e 's/^>Chr6D/>NC_057811.1/' \
# -e 's/^>Chr7A/>NC_057812.1/' \
# -e 's/^>Chr7B/>NC_057813.1/' \
# -e 's/^>Chr7D/>NC_057814.1/' \
# -e 's/^ChrMT/NC_036024.1/' "$reference_fa_renamed"
#echo "Finished rename"
samtools faidx "$reference_fa_renamed"
samtools view -b -T "$reference_fa_renamed" "$direct_rna_sam" > "$direct_rna_bam"



#samtools view -@ 16 -Sb "$direct_rna_sam" > "$direct_rna_bam"
#samtools sort -@ 16 -o "$direct_rna_bam_sorted" "$direct_rna_bam"
#samtools index -b -@ 16 "$direct_rna_bam_sorted"

## The conda environment includes the following packages
## rseqc v 5.0.4 # Note this is the package that contains the script geneBody_coverage.py
## ucsc-gtftogenepred ucsc-genepredtobed
#
#
## Convert the GTF file to BED12 using gtfToGenePred and genePredToBed
#source ~/.bashrc
#conda activate /hpc-home/mahony/miniconda3
#annotations_gff=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3
#reference_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3.bed
#annotations_genepred=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.genepred
#annotations_12_bed=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.12.bed
#annotations_12_bed_NC_names=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/iwgsc_refseqv2.1_annotation_200916_HC.12.NC.bed

#echo "gffread"
#gffread "$annotations_gff" -T -o "$reference_bed"
#head "$reference_bed"
#echo "gtfToGenePred"
#gtfToGenePred "$reference_bed" "$annotations_genepred"
#echo "genePredToBed"
#genePredToBed "$annotations_genepred" "$annotations_12_bed"
# conda deactivate

#echo "sed rename gff file chromosome names"
#sed -e 's/^1A/NC_057794.1/' \
#    -e 's/^1B/NC_057795.1/' \
#    -e 's/^1D/NC_057796.1/' \
#    -e 's/^2A/NC_057797.1/' \
#    -e 's/^2B/NC_057798.1/' \
#    -e 's/^2D/NC_057799.1/' \
#    -e 's/^3A/NC_057800.1/' \
#    -e 's/^3B/NC_057801.1/' \
#    -e 's/^3D/NC_057802.1/' \
#    -e 's/^4A/NC_057803.1/' \
#    -e 's/^4B/NC_057804.1/' \
#    -e 's/^4D/NC_057805.1/' \
#    -e 's/^5A/NC_057806.1/' \
#    -e 's/^5B/NC_057807.1/' \
#    -e 's/^5D/NC_057808.1/' \
#    -e 's/^6A/NC_057809.1/' \
#    -e 's/^6B/NC_057810.1/' \
#    -e 's/^6D/NC_057811.1/' \
#    -e 's/^7A/NC_057812.1/' \
#    -e 's/^7B/NC_057813.1/' \
#    -e 's/^7D/NC_057814.1/' \
#    -e 's/^MT/NC_036024.1/' "$annotations_12_bed" > "$annotations_12_bed_NC_names"
#echo "Sed Finished"
#
#bam_1=../intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_sorted.bam
#bam_2=../intermediate_data/CAD_R.PacBio_sorted.bam
# NC_057794.1
## It requires a bai index file not a csi index file but just altered the file ending
## cp  ../intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_sorted.bam.csi ../intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_sorted.bam.bai
#
#
#source ~/.bashrc
#conda activate /hpc-home/mahony/miniconda3
#geneBody_coverage.py -r "$annotations_12_bed" -i  "$bam_1" "$bam_2" -o out_RSEM_default_bams
#conda deactivate#