#!/bin/bash
#SBATCH -p ei-medium
#SBATCH -c 16	
#SBATCH --mem 50G				# memory pool for all cores
#SBATCH --time=0-01:01:00				# time limit
#SBATCH --output output_%x		# STDOUT and STDERR
#SBATCH --mail-type=END,FAIL			# notifications for job done & fail
#SBATCH --mail-user=lucy.mahony@earlham.ac.uk	# send-to address

# Convert BAM files to useful formats for downstream analysis:
#- sorted bam files
#- indexed bam files
#- sam files 
#- bed files 
#- sorted and indexed bam files for specific genomic region. 

# File Paths
# Input
reference=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna
reference_gff=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data/iwgsc_refseqv2.1_annotation_200916_HC.gff3
direct_rna_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.bam
iso_seq_bam=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/mapped_iso_seq_transcripts/CAD_R.PacBio.bam
cage_tag_clusters=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/intermediate_data/snakemake_intermediate_data/R9_TC.bed
# Output 
intermediate_data_directory=/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data
reference_index=GCF_018294505.1_IWGSC_CS_RefSeq_v2.1_genomic.fna.mmi
bedcoverage=bed_coverage.txt

#Set parameters
chr=1A
other_chr=NC_057794.1d # Corresponding chromosome in the reference format
start=299443927
stop=29945082

# Source packages
source package bedops-2.4.28
source package b0ed0698-358b-4c9b-9d21-603ea8d6e478 # bedtools 2.31.0 
source package c92263ec-95e5-43eb-a527-8f1496d56f1a # samtools Samtools 1.18 

# Set up output directory for the specific region
region_file="${chr}_${start}_${stop}"
mkdir -p "${intermediate_data_directory}/${region_file}"
# Extract region sequence of assembly
sed -n "/^>${other_chr}\b/,/^>/p" "$reference" | sed '$d' > "${intermediate_data_directory}/${region_file}/region.fna"
# Filter GFF for region of interest
awk -v chr=Chr"$chr" -v start="$start" -v stop="$stop" \
    '$1 == chr && $4 >= start && $5 <= stop' "$reference_gff" \
    > "${intermediate_data_directory}/${region_file}/annotation.gff"

# Filter CAGE tag clusters for region of interest
awk -v chr=Chr"$chr" -v start="$start" -v stop="$stop" \
    '$1 == chr && $2 >= start && $3 <= stop' "$cage_tag_clusters" \
    > "${intermediate_data_directory}/${region_file}/cage_tag_clusters_region.bed"

# Update chromosome name to match reference format in the filtered CAGE file
sed -i -e "s/\b${chr}\b/${other_chr}/g" "${intermediate_data_directory}/${region_file}/cage_tag_clusters_region.bed"

# Filter ISO Seq BED file for region of interest and Sort and index the bam files 
bam_files=("direct_rna_bam" "iso_seq_bam")
sorted_output_names=("direct_rna_sorted.bam" "iso_seq_sorted.bam")
bed_output_names=("direct_rna.bed" "iso_seq.bed")
region_bed_output_names=("direct_rna_region.bed" "iso_seq_region.bed")

# Loop through BAM files for sorting, indexing, converting to BED, and filtering by region
for i in "${!bam_files[@]}"; do
    bam_file=${bam_files[$i]}
    sorted_bam="${intermediate_data_directory}/${region_file}/${sorted_output_names[$i]}"
    bed_file="${intermediate_data_directory}/${region_file}/${bed_output_names[$i]}"
    region_bed="${intermediate_data_directory}/${region_file}/${region_bed_output_names[$i]}"
    
    # Sort BAM
    samtools sort -o "$sorted_bam" "${!bam_file}"
    
    # Index BAM
    samtools index -c "$sorted_bam"
    
    # Convert BAM to BED
    bedtools bamtobed -i "${!bam_file}" > "$bed_file"
    
    # Filter BED for region of interest
    awk -v chr="$other_chr" -v start="$start" -v stop="$stop" \
        '$1 == chr && $2 >= start && $3 <= stop' "$bed_file" > "$region_bed"
    
    echo "Processing completed for ${!bam_file}: sorted, indexed, converted to BED, and filtered."
done
