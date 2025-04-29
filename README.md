# Direct_RNA
Analysing ONT direct RNA sequencing data.

# Input data
/ei/data/reads/PDEV-88/2024_10_24_PSEQ2766
# Genome data 
/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/CAGE/input_data/chinese_spring_genome_data
# RNAseq data
/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/input_data/241114_CAD_CSref_TPM_root_rm0_from_BW.txt

# Anaylysis
1. Re-base calling and alignment perfomed using Dorado v0.7.2 . which is installed on the HPC already.
Using the simplex model - rna004_130bps_sup@v3.0.1 - as well as the modification model - rna004_130bps_sup@v3.0.1_m6A_DRACH@v1
With parameter to measure the polyA tail
run_dorado_model.sh, mapping results from samtools_flagstat.sh Splice information extracted with run_minimap2_splice.sh. Region of interest extracted with subsection_sam_file.sh

2. Mapping of ISO-Seq data to Chinese Spring with minimap2 for comparison with the direct RNA data.
run_mapping_iso_seq_data.sh

3. Plot the distribution of polyA tails and correlates with expression values.  
run_polyA.sh

4. Conversion of file formats and extracting a region of interest. run_file_conversions.sh 

5. Extract RNA modification information. run_m6A_modifications.sh

6. Bed tools interset to identify the distribution of complete transcripts. run_bedtools_intersect.sh

7. run_rseqc.sh

# Transcript visualisation 
Whole chromosomes or regions of interest  extracted in 4. can then be visualised in IGV.  

