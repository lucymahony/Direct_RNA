import pandas as pd


def count_polyA_lengths(sam_file, output_file, quality_threshold):
    """
    Input:
        sam_file: File of the mapped ONT data from the modified model base caller
    Output:
        output_file: File path to save the counts of different lengths of polyAs, but rather than count it is 
        a ratio of number of counts out of all transcripts. 
    """
    # Initialize a dictionary for counting polyA tail lengths and a variable for total count
    polyA_counts = {}
    total_above_threshold = 0

    # Read through the SAM file
    with open(sam_file, "r") as f:
        for line in f:
            if line.startswith("@"):
                continue  # Skip header lines

            fields = line.strip().split("\t")

            # Extract FLAG and check if it's a unique hit (not secondary or supplementary)
            flag = int(fields[1])
            is_unique_hit = flag & 0x900 == 0  # Filters out secondary (2048) and supplementary (256) alignments

            if is_unique_hit:
                polyA_length = None
                quality_score = None

                # Manually search for pt:i (polyA tail length) and qs:f (quality score) tags
                for field in fields:
                    if field.startswith("pt:i:"):
                        polyA_length = int(field.split(":")[2])  # Extract polyA length
                    elif field.startswith("qs:f:"):
                        quality_score = float(field.split(":")[2])  # Extract quality score

                # If both tags are found and quality meets the threshold, count the polyA length
                if polyA_length is not None and quality_score is not None and quality_score >= quality_threshold:
                    total_above_threshold += 1  # Increment the total count above threshold
                    if polyA_length in polyA_counts:
                        polyA_counts[polyA_length] += 1
                    else:
                        polyA_counts[polyA_length] = 1

    # Calculate the ratio of reads for each polyA length
    polyA_ratios = {length: count / total_above_threshold for length, count in polyA_counts.items()}

    # Save results to a file
    with open(output_file, "w") as out:
        for length, ratio in sorted(polyA_ratios.items(), key=lambda x: x[1], reverse=True):
            out.write(f"{ratio:.4f}\t{length}\n")

    # Display the result
    print(f"Results saved to {output_file}")


def read_in_rna_seq_data(rna_seq_file_path, tpm_threshold):
    """
    Read the TPM values from the rna seq file.
    Return a df of (gene name, tpm) for genes  with an average TPM > 0.5 
    """
    # Skip the first row, column names are Gene, CADR1, CADR2, CADR3
    tpm = pd.read_csv(rna_seq_file_path, sep='\t', skiprows=1, names = ['Gene', 'CADR1', 'CADR2', 'CADR3'])
    print(tpm.head())
    # Add column of averages 
    tpm['Mean_TPM'] = ((tpm['CADR1'] + tpm['CADR2'] + tpm['CADR3']) / 3).astype(float)
    original_number_of_genes = tpm.shape[0]
    # Filter with tpm threshold
    tpm = tpm[tpm['Mean_TPM'] >= 0.5]
    filtered_number_of_genes = tpm.shape[0]
    print(f'Filtering removes {original_number_of_genes - filtered_number_of_genes} genes from the dataset of {filtered_number_of_genes} genes')
    list_genes = tpm['Gene'].to_list()
    return tpm, list_genes

def find_gene_in_sam_file(sam_file, gene):
    with open(sam_file, "r") as f:
            for line in f:
                if line.startswith("@"):
                    continue  # Skip header lines

                fields = line.strip().split("\t")

                # Extract FLAG and check if it's a unique hit (not secondary or supplementary)
                flag = int(fields[1])
                is_unique_hit = flag & 0x900 == 0  # Filters out secondary (2048) and supplementary (256) alignments

                if is_unique_hit:
                    print(fields)
                    exit()





if __name__ == "__main__":
    # Define the input SAM file and quality threshold
    sam_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_simplex.sam"
    output_file= "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/polyA_ratios.txt"
    quality_threshold = 0.0  # Set the quality threshold
    # File path to RNAseq values and set threshold of average TPM
    rna_seq= "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/input_data/241114_CAD_CSref_TPM_root_rm0_from_BW.txt"
    tpm_threshold = 0.5

    #tpm, list_genes = read_in_rna_seq_data(rna_seq, tpm_threshold)
    #print(list_genes[0]) # TraesCS1A02G000200 
    find_gene_in_sam_file(sam_file, 'TraesCS1A02G000200')