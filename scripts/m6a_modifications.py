import re
import matplotlib.pyplot as plt
import seaborn as sns 
import numpy as np

# Define the path to your SAM file
sam_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/20241024_1559_MN32565_FBA19032_d1057abb_test.sam"
another_sam = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/output_dorado/modification/20241024_1559_MN32565_FBA19032_d1057abb.minimap2.sam"
plot_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/m6a.png"

# Initialize variables
normalized_positions = []  # To store normalized positions of m6A modifications
modifications_per_transcript = []  # To store the number of modifications per transcript

# Pattern for MM tag to extract m6A modification positions and CIGAR for transcript length

mm_pattern = re.compile(r"MM:Z:A\+[^\s,]+,([0-9,]+);")  # Matches MM tag for modifications on adenine

cigar_pattern = re.compile(r"\d+[MIDNSHP=X]")  # CIGAR pattern to extract transcript length

# Parse SAM file
with open(sam_file, "r") as file:
    for line in file:
        if line.startswith("@"):  # Skip header lines
            continue

        # Split line into fields
        fields = line.strip().split("\t")

        # Get CIGAR string and calculate transcript length
        cigar = fields[5]
        cigar_elements = re.findall(cigar_pattern, cigar)

        # Calculate transcript length by summing the lengths of 'M', 'I', and '=' (aligned bases)
        transcript_length = sum(int(element[:-1]) for element in cigar_elements if element[-1] in "MI=")

        # Sanity check for transcript length
        if transcript_length <= 0:
            print(f"Unrealistic transcript length: {transcript_length} in line:\n{line}")
            continue

        # Search for MM tag to get modification positions
        mm_match = mm_pattern.search(line)

        if mm_match:
            # Extract positions and calculate the normalized position
            positions = list(map(int, mm_match.group(1).split(",")))
            normalised_pos = [pos / transcript_length for pos in positions]
            
            # Check for normalized positions > 1
            if any(pos > 1 for pos in normalised_pos):
                print(f"Error: Normalized position > 1 detected. Skipping line.")
                print(f"Line: {line}")
                print(f"CIGAR: {cigar}, Elements: {cigar_elements}, Transcript Length: {transcript_length}")
                print(f"MM Match: {mm_match.group(1)}, Positions: {positions}")
                print(f"Normalized Positions: {normalised_pos}")
                continue

            # Store valid normalized positions
            normalized_positions.extend(normalised_pos)

            # Store the number of modifications for this transcript
            modifications_per_transcript.append(len(positions))

# Summary statistics
total_modifications = len(normalized_positions)
average_modifications = np.mean(modifications_per_transcript) if modifications_per_transcript else 0
median_modifications = np.median(modifications_per_transcript) if modifications_per_transcript else 0
average_position = np.mean(normalized_positions) if normalized_positions else 0
median_position = np.median(normalized_positions) if normalized_positions else 0
min_position = np.min(normalized_positions) if normalized_positions else 0
max_position = np.max(normalized_positions) if normalized_positions else 0

# Print results
print(f"Total modifications: {total_modifications}")
print(f"Average modifications per transcript: {average_modifications:.2f}")
print(f"Median modifications per transcript: {median_modifications:.2f}")
print(f"Average normalized position: {average_position:.2f}")
print(f"Median normalized position: {median_position:.2f}")
print(f"Min normalized position: {min_position:.2f}")
print(f"Max normalized position: {max_position:.2f}")

# Visualization
plt.figure(figsize=(10, 6))
sns.histplot(normalized_positions, bins=50, kde=True, color='blue')
plt.title("Distribution of Normalized m6A Modification Positions")
plt.xlabel("Normalized Position")
plt.ylabel("Frequency")
plt.savefig(plot_file)

exit()
# Step 3: Plot the Normalized Positions of m6A Modifications
plt.figure(figsize=(8, 6))
sns.violinplot(data=normalized_positions, color="skyblue")

# Labeling and title
plt.xlabel("Normalized Position of m6A Modification")
plt.ylabel("Density")
plt.xlim(0, 1)
plt.title("Distribution of Normalized m6A Modification Positions")

# Save the plot
plt.savefig(plot_file)

# Step 4: Calculate Statistics
total_transcripts = len(modifications_per_transcript)
total_modifications = sum(modifications_per_transcript)
average_modifications_per_transcript = total_modifications / total_transcripts if total_transcripts > 0 else 0

# Output the statistics
print(f"Total Transcripts: {total_transcripts}")
print(f"Total m6A Modifications: {total_modifications}")
print(f"Average m6A Modifications per Transcript: {average_modifications_per_transcript:.2f}")
print(f'file saved at {plot_file}')
