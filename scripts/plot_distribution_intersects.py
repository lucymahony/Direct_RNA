import matplotlib.pyplot as plt

# File path
file_path = "../intermediate_data/complete_reads.bed"
fig_file_path="../intermediate_data/complete_reads_distribution.pdf"
# Initialize a list to store the overlap-to-annotation length ratios
ratios = []

# Read the file and extract the last column (overlap/annotation length)
with open(file_path, "r") as file:
    for line in file:
        fields = line.strip().split()
        # The last column contains the overlap/annotation length value
        ratio = float(fields[-1])
        ratios.append(ratio)

# Plot the distribution as a histogram
plt.figure(figsize=(10, 6))
plt.hist(ratios, bins=1000, edgecolor='black', alpha=0.7)
plt.xlabel("Overlap / Annotation Length Ratio")
plt.xlim(0,1)
plt.ylabel("Frequency")
plt.title("Distribution of Overlap / Annotation Length Ratios")
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Show the plot
plt.savefig(fig_file_path, dpi=900)

