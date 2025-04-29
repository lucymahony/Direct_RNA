import matplotlib.pyplot as plt

input_file= "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/polyA_ratios.txt"
plot_file_xlims="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/plotA_xlims.png"
plot_file="/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/Direct_RNA/intermediate_data/plotA.png"

data = []

with open(input_file, "r") as file:
    for line in file:
        count, length = line.strip().split()
        data.append((int(length), float(count)))  # Convert to int and float for sorting and plotting

# Sort the data by length
data.sort()  # Sorts by the first element (length) in each tuple by default
lengths, counts = zip(*data)  # Unpack sorted data into separate lists

# Plot
#plt.plot(lengths, counts, marker='o', color='#0f62f3')  # Adding markers for clarity
#plt.xlabel("PolyA Tail Length (bp)")
#plt.ylabel("Ratio of reads")
#plt.xlim(0, 100)
#plt.title("Distribution of Poly(A) Tail Lengths")
#plt.savefig(plot_file_xlims, dpi=900)

# Plot No X lims
plt.plot(lengths, counts, marker='o', color='#0f62f3')  # Adding markers for clarity
plt.xlabel("PolyA Tail Length (bp)")
plt.ylabel("Ratio of reads")
plt.title("Distribution of Poly(A) Tail Lengths")
plt.savefig(plot_file, dpi=900)


print(f"File saved to {plot_file}")

