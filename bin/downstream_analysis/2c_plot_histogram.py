import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the CSV file
file_path = "Gobi2025microRNA/shared_targets.csv"
df = pd.read_csv(file_path, header=None)

# Extract the Overlap Size (7th column, index 6 since Python is 0-based)
overlap_sizes = df.iloc[:, 6]

# Define 10 bins (auto-ranges from min to max overlap size)
bins = np.linspace(overlap_sizes.min(), overlap_sizes.max(), 11)

# Plot histogram
plt.hist(overlap_sizes, bins=bins, edgecolor='black', alpha=0.7)
plt.yscale('log')


# Labels and title
plt.xlabel("Overlap Size")
plt.ylabel("Count")
plt.title("Histogram of Overlap Sizes log-distribution (Divided into 10 Bins)")

# Show the plot
plt.show()