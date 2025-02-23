import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Get files
files = {
    "MOORE": "moore.isomir.peaks.family.filtered.bed",  # Ersetze mit dem Dateinamen des ersten Datensets
    "BJERKE": "bjerke.isomir.peaks.family.filtered.bed"  # Ersetze mit dem Dateinamen des zweiten Datensets
}

# Read in Files with different format
def load_bed_file(file_path):
    df = pd.read_csv(file_path, sep="\t", header=None)

    # Check numer of columns
    if df.shape[1] >= 12:
        df.columns = ["chr", "start", "end", "miRNA", "score", "strand",
                      "ref_chr", "ref_start", "ref_end", "gen", "region", "extra"]
    else:
        df.columns = ["chr", "start", "end", "miRNA", "score", "strand"]

    # Make numeric values out of column score
    df["score"] = pd.to_numeric(df["score"], errors="coerce")

    return df

# Classification of miRNA types
def classify_miRNA(miRNA):
    if "3p" in miRNA and "5p" in miRNA:
        return "3p5p"
    elif "3p" in miRNA:
        return "3p"
    elif "5p" in miRNA:
        return "5p"
    else:
        return "canonical"

# List to store data for combined plot2
combined_plot_data = []

for dataset_name, file_path in files.items():

    # Read File
    df = load_bed_file(file_path)

    # plot 1: Top 10 miRNAs with most Interactions/expressions
    miRNA_counts = df.groupby("miRNA")["score"].sum().reset_index()
    miRNA_counts = miRNA_counts.sort_values(by="score", ascending=False).head(10)

    # Create plot for both datasets
    plt.figure(figsize=(10, 6))
    plt.barh(miRNA_counts["miRNA"], miRNA_counts["score"], color=plt.get_cmap("tab10").colors[:10])
    plt.xlabel("Count")
    plt.ylabel("miRNA ID")
    plt.title(f"Top 10 miRNA with most Interactions ({dataset_name})")
    plt.gca().invert_yaxis()

    # Delete Tick Marks
    plt.tick_params(axis="x", which="both", bottom=False, top=False)
    plt.tick_params(axis="y", which="both", left=False, right=False)

    # Save the first plots
    plt.savefig(f"miRNA_expression_plot_{dataset_name}.png", dpi=300, bbox_inches="tight")
    plt.show()

    # plot 2: miRNA-types for both data sets
    df["dataset"] = dataset_name
    df["isomir_type"] = df["miRNA"].apply(classify_miRNA)

    # Aggregation of miRNA types
    plot_data = df.groupby(["dataset", "isomir_type"])["score"].sum().reset_index()

    # Save to be able to combine the plots
    combined_plot_data.append(plot_data)

# plot 2: combined plot for both data sets
final_plot_data = pd.concat(combined_plot_data)

plt.figure(figsize=(8, 6))
sns.barplot(data=final_plot_data, x="dataset", y="score", hue="isomir_type", palette="tab10")

plt.xlabel("Dataset")
plt.ylabel("miRNA Expression Score")
plt.title("miRNA type per dataset")


# save combined plot
plt.savefig("combined_miRNA_type_per_dataset.png", dpi=300, bbox_inches="tight")
plt.show()