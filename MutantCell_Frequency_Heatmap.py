##################################################################################################
##        Mutant Cell Frequency Visualization (Heatmap)                                         ##
##                                                                                              ##
##  This script generates heatmap of mutant cell frequencies of all targets for a sample        ##
##                                                                                              ##
##################################################################################################


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.patheffects import withStroke
from matplotlib import rcParams

rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.labelsize"] = "7"
rcParams["axes.titlesize"] = "7"
rcParams["font.size"] = "7"
rcParams["xtick.labelsize"] = "7"
rcParams["ytick.labelsize"] = "7"
rcParams["axes.linewidth"] = "0.5"
rcParams["xtick.major.width"] = "0.5"
rcParams["ytick.major.width"] = "0.5"
rcParams["xtick.major.size"] = "3"
rcParams["ytick.major.size"] = "3"
rcParams["mathtext.fontset"] = "custom"
rcParams["mathtext.rm"] = "serif"
rcParams["pdf.fonttype"] = 42
rcParams["xtick.direction"] = "out"
rcParams["ytick.direction"] = "out"

datadir = "Data"
outdir = "Results"
os.makedirs(outdir, exist_ok=True)

mut_status = "best_pred_adj"
sample = "CLL6"  # 'CLL3', 'CLL4', 'CLL5', 'CLL6', 'CLL7'
non_mut_group = "is_non_mut"
absolute_cell_count_threshold = 0
reassign_unprofiled_to_closest_genotype = True

gen_df = pd.read_csv(
    os.path.join(
        datadir,
        "CLL_RT_GenotypingProfiles.csv",
    ),
)

# === Preprocess === #
# ! Filter: No cell groups where mutants are not expected
gen_df = gen_df[gen_df[non_mut_group] == False]
# ! Filter: Cells in the expected sample
gen_df = gen_df[gen_df["sample"] == gen_df["expected_sample"]]
samples = gen_df["expected_sample"].unique()
print("Sample:", sample)

# ===================================================================================================
# Mutant cell frequency calculation & Visualization
# ===================================================================================================
targets = gen_df[gen_df["expected_sample"] == sample]["target"].unique()
if len(targets) <= 1:
    raise ValueError(f"{sample} has only one target: {targets}")
print("Targets:", targets)
genotype_data = gen_df.query("expected_sample == @sample & target in @targets")

# * === 1-1) Replace unprofiled according to unique closest genotype === * #
# Assign genotype if there is a unique genotype that is the closest
genotype_data[mut_status] = genotype_data[mut_status].str.replace(
    "FalsePositive", "Unprofiled"
)
"""
If there are cells that are Unprofiled for any of the targets,
assign the genotype to the closest (most similar) genotype, 
only if there is a single genotype that is most similar. (Unique closest genotype)
"""
# Assign Closest Genotype if Unprofiled Exists
# Convert to cell-by-gene matrix where each entry indicates the genotype value (MUT, WT, Unprofiled)
cell_gene_matrix = pd.pivot_table(
    genotype_data, index="cell_id", columns="target", values=mut_status, aggfunc="first"
).fillna("Unprofiled")
assert cell_gene_matrix.isna().sum().sum() == 0

# Define MUT=1, WT=0, Unprofiled=-1 for similarity calculation
cell_gene_matrix = cell_gene_matrix.applymap(
    lambda x: {"MUT": 1, "WT": 0, "Unprofiled": -1}[x]
)
print(
    "# of cells with any unprofiled target:", (cell_gene_matrix == -1).any(axis=1).sum()
)

# Assign closest genotype if Unprofiled exists
if reassign_unprofiled_to_closest_genotype:
    for idx, row in cell_gene_matrix.iterrows():
        if -1 in row.values:
            # Calculate similarity with other genotypes
            # rows that do not have any Unprofiled (fetch rows that don't have -1)
            candidate_genotypes = cell_gene_matrix.drop(index=idx)[
                ~(cell_gene_matrix == -1).any(axis=1)
            ]
            candidate_genotypes = candidate_genotypes.drop_duplicates()

            similarity_scores = candidate_genotypes.apply(
                lambda x: np.sum((x == row) & (x != -1)), axis=1
            )
            max_similarity = similarity_scores.max()
            # Assign genotype only if there is a single most similar genotype
            if (similarity_scores == max_similarity).sum() == 1:
                closest_genotype_idx = similarity_scores.idxmax()
                closest_genotype = cell_gene_matrix.loc[closest_genotype_idx]
                for target in cell_gene_matrix.columns:
                    if row[target] == -1:
                        mask = (genotype_data["cell_id"] == idx) & (
                            genotype_data["target"] == target
                        )
                        if mask.sum() == 0:
                            continue
                        genotype_data.loc[
                            mask,
                            mut_status,
                        ] = genotype_data.loc[
                            (genotype_data["cell_id"] == closest_genotype_idx)
                            & (genotype_data["target"] == target),
                            mut_status,
                        ].values[0]
                        cell_gene_matrix.loc[idx, target] = closest_genotype[target]
    print(
        "# of cells with any unprofiled target:",
        (cell_gene_matrix == -1).any(axis=1).sum(),
    )

# * === 1-2) Remaining 'Unprofiled' becomes 'WT' === * #
# === Get Genotypes that exist within the single-cell genotype data === #
# Rename columns that contain '_' and remove them
cell_gene_matrix.columns = cell_gene_matrix.columns.str.replace("_", "")

# Replace -1 (Unprofiled) with 0 (WT)
cell_gene_matrix = cell_gene_matrix.replace(-1, 0)

# Rename Index to Genotype instead of cell barcode
# create 'genotype' column for each cell: concatenate target names where the value is 1 with '_'
cell_gene_matrix["genotype"] = cell_gene_matrix.apply(
    lambda x: "_".join(x.index[x == 1]), axis=1
)

# if genotype empty(all targets WT), assign 'WT'
cell_gene_matrix["genotype"] = cell_gene_matrix["genotype"].replace("", "WT")
genotype_matrix = cell_gene_matrix.reset_index(drop=True).set_index("genotype")

# count each genotype
genotype_counts = genotype_matrix.drop("WT").index.value_counts()
genotype_ratio = genotype_matrix.drop("WT").index.value_counts(normalize=True)

# Drop duplicate genotypes
genotype_matrix = genotype_matrix[~genotype_matrix.index.duplicated(keep="first")]
genotype_matrix = genotype_matrix.sort_values(by=genotype_matrix.columns.tolist())
cell_gene_matrix = cell_gene_matrix[cell_gene_matrix["genotype"] != "WT"]

# * === Draw cell / genotype counts heatmap === * #
gene_by_genotype = pd.pivot_table(
    cell_gene_matrix,
    index="genotype",
    aggfunc="sum",
)  # .drop("WT")

# Create a new gene-by-gene matrix filled with zeros
all_genes = list(gene_by_genotype.columns)
gene_by_gene = pd.DataFrame(0, index=all_genes, columns=all_genes)

# First, process diagonal elements (cells with only one gene mutated)
for gene in all_genes:
    # Look for the genotype that consists of only this gene
    if gene in gene_by_genotype.index:  # If the single gene exists as a genotype
        gene_by_gene.at[gene, gene] = gene_by_genotype.loc[gene, gene]

# Then process co-occurrences from multi-gene genotypes
for genotype, row in gene_by_genotype.iterrows():
    # Parse the genotype to get individual genes
    genotype_genes = [g for g in genotype.split("_") if g in all_genes]

    if len(genotype_genes) > 1:  # Only process multi-gene genotypes
        count = row[
            genotype_genes[0]
        ]  # Count should be the same for all genes in the genotype

        # For each pair of genes in this genotype, add the count to their co-occurrence entries
        for i, gene_i in enumerate(genotype_genes):
            for gene_j in genotype_genes[i + 1 :]:
                gene_by_gene.at[gene_i, gene_j] += count
                gene_by_gene.at[gene_j, gene_i] += count  # Matrix is symmetric


# Convert to percentage
total_n_cells = genotype_counts.sum()
print(total_n_cells)
gene_by_gene_percent = gene_by_gene / total_n_cells * 100

# Sum of upper triangle (excluding diagonal)
upper_triangle_with_diag_sum = np.triu(gene_by_gene_percent.values).sum()
print(f"Sum including diagonal: {upper_triangle_with_diag_sum:.2f}%")

gene_by_genotype = gene_by_gene_percent.copy()

# Calculate the sum of each gene across all genotypes
gene_sums = gene_by_gene.sum(axis=0)

# Calculate the sum of each genotype across all genes
genotype_sums = gene_by_gene.sum(axis=1)

# Add the sums as margins to the dataframe
gene_by_genotype_with_margins = gene_by_gene.copy()
gene_by_genotype_with_margins.loc["Total"] = gene_sums
gene_by_genotype_with_margins["Total"] = genotype_sums
gene_by_genotype_with_margins.loc["Total", "Total"] = gene_sums.sum()


fontsize = 12

# Calculate the sum of each gene across all genotypes
gene_sums = gene_by_genotype.sum(axis=0)

# Calculate the sum of each genotype across all genes
genotype_sums = gene_by_genotype.max(axis=1)

# Create the figure and gridspec layout
fig = plt.figure(figsize=(8, 6))
gs = gridspec.GridSpec(
    2, 2, width_ratios=[4, 1], height_ratios=[1, 4], wspace=0.05, hspace=0.05
)

# Create the heatmap
ax = plt.subplot(gs[1, 0])
ax = sns.heatmap(
    gene_by_genotype,
    cmap="Blues",  # "viridis",
    ax=ax,
    cbar=True,
    vmax=15,
    vmin=0,
    annot=False,
    fmt=".2f",  # fmt="g",
    annot_kws={
        "fontsize": fontsize,
        "weight": "bold",
        "color": "0.2",
        "path_effects": [withStroke(linewidth=2.0, foreground="white")],
        "clip_on": False,
    },
    linewidths=1.0,
    linecolor="black",
)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(0.5)
    spine.set_color("black")

plt.xticks(
    size=fontsize,
    rotation=90,
)
plt.yticks(
    size=fontsize,
    rotation=0,
)
# Add labels and title
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title(sample)

fpath = os.path.join(
    outdir, f"{sample}_MUT_frequency_heatmap.pdf"
)  # #muts out of all mutants
plt.savefig(fpath, dpi=300, bbox_inches="tight")
plt.show()
