##################################################################################################
##        Genotyping Across Cell States Visualization (Distributions)                           ##
##################################################################################################
##  This script generates heatmaps and stacked barplots showing per-target MUT vs WT            ##
##  distributions across defined biological groupings (cell states / cell types)                ##
##                                                                                              ##
## Input Data: CLL_RT_GenotypingProfiles.csv (Genotyping profiles refined by ML)                ##
##                                                                                              ##
## Generated Outputs (saved under `Results/GEN_Visualizations/`)                                ##
##  - Heatmap_MUT_ratio_per_CellState.pdf           : Heatmap of normalized mutation ratios by  ##
##                                                    cell state (first analysis block)         ##
##  - Barplot_MUT_WT_ratio_<TARGET>.pdf (per target): Stacked barplots MUT vs WT counts/ratios  ##
##  - Heatmap_MUT_ratio_per_CellType.pdf           : Heatmap across cell types (CLL / LBCL)    ##
##                                                                                              ##
##################################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from matplotlib.patheffects import withStroke

datadir = "Data/"
outdir = "Results"
os.makedirs(outdir, exist_ok=True)

# ===================================================================================================
# Mutant cell proportion across cell states and targets
# ===================================================================================================
cell_eda = pd.read_csv(
    os.path.join(
        datadir,
        "CLL_RT_GenotypingProfiles.csv",
    ),
)

# ! Remove unwanted targets
targets_to_remove = [
    "TNFAIP3_Mut1",
    "DDX3X_Mut1",
    "DDX3X_Mut2",
    "HIST1H1C",
]
cell_eda = cell_eda[~cell_eda["target"].isin(targets_to_remove)]

# ! Remove unwanted cell types
cell_group_col = "cell_state"
cell_group_to_exclude = ["CD4 T", "CD8 T", "Myeloid"]
cell_eda_filt = cell_eda.query(f"{cell_group_col} not in {cell_group_to_exclude}")

# ! Keep only the genotyped cells
cell_eda_filt = cell_eda_filt[cell_eda_filt["best_pred_adj"].isin(["WT", "MUT"])]

# Count #genotyped cells in each cell state
cell_state_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col])
    .size()
    .reset_index()
    .rename(columns={0: "n_cellstate_cells"})
)

# Count #MUT cells in each cell state
cell_state_mut_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col, "best_pred_adj"])
    .size()
    .reset_index()
    .query('best_pred_adj == "MUT"')
    .rename(columns={0: "n_cellstate_mut_cells"})
)
# Merge Total Cell Count and MUT Cell Count
cell_state_count = cell_state_count.merge(
    cell_state_mut_count.drop(columns="best_pred_adj"),
    on=[
        "target",
        "sample",
        cell_group_col,
    ],
    how="left",
).fillna(0)

# Count WT cells in each cell state
cell_state_wt_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col, "best_pred_adj"])
    .size()
    .reset_index()
    .query('best_pred_adj == "WT"')
    .rename(columns={0: "n_cellstate_wt_cells"})
)
# Merge Total Cell Count and WT Cell Count
cell_state_count = cell_state_count.merge(
    cell_state_wt_count.drop(columns="best_pred_adj"),
    on=[
        "target",
        "sample",
        cell_group_col,
    ],
    how="left",
).fillna(0)

# Calculate the mutation ratio
cell_state_count["mut_ratio"] = (
    cell_state_count["n_cellstate_mut_cells"] / cell_state_count["n_cellstate_cells"]
)
# Calculate the WT ratio
cell_state_count["wt_ratio"] = (
    cell_state_count["n_cellstate_wt_cells"] / cell_state_count["n_cellstate_cells"]
)

cell_state_count.sort_values(["target", "sample", "cell_state"])

# Make targets unique
tmp = cell_state_count.groupby(["target", "sample"]).size().reset_index()
assert tmp["target"].nunique() == tmp.shape[0]
assert (cell_state_count.value_counts(["target", "sample", "cell_state"]) == 1).all()

# Normalize MUT ratio to the total MUT ratio for that target & cell state
mut_ratio_sum = cell_state_count.groupby(["target"])["mut_ratio"].sum().reset_index()
# map mut ratio sum to target-cell state pairs
cell_state_count["mut_ratio_norm"] = cell_state_count.apply(
    lambda x: x["mut_ratio"]
    / mut_ratio_sum[(mut_ratio_sum["target"] == x["target"])]["mut_ratio"].values[0],
    axis=1,
)
cell_state_count.fillna(0, inplace=True)

# Pivot the table to get the cell state count matrix
cell_state_count_matrix = (
    cell_state_count.pivot_table(
        index=["target"],
        columns=cell_group_col,
        values="mut_ratio_norm",
    )
    .reset_index()
    .fillna(0)
)
# Pivot the table to get raw mut_ratio values (for annotations)
cell_state_count_matrix_raw = (
    cell_state_count.pivot_table(
        index=["target"],
        columns=cell_group_col,
        values="mut_ratio",
    )
    .reset_index()
    .fillna(0)
)

data = cell_state_count_matrix.set_index(["target"]).T
data_raw = cell_state_count_matrix_raw.set_index(["target"]).T  # For annotations


target2sample = cell_eda[["target", "expected_sample"]].drop_duplicates()
assert target2sample["target"].nunique() == target2sample.shape[0]
target2sample = target2sample.set_index("target")["expected_sample"].to_dict()

# Create the ordered list of columns based on samples
unique_samples = sorted(
    list(dict.fromkeys([target2sample[col] for col in data.columns]))
)

# # * === Order Input Matrix === * #
# new_order = []
# for sample in unique_samples:
#     sample_targets = [
#         target
#         for target, samp in target2sample.items()
#         if samp == sample and target in data.columns
#     ]
#     new_order.extend(sorted(sample_targets))

# # Reorder the data before creating the clustermap
# data_reordered = data[new_order]
# data_raw_reordered = data_raw[new_order]  # For annotations


# * === Order Input Matrix (Manual) === * #
col_order = [
    "EMD",
    "B2M",
    "POU2F2_Mut1",
    "MCL1",
    "PLCG2_Mut1",
    "PLCG2_Mut2",
    "PLCG2_Mut3",
    "JAK1",
    "PRKDC",
    "SRRM2",
    "JUNB",
    "IKZF3",
    "SF3B1",
    "IRF8_Mut1",
    # "HIST1H1C",
    "BTK",
]
row_order = [
    "DZ.LZ GCB",
    "LZ GCB",
    "Proliferating",
    "Stress Response",
    "Inflammatory",
    "CCND2$^{hi}$",
    "Immune Response",
    "Plasma B",
]

# reorder the data
data_reordered = data.reindex(row_order)
data_reordered = data_reordered[col_order]
data_raw_reordered = data_raw.reindex(row_order)  # For annotations
data_raw_reordered = data_raw_reordered[col_order]  # For annotations

# Create column colors based on the reordered data
sample_order = [target2sample[col] for col in data_reordered.columns]
colors = sns.color_palette("hls", len(unique_samples))
sample_colors = dict(zip(unique_samples, colors))
col_colors = [sample_colors[sample] for sample in sample_order]

fontsize = 16


# Create the clustermap with the reordered data
ax = sns.clustermap(
    data_reordered,  # Use pre-reordered data
    cmap="Blues",
    cbar_pos=(1.12, 0.45, 0.03, 0.23),
    cbar_kws={"label": "%Mutation cells"},
    figsize=(12, 7),
    row_cluster=False,
    col_cluster=False,
    dendrogram_ratio=(0.0, 0.00001),
    col_colors=col_colors,
    col_linkage=None,
    # annot=data_raw_reordered.values, #! If you want to show raw mutation ratios
    # fmt=".2f",
    # annot_kws={
    #     "fontsize": fontsize,
    #     "weight": "bold",
    #     "color": "grey",
    #     "path_effects": [withStroke(linewidth=1, foreground="white")],
    # },
    vmax=0.4,
)

# Set the Colorbar label font size
cbar = ax.ax_cbar
cbar.set_ylabel(
    cbar.get_ylabel(), fontsize=fontsize + 2, fontweight="bold", labelpad=20
)  # Set font size and weight
cbar.set_yticklabels(
    cbar.get_yticklabels(), fontsize=fontsize, weight="bold"
)  # Set the colorbar tick labels font sizes

# Legend
legend_elements = [
    Patch(facecolor=sample_colors[sample], label=sample) for sample in unique_samples
]
leg = plt.legend(
    handles=legend_elements,
    title="Samples",
    title_fontsize=fontsize + 2,
    bbox_to_anchor=(3.46, 2.4),
    loc=1,
    ncol=1,
    frameon=False,
    labelcolor="0.2",
    prop={"weight": "bold", "size": fontsize},
)
leg.get_title().set_weight("bold")

# Remove x and y labels
ax.ax_heatmap.set_xlabel("Target", fontsize=fontsize + 2, weight="bold", labelpad=10)
ax.ax_heatmap.set_ylabel(
    "Cell State", fontsize=fontsize + 2, weight="bold", labelpad=10
)

# Set the axes tick labels font sizes
ax.ax_heatmap.set_yticklabels(
    ax.ax_heatmap.get_yticklabels(), fontsize=fontsize, weight="bold"
)
ax.ax_heatmap.set_xticklabels(
    ax.ax_heatmap.get_xticklabels(), fontsize=fontsize, weight="bold"
)
# y축 눈금 레이블 회전
for label in ax.ax_heatmap.get_yticklabels():
    label.set_rotation(0)
    label.set_ha("left")
# x축 눈금 레이블 회전
for label in ax.ax_heatmap.get_xticklabels():
    label.set_rotation(90)
    # label.set_ha("right")
# Title
ax.ax_heatmap.set_title(
    "Mutation Ratio in Cell States",
    fontsize=fontsize + 2,
    weight="bold",
    pad=20,
    color="0.2",
)
fig = ax.figure

fig.savefig(
    os.path.join(outdir, "Heatmap_MUT_ratio_per_CellState.pdf"),
    bbox_inches="tight",
    dpi=300,
)

plt.show()


# ===================================================================================================
# MUT/WT proportion across cell states for each target (Bar Plot)
# ===================================================================================================
# ====== Plotting Ratio in Stacked Bar Plot with Raw Values Annotations (Ratios, Counts Precalculated) ====== #
def plot_stacked_barplot_precalculated(
    df,
    x="sample",
    y_cols=["mut_ratio", "wt_ratio"],
    annot_cols=["mut_counts", "wt_counts"],
    figsize=(6, 5),
    fontsize=10,
    linewidth=1.5,
    palette=None,
    legend_labels=None,
    legend_title="",
    title="",
    annotate=True,
    ax=None,
):
    """
    Create a stacked bar plot showing the distribution of data and the raw counts across a categorical variable.
    Designed for df with ratio and counts already calculated.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame containing the data to visualize with columns for cell states,
        cell type ratios, and cell counts.
    x : str, default="cell_state"
        Column name in df for x-axis categories (cell states).
    y_cols : list, default=["mut_ratio", "wt_ratio"]
        Column names in df for the stacked bar values (proportions).
    annot_cols : list, default=["n_cellstate_mut_cells", "n_cellstate_wt_cells"]
        Column names in df for annotation values (cell counts).
    legend_labels : list, default=None
        Labels to display in the legend. If None, uses the labels from the plot.
    fontsize : int, default=10
        Base font size for text elements.
    palette : str, list, or None, default=None
        Colors for the bars. Can be:
        - None: Uses default seaborn color palette
        - str: Uses seaborn color palette with this name
        - list/iterable: Uses these colors directly
    figsize : tuple, default=(6, 5)
        Figure dimensions (width, height) in inches.
    title : str, default=""
        Title for the plot.
    legend_title : str, default=""
        Title for the legend.

    Returns:
    --------
    ax : axis object
        Axis object for further customization if needed.
    """
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as PathEffects
    import seaborn as sns
    import collections

    # Handle palette options
    if palette is None:
        # Use default seaborn palette if None is provided
        palette = sns.color_palette("tab10", len(y_cols))
    elif isinstance(palette, str):
        # If palette is a string, use it as a seaborn palette name
        palette = sns.color_palette(palette, len(y_cols))
    elif not isinstance(palette, str) and isinstance(palette, collections.abc.Iterable):
        # If palette is a non-string iterable, subset it to y_cols length
        palette = list(palette)[: len(y_cols)]

    # Create figure and axis
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    plt.grid(False)

    # Extract the x-axis values
    x_vals = df[x]

    # Plot the stacked bars
    bottom = df[y_cols[0]].copy()
    for i, y_col in enumerate(y_cols):
        y_vals = df[y_col]
        if i == 0:
            ax.bar(
                x_vals,
                y_vals,
                label=y_col,
                color=palette[i],
                edgecolor="0.2",
                lw=linewidth,
                width=0.85,
            )
        else:
            ax.bar(
                x_vals,
                y_vals,
                bottom=bottom,
                label=y_col,
                color=palette[i],
                edgecolor="0.2",
                lw=linewidth,
                width=0.85,
            )
            bottom += y_vals

    # Add annotations for cell counts
    if annotate:
        for i, row in df.iterrows():
            cell_state = row[x]

            assert len(y_cols) == len(annot_cols), (
                "y_cols and annot_cols must have the same length"
            )
            bottom = 0
            for j in range(len(y_cols)):
                y_col = y_cols[j]
                annot_col = annot_cols[j]

                text = ax.text(
                    x=cell_state,
                    y=bottom + (row[y_col] / 2),
                    s=f"{row[annot_col]:.0f}",
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=fontsize,
                    # weight="bold",
                    color="white",  # "0.2",
                    rotation=90,
                    ha="center",
                )
                bottom += row[y_col]

    # Set the y-axis limit to 0-1 since we are plotting proportions
    ax.set_ylim(0, 1.001)

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    leg = plt.legend(
        handles,
        legend_labels if legend_labels is not None else labels,
        ncol=1,
        bbox_to_anchor=(1.0, 1.0),
        frameon=False,
        labelcolor="0.2",
        prop={"weight": "bold", "size": fontsize},
        title=legend_title,
        title_fontsize=fontsize,
    )
    leg.get_title().set_weight("bold")

    # Axis ticks
    ax.tick_params(width=linewidth, color="0.2")
    plt.xticks(
        size=fontsize,
        rotation=90,
        # rotation_mode="anchor",
        # ha="right",
        # weight="bold",
        color="0.2",
    )
    plt.yticks(size=fontsize, color="0.2")  # , weight="bold"

    # Add labels and title
    plt.xlabel("", size=fontsize + 2, weight="bold", color="0.2")
    plt.ylabel("", size=fontsize + 2, weight="bold", color="0.2")

    # Axis styling
    n_bars = len(x_vals)
    plt.xlim(-0.60, n_bars - 0.37)  # Tighten x-axis limits to reduce edge spacing

    for axis in ["bottom", "left"]:
        ax.spines[axis].set_linewidth(linewidth)
        ax.spines[axis].set_color("0.2")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.title(
        title,
        size=fontsize + 3,
        weight="bold",
        color="0.2",
        pad=10,
    )

    plt.tight_layout()

    return ax


targets = cell_state_count["target"].unique()
for target in targets:
    target_df = cell_state_count.query(f'target == "{target}"')
    target_df = target_df.sort_values("mut_ratio")  # "cell_state"
    # ax = mvz.plot_stacked_barplot_precalculated(
    ax = plot_stacked_barplot_precalculated(
        target_df,
        x="cell_state",
        y_cols=["mut_ratio", "wt_ratio"],
        annot_cols=["n_cellstate_mut_cells", "n_cellstate_wt_cells"],
        legend_labels=["MUT", "WT"],
        palette=["#081ac2", "#88ccee"],
        fontsize=14,
        linewidth=1.5,
        figsize=(5, 5),
        title=target,
        annotate=True,
    )
    fig = ax.figure
    outdir2 = os.path.join(outdir, "Genotype_Barplots")
    os.makedirs(outdir2, exist_ok=True)
    fig.savefig(
        os.path.join(
            outdir2,
            f"Barplot_MUT_WT_ratio_{target}.pdf",
        ),
        bbox_inches="tight",
        dpi=300,
    )


# ===================================================================================================
# Mutant cell proportion across cell types (CLL vs. LBCL) and targets
# ===================================================================================================
cell_group_col = "cell_type"

# Count #genotyped cells in each cell state
cell_state_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col])
    .size()
    .reset_index()
    .rename(columns={0: "n_cellstate_cells"})
)

# Count #MUT cells in each cell state
cell_state_mut_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col, "best_pred_adj"])
    .size()
    .reset_index()
    .query('best_pred_adj == "MUT"')
    .rename(columns={0: "n_cellstate_mut_cells"})
)
# Merge Total Cell Count and MUT Cell Count
cell_state_count = cell_state_count.merge(
    cell_state_mut_count.drop(columns="best_pred_adj"),
    on=[
        "target",
        "sample",
        cell_group_col,
    ],
    how="left",
).fillna(0)

# Count WT cells in each cell state
cell_state_wt_count = (
    cell_eda_filt[(~cell_eda_filt[cell_group_col].isna()) & (cell_eda_filt["expected"])]
    .groupby(["target", "sample", cell_group_col, "best_pred_adj"])
    .size()
    .reset_index()
    .query('best_pred_adj == "WT"')
    .rename(columns={0: "n_cellstate_wt_cells"})
)
# Merge Total Cell Count and WT Cell Count
cell_state_count = cell_state_count.merge(
    cell_state_wt_count.drop(columns="best_pred_adj"),
    on=[
        "target",
        "sample",
        cell_group_col,
    ],
    how="left",
).fillna(0)

# Calculate the mutation ratio
cell_state_count["mut_ratio"] = (
    cell_state_count["n_cellstate_mut_cells"] / cell_state_count["n_cellstate_cells"]
)
# Calculate the WT ratio
cell_state_count["wt_ratio"] = (
    cell_state_count["n_cellstate_wt_cells"] / cell_state_count["n_cellstate_cells"]
)

cell_state_count.sort_values(["target", "sample", cell_group_col])

# Make targets unique
tmp = cell_state_count.groupby(["target", "sample"]).size().reset_index()
assert tmp["target"].nunique() == tmp.shape[0]
assert (cell_state_count.value_counts(["target", "sample", cell_group_col]) == 1).all()

# Normalize MUT ratio to the total MUT ratio for that target & cell type
mut_ratio_sum = cell_state_count.groupby(["target"])["mut_ratio"].sum().reset_index()
# map mut ratio sum to target-cell type pairs
cell_state_count["mut_ratio_norm"] = cell_state_count.apply(
    lambda x: x["mut_ratio"]
    / mut_ratio_sum[(mut_ratio_sum["target"] == x["target"])]["mut_ratio"].values[0],
    axis=1,
)
cell_state_count.fillna(0, inplace=True)

# Pivot the table to get the cell state count matrix
cell_state_count_matrix = (
    cell_state_count.pivot_table(
        index=["target"],
        columns=cell_group_col,
        values="mut_ratio_norm",
    )
    .reset_index()
    .fillna(0)
)
# Pivot the table to get raw mut_ratio values (for annotations)
cell_state_count_matrix_raw = (
    cell_state_count.pivot_table(
        index=["target"],
        columns=cell_group_col,
        values="mut_ratio",
    )
    .reset_index()
    .fillna(0)
)

data = cell_state_count_matrix.set_index(["target"]).T
data_raw = cell_state_count_matrix_raw.set_index(["target"]).T  # For annotations

target2sample = cell_eda[["target", "expected_sample"]].drop_duplicates()
assert target2sample["target"].nunique() == target2sample.shape[0]
target2sample = target2sample.set_index("target")["expected_sample"].to_dict()

# Create the ordered list of columns based on samples
unique_samples = sorted(
    list(dict.fromkeys([target2sample[col] for col in data.columns]))
)

# # * === Order Input Matrix === * #
# new_order = []
# for sample in unique_samples:
#     sample_targets = [
#         target
#         for target, samp in target2sample.items()
#         if samp == sample and target in data.columns
#     ]
#     new_order.extend(sorted(sample_targets))

# # Reorder the data before creating the clustermap
# data_reordered = data[new_order]
# data_raw_reordered = data_raw[new_order]  # For annotationso


# * === Order Input Matrix (Manual) === * #
col_order = [
    "EMD",
    "B2M",
    "POU2F2_Mut1",
    "MCL1",
    "PLCG2_Mut1",
    "PLCG2_Mut2",
    "PLCG2_Mut3",
    "JAK1",
    "PRKDC",
    "SRRM2",
    "JUNB",
    "IKZF3",
    "SF3B1",
    "IRF8_Mut1",
    # "HIST1H1C",
    "BTK",
]
# reorder the data
data_reordered = data[col_order]
# data_reordered = data_reordered.reindex(row_order)
data_raw_reordered = data_raw[col_order]  # For annotations
# data_raw_reordered = data_raw_reordered.reindex(row_order)  # For annotations


# Create column colors based on the reordered data
sample_order = [target2sample[col] for col in data_reordered.columns]
colors = sns.color_palette("hls", len(unique_samples))
sample_colors = dict(zip(unique_samples, colors))
col_colors = [sample_colors[sample] for sample in sample_order]

fontsize = 16

# Create the clustermap with the reordered data
ax = sns.clustermap(
    data_reordered,  # Use pre-reordered data
    cmap="Blues",
    cbar_pos=(1.15, 0.27, 0.03, 0.23),
    cbar_kws={"label": "Mutation ratio"},
    figsize=(10, 2.4),
    row_cluster=True,
    col_cluster=False,
    dendrogram_ratio=(0.0, 0.00001),
    col_colors=col_colors,
    colors_ratio=0.13,
    col_linkage=None,
    # annot=data_raw_reordered.values,  # Use raw mutation ratios for annotations
    fmt=".2f",
    annot_kws={
        "fontsize": fontsize,
        "weight": "bold",
        "color": "grey",
        "path_effects": [withStroke(linewidth=1, foreground="white")],
    },
)

# Set the Colorbar label font size
cbar = ax.ax_cbar
cbar.set_ylabel(
    cbar.get_ylabel(), fontsize=fontsize + 2, fontweight="bold", labelpad=15
)  # Set font size and weight
cbar.set_yticklabels(
    cbar.get_yticklabels(), fontsize=fontsize, weight="bold"
)  # Set the colorbar tick labels font sizes

# Legend
legend_elements = [
    Patch(facecolor=sample_colors[sample], label=sample) for sample in unique_samples
]
leg = plt.legend(
    handles=legend_elements,
    title="Samples",
    title_fontsize=fontsize + 2,
    bbox_to_anchor=(3.36, 5.3),
    loc=1,
    ncol=1,
    frameon=False,
    labelcolor="0.2",
    prop={"weight": "bold", "size": fontsize},
)
leg.get_title().set_weight("bold")

# Remove x and y labels
ax.ax_heatmap.set_xlabel("Target", fontsize=fontsize + 2, weight="bold", labelpad=10)
ax.ax_heatmap.set_ylabel(
    "Cell State", fontsize=fontsize + 2, weight="bold", labelpad=20
)

# Set the axes tick labels font sizes
ax.ax_heatmap.set_yticklabels(
    ax.ax_heatmap.get_yticklabels(), fontsize=fontsize, weight="bold"
)
ax.ax_heatmap.set_xticklabels(
    ax.ax_heatmap.get_xticklabels(), fontsize=fontsize, weight="bold"
)
# y축 눈금 레이블 회전
for label in ax.ax_heatmap.get_yticklabels():
    label.set_rotation(0)  # 45도 회전
    label.set_ha("left")  # 오른쪽으로 정렬
# x축 눈금 레이블 회전
for label in ax.ax_heatmap.get_xticklabels():
    label.set_rotation(90)  # 45도 회전
    # label.set_ha("right")  # 오른쪽으로 정렬
# Title
ax.ax_heatmap.set_title(
    "Mutation Ratio (CLL vs. Large B-cell)",
    fontsize=16,
    weight="bold",
    pad=20,
    color="0.2",
)
fig = ax.figure
fig.savefig(
    os.path.join(outdir, "Heatmap_MUT_ratio_per_CellType.pdf"),
    bbox_inches="tight",
    dpi=300,
)

plt.show()
