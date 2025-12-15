###################################################################################################
##                        Genotyping Data Visualization                                          ##
###################################################################################################
## This script generates two visualization outputs for GoT-Multi genotyping data:                ##
##  1. UMAP projection colored by genotype (MUT / WT / Unprofiled).                              ##
##  2. Scatter plot comparing per-cell GEN UMI counts vs GEX UMI counts for each target          ##
##                                                                                               ##
## Inputs                                                                                        ##
##  - CellLineMixingGenotypingProfiles.csv: Genotyping (refined) profiles                        ##
##      with columns including: cell_id, target, best_pred_adj (MUT/WT labels).                  ##
##  - CellLineMixingGEX.h5ad: AnnData object with UMAP embedding precomputed                     ##
##      column 'cell_id' in adata.obs linking to the genotyping profile.                         ##
##  - CellLineMixingIronThroneGenotypingProfiles.csv: Genotyping (from IronThrone) profiles      ##
##      with columns: *_calls_*, expected, target.                                               ##
##                                                                                               ##
## Generated Outputs:                                                                            ##
##  - umap_<TARGET>_cell_line_mixing_specificity.pdf : UMAP colored by MUT / WT / Unprofiled     ##
##  - GENumi_vs_GEXumi.pdf : Scatter plot of GEN vs GEX UMI counts + Pearson correlation         ##
##                                                                                               ##
## Usage: Run `python 2_GenotypingDataVisualization.py` after placing inputs in Data/.           ##
###################################################################################################

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import collections

warnings.filterwarnings("ignore")


def _load_inputs(datadir):
    """Internal helper to load prediction dataframe and AnnData object.

    Parameters
    ----------
    datadir : str
        Directory containing required CSV / h5ad inputs.

    Returns
    -------
    pred_df : pandas.DataFrame
        Long-format predictions table with columns at least: cell_id, target, best_pred_adj.
    gex : AnnData
        Expression AnnData containing UMAP embedding and 'cell_id' column in `.obs`.
    """
    pred_df = pd.read_csv(os.path.join(datadir, "CellLineMixingGenotypingProfiles.csv"))
    gex = sc.read_h5ad(os.path.join(datadir, "CellLineMixingGEX.h5ad"))
    return pred_df, gex


# ===================================================================================================
# Visualization: Genotyping on UMAP
# ===================================================================================================
def _annotate_predictions(gex, pred_df):
    """Annotate AnnData with per-target prediction columns.

    Adds a metadata column '<target>_pred' for each unique target labeling cells as
    MUT / WT (from 'best_pred_adj') or 'Unprofiled' if missing.

    Filters cells with no annotation.
    """
    for target in pred_df["target"].unique():
        tmp = pred_df[pred_df["target"] == target]
        gex.obs[f"{target}_pred"] = gex.obs["cell_id"].map(
            tmp.set_index("cell_id")["best_pred_adj"]
        )
        gex.obs[f"{target}_pred"] = gex.obs[f"{target}_pred"].fillna("Unprofiled")
    # Remove cells with no annotation
    gex = gex[gex.obs.isna().sum(axis=1) < 2]
    return gex


def plot_umap(
    adata,
    figsize=(10, 10),
    dpi=150,
    order_by=None,
    order_by_vals=None,
    return_copy=False,
    random_shuffle=True,
    ax=None,
    save_fpath=None,
    **kwargs,
):
    """Plot UMAP embedding with optional ordering and reproducible random shuffling.

    Parameters
    ----------
    adata : AnnData
        AnnData object containing UMAP coordinates in `.obsm['X_umap']` (assumed precomputed).
    figsize : tuple, default (10, 10)
        Figure size (ignored if external `ax` is passed and `color` has >1 elements).
    dpi : int, default 150
        Resolution of the figure when creating a new one or saving.
    order_by : str or None
        Metadata column to group/order cells prior to plotting (e.g., 'is_profiled').
    order_by_vals : list or None
        Explicit ordering of categories in `order_by`; must cover all cells if provided.
    return_copy : bool, default False
        If True returns (ax, reordered_adata); otherwise returns ax only.
    random_shuffle : bool, default True
        Randomly permute plotting order (helps reduce overplot bias).
    ax : matplotlib.axes.Axes or None
        Axis to draw on; if None a new figure/axis is created (unless multi-color list triggers grid).
    save_fpath : str or None
        If provided, path to save figure (PDF/PNG/etc inferred from extension).
    **kwargs : dict
        Additional parameters forwarded to `scanpy.pl.umap` (e.g., color, size, alpha, frameon...).

    Returns
    -------
    matplotlib.axes.Axes or (Axes, AnnData)
        Axis handle; optionally also the reordered AnnData if `return_copy=True`.
    """
    if order_by_vals is not None and order_by is None:
        raise ValueError("order_by must be provided")

    if random_shuffle:
        np.random.seed(1)
        ri = np.random.permutation(list(range(adata.shape[0])))
        adata = adata[ri, :]

        # Shuffle size and alpha if they are list-like
        for key in ["size", "alpha"]:
            if (
                key in kwargs
                and hasattr(kwargs[key], "__len__")
                and len(kwargs[key]) == adata.shape[0]
            ):
                kwargs[key] = np.array(kwargs[key])[ri]

    # order by the given metadata column
    if order_by is not None:
        # in the order of the given values
        if order_by_vals is None:
            order_by_vals = list(adata.obs[order_by].unique())

        idx_list = []
        for order_val in order_by_vals:
            idx_list.append(adata.obs[adata.obs[order_by] == order_val].index.tolist())
        assert sum([len(x) for x in idx_list]) == adata.shape[0], (
            "All values in the `order_by` column must be included"
        )
        order = [i for idxs in idx_list for i in idxs]
        order_indices = [
            adata.obs.index.get_loc(i) for i in order
        ]  # Get positional indices
        adata = adata[order, :]

        # Reorder size and alpha if they are list-like
        for key in ["size", "alpha"]:
            if (
                key in kwargs
                and hasattr(kwargs[key], "__len__")
                and len(kwargs[key]) == adata.shape[0]
            ):
                kwargs[key] = np.array(kwargs[key])[order_indices]

    if (
        kwargs.get("color") is not None
        and isinstance(kwargs["color"], collections.abc.Iterable)
        and len(kwargs["color"]) > 1
    ):
        ax = None
    else:
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
        else:
            fig = ax.get_figure()

    ax = sc.pl.umap(adata, **kwargs, ax=ax, show=False)

    if save_fpath is not None:
        fig.savefig(save_fpath, dpi=dpi, bbox_inches="tight")

    if return_copy:
        return ax, adata
    return ax


def _plot_per_target_umaps(gex, pred_df, outdir):
    """Iterate over targets to produce UMAP."""
    for target in pred_df["target"].unique():
        # (a) Binary profiling indicator for ordering layer priority
        gex.obs["is_profiled"] = gex.obs[f"{target}_pred"] != "Unprofiled"

        # (b) Stable ordering: unprofiled first, profiled randomized internally for visual dispersion
        order_by_vals = [False, True]
        idx_list = []
        for order_val in order_by_vals:
            if order_val:  # randomize profiled
                idx_list.append(
                    gex.obs[gex.obs["is_profiled"] == order_val]
                    .sample(frac=1, random_state=7)
                    .index
                )
            else:
                idx_list.append(gex.obs[gex.obs["is_profiled"] == order_val].index)
        assert sum(len(idx) for idx in idx_list) == gex.shape[0]
        gex = gex[pd.Index(np.concatenate(idx_list))]

        # (c) Draw & save
        ax = plot_umap(
            gex,
            color=[f"{target}_pred"],
            figsize=(5, 5),
            size=60,
            palette={"MUT": "#081ac2", "WT": "#8ce", "Unprofiled": "lightgrey"},
            alpha=1.0,
            dpi=300,
            frameon=False,
            random_shuffle=False,
        )
        ax.set_title(target, size=12, weight="bold", color="0.2", pad=1)
        plt.legend(
            bbox_to_anchor=(1.35, 0.75),
            frameon=False,
            fontsize=14,
            labelcolor="0.2",
            prop={"weight": "bold"},
        )
        plt.savefig(
            os.path.join(outdir, f"umap_{target}_cell_line_mixing_specificity.pdf"),
            bbox_inches="tight",
            dpi=300,
        )


# ===================================================================================================
# Visualization: GEN UMI Expression vs GEX UMI Expression
# ===================================================================================================
def _plot_gen_vs_gex(datadir, gex, outdir):
    """Generate scatter plot comparing GEN vs GEX UMI counts with correlation statistics."""
    gen_df = pd.read_csv(
        os.path.join(
            datadir,
            "CellLineMixingIronThroneGenotypingProfiles.csv",
        ),
        index_col=0,
    )

    # Total genotyping UMIs per cell: sum of wildcard + mutant contributions
    gen_df["gen_umi_count"] = (
        gen_df["wt_calls_geneseq"]
        + gen_df["wt_calls_probebc"]
        + gen_df["mut_calls_geneseq"]
        + gen_df["mut_calls_probebc"]
    )

    # Dense conversion (warn if large)
    if gex.n_obs * gex.n_vars > 5e7:
        print("[WARN] Converting a large matrix to dense; memory usage may be high.")
    count_df = pd.DataFrame(gex.X.toarray(), index=gex.obs.index, columns=gex.var.index)
    gen_df["gex_umi_count"] = gen_df.apply(
        lambda row: count_df.loc[row.name, row["target"]], axis=1
    )

    from scipy.stats import pearsonr  # local import to keep top-level lean

    df = gen_df[gen_df["expected"] == 1]
    x = "gen_umi_count"
    y = "gex_umi_count"
    hue = "target"

    # Randomize plotting order for aesthetic
    df = df.sample(frac=1, random_state=42).reset_index(drop=True)
    plt.figure(figsize=(5, 5))
    sns.scatterplot(
        df,
        x=x,
        y=y,
        hue=hue,
        alpha=0.6,
        edgecolor="none",
        palette={"HIST1H1C": "#154c79", "S100A10": "#e28743"},
    )
    ax = plt.gca()
    for axis in ["bottom", "left", "top", "right"]:
        ax.spines[axis].set_linewidth(2)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(width=2)
    plt.xticks(size=10, weight="bold")
    plt.yticks(size=10, weight="bold")
    plt.legend(loc=1, bbox_to_anchor=(1.25, 1), frameon=True, prop={"weight": "bold"})

    pcc, p = pearsonr(df[x], df[y])
    plt.title(
        f"GEN UMI Count vs GEX UMI Count \n(PCC: {pcc:.4f})",
        weight="bold",
        color="0.2",
        size=15,
    )
    plt.xlabel("#GEN UMI count", weight="bold", color="0.2", size=15)
    plt.ylabel("#GEX UMI count", weight="bold", color="0.2", size=15)
    plt.savefig(
        os.path.join(outdir, "GENumi_vs_GEXumi.pdf"),
        bbox_inches="tight",
        dpi=300,
    )


# ===================================================================================================
# Main Execution (guard for importability)
# ===================================================================================================
if __name__ == "__main__":
    datadir = "Data/"
    outdir = "Results/"
    os.makedirs(outdir, exist_ok=True)

    # Step 1. Load inputs
    pred_df, gex = _load_inputs(datadir)

    # Step 2. Annotate predictions
    gex = _annotate_predictions(gex, pred_df)

    # Step 3 & 4. Per-target UMAPs
    _plot_per_target_umaps(gex, pred_df, outdir)

    # Step 5. GEN vs GEX scatter & correlation
    _plot_gen_vs_gex(datadir, gex, outdir)
    print(f"Outputs written to: {outdir}")



# ===================================================================================================
# Mutant frequency via UMI downsampling analysis
# ===================================================================================================
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Input: GoT-Multi-ML genotyping profile


gen = pd.read_csv(
    os.path.join(inputdir, "ML_genotyping/aggregated/final_prediction_results.csv")
)
targets = gen_pred["target"].unique()

# TODO: === Downsampling to a single amplicon UMI per cell === #
def perform_single_amplicon_umi_downsampling_analysis(
    data: pd.DataFrame,
    cells_to_sample_per_sample: int,
    num_iterations: int = 100,
    cluster_column: str = "cluster_id",
    analyze_per_sample: bool = True,
) -> pd.DataFrame:
    """
    Performs the "Downsampling to a single amplicon UMI per cell" analysis
    as described in the GoT paper to assess robustness of mutant cell frequencies.

    Args:
        data (pd.DataFrame): Input DataFrame with columns:
                             'cell_id', 'sample', 'target', 'wt_calls', 'mut_calls',
                             and a column for cell clusters (specified by cluster_column).
                             'wt_calls' and 'mut_calls' are the number of wild-type and mutant UMIs.
        cells_to_sample_per_sample (int): The number of cells to randomly subsample
                                         from each unique 'sample' (patient).
        num_iterations (int): The number of times to repeat the downsampling process.
        cluster_column (str): Name of the cluster column.
        analyze_per_sample (bool): If True, analyze per sample; otherwise integrated.

    Returns:
        (final_summary, results_df)
    """

    results = []
    if cluster_column not in data.columns:
        raise ValueError(f"'{cluster_column}' not found in data columns")

    print(f"Starting downsampling analysis with {num_iterations} iterations...")
    print(
        f"Analysis mode: {'Per sample' if analyze_per_sample else 'Integrated across samples'}"
    )

    samples_to_analyze = data["sample"].unique() if analyze_per_sample else [None]

    for sample_name in samples_to_analyze:
        print(
            f"Processing sample: {sample_name if sample_name else 'All samples (integrated)'}"
        )
        for iteration in range(num_iterations):
            if analyze_per_sample:
                sample_data = data[data["sample"] == sample_name]
                n_cells = min(cells_to_sample_per_sample, len(sample_data))
                if n_cells == 0:
                    continue
                sampled_data = sample_data.sample(
                    n=n_cells, random_state=iteration, replace=False
                )
            else:
                sampled_cells_list = []
                for samp in data["sample"].unique():
                    sample_data = data[data["sample"] == samp]
                    n_cells = min(cells_to_sample_per_sample, len(sample_data))
                    if n_cells > 0:
                        sampled_cells_list.append(
                            sample_data.sample(
                                n=n_cells, random_state=iteration, replace=False
                            )
                        )
                if not sampled_cells_list:
                    continue
                sampled_data = pd.concat(sampled_cells_list)

            sampled_data = sampled_data.copy()
            sampled_data["genotype_call"] = "NA"

            for idx, row in sampled_data.iterrows():
                total_umis = row["wt_calls"] + row["mut_calls"]
                if total_umis == 0:
                    continue
                elif total_umis == 1:
                    sampled_data.loc[idx, "genotype_call"] = (
                        "mutant" if row["mut_calls"] == 1 else "wildtype"
                    )
                else:
                    prob_mut = row["mut_calls"] / total_umis
                    np.random.seed(iteration * 10000 + idx)
                    sampled_data.loc[idx, "genotype_call"] = (
                        "mutant" if np.random.rand() < prob_mut else "wildtype"
                    )

            genotyped_cells = sampled_data[
                sampled_data["genotype_call"].isin(["mutant", "wildtype"])
            ].copy()
            if genotyped_cells.empty:
                continue

            cluster_frequencies = {}
            for cluster in genotyped_cells[cluster_column].unique():
                cluster_data = genotyped_cells[
                    genotyped_cells[cluster_column] == cluster
                ]
                mutant_cells_in_cluster = (
                    cluster_data["genotype_call"] == "mutant"
                ).sum()
                total_genotyped_cells_in_cluster = len(cluster_data)
                cluster_frequencies[cluster] = (
                    (mutant_cells_in_cluster / total_genotyped_cells_in_cluster)
                    if total_genotyped_cells_in_cluster > 0
                    else 0
                )

            total_freq_sum = sum(cluster_frequencies.values())
            normalized_cluster_frequencies = {
                c: (f / total_freq_sum if total_freq_sum > 0 else 0)
                for c, f in cluster_frequencies.items()
            }

            for cluster, norm_freq in normalized_cluster_frequencies.items():
                result_entry = {
                    "iteration": iteration,
                    "cluster_id": cluster,
                    "normalized_mutant_frequency": norm_freq,
                    "mutant_frequency": cluster_frequencies[cluster],
                }
                if analyze_per_sample:
                    result_entry["sample"] = sample_name
                results.append(result_entry)

    if not results:
        print("No valid results accumulated across all iterations.")
        return pd.DataFrame()

    results_df = pd.DataFrame(results)
    grouping_cols = ["sample", "cluster_id"] if analyze_per_sample else ["cluster_id"]
    final_summary = (
        results_df.groupby(grouping_cols)
        .agg(
            {
                "normalized_mutant_frequency": ["mean", "std"],
                "mutant_frequency": ["mean", "std"],
            }
        )
        .reset_index()
    )
    final_summary.columns = grouping_cols + [
        "mean_normalized_freq",
        "std_normalized_freq",
        "mean_mutant_freq",
        "std_mutant_freq",
    ]
    final_summary["std_normalized_freq"] = final_summary["std_normalized_freq"].fillna(
        0
    )
    final_summary["std_mutant_freq"] = final_summary["std_mutant_freq"].fillna(0)
    print("Analysis complete.")
    return final_summary, results_df


def calculate_actual_mutant_frequencies(
    data: pd.DataFrame,
    cluster_column: str = "cell_state",
    genotype_column: str = "best_pred_adj",
) -> pd.DataFrame:
    genotyped_data = data[data[genotype_column].isin(["MUT", "WT"])].copy()
    results = []
    for sample in genotyped_data["sample"].unique():
        sample_data = genotyped_data[genotyped_data["sample"] == sample]
        for cluster in sample_data[cluster_column].unique():
            cluster_data = sample_data[sample_data[cluster_column] == cluster]
            total_cells = len(cluster_data)
            mutant_cells = (cluster_data[genotype_column] == "MUT").sum()
            mutant_freq = (mutant_cells / total_cells) if total_cells > 0 else 0
            results.append(
                {
                    "sample": sample,
                    "cluster_id": cluster,
                    "total_cells": total_cells,
                    "mutant_cells": mutant_cells,
                    "actual_mutant_freq": mutant_freq,
                }
            )
    results_df = pd.DataFrame(results)
    for sample in results_df["sample"].unique():
        sample_data = results_df[results_df["sample"] == sample]
        total_freq_sum = sample_data["actual_mutant_freq"].sum()
        if total_freq_sum > 0:
            results_df.loc[results_df["sample"] == sample, "actual_normalized_freq"] = (
                results_df.loc[results_df["sample"] == sample, "actual_mutant_freq"]
                / total_freq_sum
            )
        else:
            results_df.loc[results_df["sample"] == sample, "actual_normalized_freq"] = 0
    return results_df


def plot_target_cellstate_heatmaps(
    all_actual_frequencies: dict,
    all_downsampled_frequencies: dict,
    title: str = "Mutant Cell Frequencies Across Targets and Cell States",
    figsize: tuple = (20, 8),
    save_path: str = None,
    use_normalized: bool = True,
    annot: bool = True,
) -> None:
    plt.style.use("default")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    if use_normalized:
        actual_freq_col = "actual_normalized_freq"
        downsampled_freq_col = "mean_normalized_freq"
        color_label = "Normalized Mutant\nCell Frequency"
    else:
        actual_freq_col = "actual_mutant_freq"
        downsampled_freq_col = "mean_mutant_freq"
        color_label = "Mutant Cell Frequency"
    actual_data_list, downsampled_data_list = [], []
    for target, actual_df in all_actual_frequencies.items():
        if target in all_downsampled_frequencies:
            downsampled_df = all_downsampled_frequencies[target]
            a = actual_df.copy()
            a["target"] = target
            d = downsampled_df.copy()
            d["target"] = target
            actual_data_list.append(a)
            downsampled_data_list.append(d)
    if not actual_data_list or not downsampled_data_list:
        print("No data available for heatmap creation.")
        return
    all_actual = pd.concat(actual_data_list, ignore_index=True)
    all_downsampled = pd.concat(downsampled_data_list, ignore_index=True)
    actual_avg = (
        all_actual.groupby(["target", "cluster_id"])[actual_freq_col]
        .mean()
        .reset_index()
    )
    downsampled_avg = (
        all_downsampled.groupby(["target", "cluster_id"])[downsampled_freq_col]
        .mean()
        .reset_index()
    )
    all_targets = sorted(set(actual_avg["target"]) | set(downsampled_avg["target"]))
    all_cell_states = sorted(
        set(actual_avg["cluster_id"]) | set(downsampled_avg["cluster_id"])
    )
    actual_matrix = pd.DataFrame(
        index=all_targets, columns=all_cell_states, dtype=float
    )
    downsampled_matrix = pd.DataFrame(
        index=all_targets, columns=all_cell_states, dtype=float
    )
    for _, row in actual_avg.iterrows():
        actual_matrix.loc[row["target"], row["cluster_id"]] = row[actual_freq_col]
    for _, row in downsampled_avg.iterrows():
        downsampled_matrix.loc[row["target"], row["cluster_id"]] = row[
            downsampled_freq_col
        ]
    actual_matrix = actual_matrix.fillna(0)
    downsampled_matrix = downsampled_matrix.fillna(0)
    all_values = pd.concat([actual_matrix.stack(), downsampled_matrix.stack()])
    vmin, vmax = all_values.min(), all_values.max()
    sns.heatmap(
        actual_matrix,
        ax=ax1,
        cmap="viridis",
        cbar_kws={"label": color_label},
        vmin=vmin,
        vmax=vmax,
        annot=annot,
        fmt=".3f",
        linewidths=0,
        linecolor="white",
    )
    ax1.set_title("Actual Genotype Calls", fontsize=14, fontweight="bold")
    ax1.set_xlabel("Cell State", fontsize=12, fontweight="bold")
    ax1.set_ylabel("Target", fontsize=12, fontweight="bold")
    ax1.tick_params(axis="x", rotation=45)
    ax1.tick_params(axis="y", rotation=0)
    sns.heatmap(
        downsampled_matrix,
        ax=ax2,
        cmap="viridis",
        cbar_kws={"label": color_label},
        vmin=vmin,
        vmax=vmax,
        annot=annot,
        fmt=".3f",
        linewidths=0,
        linecolor="white",
    )
    ax2.set_title("Downsampled to Single UMI", fontsize=14, fontweight="bold")
    ax2.set_xlabel("Cell State", fontsize=12, fontweight="bold")
    ax2.set_ylabel("Target", fontsize=12, fontweight="bold")
    ax2.tick_params(axis="x", rotation=45)
    ax2.tick_params(axis="y", rotation=0)
    plt.suptitle(title, fontsize=16, fontweight="bold", y=1.02)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Target-cellstate heatmaps saved to: {save_path}")
    plt.show()
    print("\nHeatmap Summary:")
    print(f"Targets: {len(all_targets)}")
    print(f"Cell states: {len(all_cell_states)}")
    print(
        f"Total target-cellstate combinations: {len(all_targets) * len(all_cell_states)}"
    )
    print(f"Combinations with actual data: {(actual_matrix > 0).sum().sum()}")
    print(f"Combinations with downsampled data: {(downsampled_matrix > 0).sum().sum()}")
    return actual_matrix, downsampled_matrix


def plot_comprehensive_scatter_comparison(
    all_actual_frequencies: dict,
    all_downsampled_frequencies: dict,
    title: str = "Comprehensive Comparison: Actual vs Downsampled Frequencies Across All Targets",
    figsize: tuple = (14, 10),
    save_path: str = None,
    use_normalized: bool = True,
    show_error_bars: bool = True,
    errorbar_capsize: float = 3,
    y_error_type: str = "iteration",  # 'iteration', 'sample', or 'both'
    iteration_sd_aggregation: str = "mean",  # 'mean' or 'rms'
    error_metric: str = "se",  # 'sd' or 'se' for iteration component
    iteration_count: int = 100,
) -> None:
    """
    Scatter comparing actual vs downsampled frequencies.

    error_metric:
        'sd' uses standard deviation.
        'se' uses standard error = sd / sqrt(iteration_count) for iteration component.
    """
    plt.style.use("default")
    fig, ax = plt.subplots(figsize=figsize)
    if use_normalized:
        actual_freq_col = "actual_normalized_freq"
        downsampled_freq_col = "mean_normalized_freq"
        downsampled_iter_sd_col = "std_normalized_freq"
        axis_label = "Normalized Mutant Cell Frequency"
    else:
        actual_freq_col = "actual_mutant_freq"
        downsampled_freq_col = "mean_mutant_freq"
        downsampled_iter_sd_col = "std_mutant_freq"
        axis_label = "Mutant Cell Frequency"
    comparison_data_list = []
    # Build aggregated dataframe
    for target, actual_df in all_actual_frequencies.items():
        if target not in all_downsampled_frequencies:
            continue
        downsampled_df = all_downsampled_frequencies[target]
        actual_stats = (
            actual_df.groupby("cluster_id")[actual_freq_col]
            .agg(["mean", "std"])
            .reset_index()
            .rename(
                columns={
                    "mean": f"{actual_freq_col}_mean",
                    "std": f"{actual_freq_col}_sample_std",
                }
            )
        )
        downsampled_stats = (
            downsampled_df.groupby("cluster_id")[downsampled_freq_col]
            .agg(["mean", "std"])
            .reset_index()
            .rename(
                columns={
                    "mean": f"{downsampled_freq_col}_mean",
                    "std": f"{downsampled_freq_col}_sample_std",
                }
            )
        )
        # iteration SD per sample
        iter_component = (
            downsampled_df.groupby(["cluster_id", "sample"])[downsampled_iter_sd_col]
            .first()
            .reset_index()
        )
        # Convert to SE if requested (per-sample)
        if error_metric == "se":
            iter_component[downsampled_iter_sd_col] = iter_component[
                downsampled_iter_sd_col
            ] / np.sqrt(iteration_count)

        def agg_iter(x):
            if iteration_sd_aggregation == "rms":
                return np.sqrt(np.mean(x**2))
            else:
                return np.mean(x)

        iter_sd = (
            iter_component.groupby("cluster_id")[downsampled_iter_sd_col]
            .apply(agg_iter)
            .reset_index()
            .rename(
                columns={
                    downsampled_iter_sd_col: f"{downsampled_freq_col}_iteration_var"
                }
            )
        )
        # Rename iteration_var suffix to indicate metric
        metric_suffix = "se" if error_metric == "se" else "sd"
        iter_sd = iter_sd.rename(
            columns={
                f"{downsampled_freq_col}_iteration_var": f"{downsampled_freq_col}_iteration_{metric_suffix}"
            }
        )
        for df_ in (actual_stats, downsampled_stats, iter_sd):
            for col in df_.columns:
                if any(
                    col.endswith(suf)
                    for suf in ["_std", "_sample_std", "_iteration_sd", "_iteration_se"]
                ):
                    df_[col] = df_[col].fillna(0)
        merged = actual_stats.merge(
            downsampled_stats, on="cluster_id", how="inner"
        ).merge(iter_sd, on="cluster_id", how="left")
        merged["target"] = target
        comparison_data_list.append(merged)
    if not comparison_data_list:
        print("No data available for comprehensive scatter plot.")
        return
    all_comparison_data = pd.concat(comparison_data_list, ignore_index=True)
    unique_targets = sorted(all_comparison_data["target"].unique())
    unique_cell_states = sorted(all_comparison_data["cluster_id"].unique())
    target_colors = dict(
        zip(unique_targets, plt.cm.tab10(np.linspace(0, 1, len(unique_targets))))
    )
    markers = [
        "o",
        "s",
        "^",
        "D",
        "v",
        "<",
        ">",
        "p",
        "*",
        "h",
        "H",
        "+",
        "x",
        "|",
        "_",
    ]
    while len(markers) < len(unique_cell_states):
        markers.extend(markers)
    cell_state_marker_map = dict(
        zip(unique_cell_states, markers[: len(unique_cell_states)])
    )
    actual_mean_col = f"{actual_freq_col}_mean"
    down_mean_col = f"{downsampled_freq_col}_mean"
    down_iteration_col = (
        f"{downsampled_freq_col}_iteration_{'se' if error_metric == 'se' else 'sd'}"
    )
    down_sample_std_col = f"{downsampled_freq_col}_sample_std"
    for target in unique_targets:
        for cell_state in unique_cell_states:
            subset = all_comparison_data[
                (all_comparison_data["target"] == target)
                & (all_comparison_data["cluster_id"] == cell_state)
            ]
            if subset.empty:
                continue
            x_val = subset[actual_mean_col].iloc[0]
            y_val = subset[down_mean_col].iloc[0]
            y_err_iteration = (
                subset[down_iteration_col].iloc[0]
                if down_iteration_col in subset
                else 0
            )
            y_err_sample = (
                subset[down_sample_std_col].iloc[0]
                if down_sample_std_col in subset
                else 0
            )
            if y_error_type == "iteration":
                y_err = y_err_iteration
            elif y_error_type == "sample":
                y_err = y_err_sample
            elif y_error_type == "both":
                y_err = np.sqrt(y_err_iteration**2 + y_err_sample**2)
            else:
                y_err = 0
            ax.scatter(
                x_val,
                y_val,
                c=[target_colors[target]],
                marker=cell_state_marker_map[cell_state],
                s=100,
                alpha=0.85,
                edgecolors="black",
                linewidth=1,
            )
            if show_error_bars and y_err > 0:
                ax.errorbar(
                    x_val,
                    y_val,
                    yerr=y_err,
                    ecolor=target_colors[target],
                    elinewidth=1.2,
                    capsize=errorbar_capsize,
                    alpha=0.75,
                )
    min_val = min(
        all_comparison_data[actual_mean_col].min(),
        all_comparison_data[down_mean_col].min(),
    )
    max_val = max(
        all_comparison_data[actual_mean_col].max(),
        all_comparison_data[down_mean_col].max(),
    )
    ax.plot([min_val, max_val], [min_val, max_val], "k--", alpha=0.7, linewidth=2)
    correlation = all_comparison_data[actual_mean_col].corr(
        all_comparison_data[down_mean_col]
    )
    ax.text(
        0.05,
        0.95,
        f"Overall Correlation: {correlation:.3f}",
        transform=ax.transAxes,
        fontsize=12,
        fontweight="bold",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )
    metric_label = "SE" if error_metric == "se" else "SD"
    ax.set_xlabel(f"Actual {axis_label} (Mean)", fontsize=12, fontweight="bold")
    ax.set_ylabel(
        f"Downsampled {axis_label} (Mean ± {metric_label})",
        fontsize=12,
        fontweight="bold",
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=20)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3)
    target_legend_elements = [
        plt.scatter(
            [],
            [],
            c=[target_colors[t]],
            marker="o",
            s=100,
            edgecolors="black",
            label=t,
            alpha=0.7,
        )
        for t in unique_targets
    ]
    cell_state_legend_elements = [
        plt.scatter(
            [],
            [],
            c="gray",
            marker=cell_state_marker_map[c],
            s=100,
            edgecolors="black",
            label=c,
            alpha=0.7,
        )
        for c in unique_cell_states
    ]
    target_legend = ax.legend(
        handles=target_legend_elements,
        title="Targets",
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
    )
    target_legend.set_title("Targets", prop={"weight": "bold"})
    cell_state_legend = ax.legend(
        handles=cell_state_legend_elements,
        title="Cell States",
        bbox_to_anchor=(1.05, 0.3),
        loc="center left",
    )
    cell_state_legend.set_title("Cell States", prop={"weight": "bold"})
    ax.add_artist(target_legend)
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        print(f"Comprehensive scatter comparison saved to: {save_path}")
    plt.show()
    print("\nComprehensive Scatter Plot Summary:")
    print(f"Total data points: {len(all_comparison_data)}")
    print(f"Targets: {len(unique_targets)} - {unique_targets}")
    print(f"Cell states: {len(unique_cell_states)} - {unique_cell_states}")
    print(f"Overall correlation: {correlation:.3f}")
    return all_comparison_data


# ========================== #
# ====== Run Analysis ====== #
# ========================== #
print("=== Available Targets ===")
available_targets = gen["target"].unique()
print(f"Number of targets: {len(available_targets)}")
print(f"Targets: {available_targets}")
num_cells_to_sample = 500
num_iterations = 100
all_actual_frequencies = {}
all_downsampled_frequencies = {}
for target_name in available_targets:
    print(f"\n{'=' * 60}")
    print(f"ANALYZING TARGET: {target_name}")
    print(f"{'=' * 60}")
    gen_df_target = gen[gen["target"] == target_name][
        ["cell_id", "sample", "target", "wt_calls", "mut_calls", "cell_state"]
    ]
    gen_orig_target = gen_orig[gen_orig["target"] == target_name]
    print(f"Data for {target_name}:")
    print(f"  - Total cells: {len(gen_df_target)}")
    print(f"  - Samples: {gen_df_target['sample'].unique()}")
    print(f"  - Cell states: {gen_df_target['cell_state'].unique()}")
    if len(gen_df_target) < 10:
        print(
            f"Skipping {target_name} - insufficient data (only {len(gen_df_target)} cells)"
        )
        continue
    target_outdir = os.path.join(outdir, f"target_{target_name}")
    createFolder(target_outdir)
    print(f"\n=== Running per-sample analysis for {target_name} ===")
    mut_frequencies_per_sample, results_df_per_sample = (
        perform_single_amplicon_umi_downsampling_analysis(
            data=gen_df_target,
            cluster_column="cell_state",
            cells_to_sample_per_sample=num_cells_to_sample,
            num_iterations=num_iterations,
            analyze_per_sample=True,
        )
    )
    print(f"\n--- Per-Sample Results for {target_name} ---")
    print(mut_frequencies_per_sample)
    print(f"\n=== Running integrated analysis for {target_name} ===")
    mut_frequencies_integrated, results_df_integrated = (
        perform_single_amplicon_umi_downsampling_analysis(
            data=gen_df_target,
            cluster_column="cell_state",
            cells_to_sample_per_sample=num_cells_to_sample,
            num_iterations=num_iterations,
            analyze_per_sample=False,
        )
    )
    print(f"\n--- Integrated Results for {target_name} ---")
    print(mut_frequencies_integrated)
    print(f"\n=== Calculating actual mutant frequencies for {target_name} ===")
    actual_frequencies = calculate_actual_mutant_frequencies(
        data=gen_orig_target,
        cluster_column="cell_state",
        genotype_column="best_pred_adj",
    )
    print(f"\n--- Actual Mutant Frequencies for {target_name} ---")
    print(actual_frequencies)
    all_actual_frequencies[target_name] = actual_frequencies
    all_downsampled_frequencies[target_name] = mut_frequencies_per_sample
    print(f"\nCompleted analysis for {target_name}")
    print(f"Results saved in: {target_outdir}")

