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
