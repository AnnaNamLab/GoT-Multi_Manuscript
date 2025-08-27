###################################################################################################
##                         scRNA-seq Doublet Detection & QC (Python)                             ##
###################################################################################################
## This script provides utility functions to:
##  A. Perform doublet detection on a raw count AnnData object using multiple methods (SOLO /
##     scVI, DoubletDetection, and Scrublet).
##  B. Add common quality-control (QC) metrics (mitochondrial, ribosomal, hemoglobin genes, etc.).
##  C. Flag outlier cells based on Median Absolute Deviation (MAD) thresholds without removing them
##     (except for a coarse high-mitochondrial-content filter and optional hard filters).
##
## Workflow Overview:
##  Step 1. (Optional) Multi-method doublet detection on raw counts.
##  Step 2. (Optional) Remove definite very-low-quality cells (e.g., low gene count, high mitochondrial percentage >25%) prior to QC.
##  Step 3. Add QC metrics (mitochondrial, ribosomal and optionally hemoglobin gene sets).
##  Step 4. Annotate potential outlier cells using MAD-based thresholds (label only, no removal).
##
## Expected Input:
##  - An AnnData object whose `.X` holds raw UMI counts (not log-transformed).
##  - For SOLO / scVI you must have `scvi` installed and GPU is optional but recommended.
##
## Outputs / Results:
##  - Adds boolean columns in `adata.obs` for each doublet method used: `doublet_solo`,
##    `doublet_dd`, `doublet_scrublet` plus an intersected consensus column `doublet_intersected`.
##  - Adds QC metrics (as produced by `scanpy.pp.calculate_qc_metrics`).
##  - Adds an `Outlier` annotation (boolean) for potential QC outliers (not removed).
##  - Filters high mitochondrial cells (> 25%).
##
###################################################################################################

"""Utility functions for doublet detection and QC on scRNA-seq AnnData objects."""

import os
import doubletdetection
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
from copy import deepcopy
from scipy.stats import median_abs_deviation as mad


# ====== QC, Preprocess ====== #
def add_qc(adata, hb_genes=False):
    """Add standard QC metrics (mitochondrial, ribosomal, optional hemoglobin) to an AnnData object.

    Parameters
    ----------
    adata : AnnData
        AnnData object with raw counts in `.X`.
    hb_genes : bool, default False
        Whether to detect and include hemoglobin genes in QC metrics.

    Returns
    -------
    AnnData
        Updated AnnData object (a shallow copy after MT gene removal) with:
          - Columns added to `adata.var`: `mt`, `ribo` (and optionally `hb`).
          - QC summary metrics in `adata.obs` from `scanpy.pp.calculate_qc_metrics`.
          - Mitochondrial counts retained in `adata.obsm['MT']` for reference.
    """
    # === Identify mitochondrial genes (prefix MT-) === #
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # ! make sure capitalization is correct!
    print("#mt genes: ", adata.var["mt"].sum())

    # === Identify ribosomal genes (download KEGG_RIBOSOME if online, else fallback local file) === #
    try:
        ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
        ribo_genes = pd.read_table(ribo_url, skiprows=2, header=None)
    except Exception as e:
        ribo_genes_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "ribo_genes.txt"
        )
        ribo_genes = pd.read_csv(ribo_genes_path, header=None)
        print(
            f"[WARN] Online KEGG_RIBOSOME download failed ({e}); using local file {ribo_genes_path}."
        )
    adata.var["ribo"] = adata.var_names.isin(ribo_genes[0].values)
    print("#ribo genes: ", adata.var["ribo"].sum())

    qc_var_list = ["mt", "ribo"]
    qc_var_remove = [
        "total_counts_mt",
        "log1p_total_counts_mt",
        "total_counts_ribo",
        "log1p_total_counts_ribo",
    ]

    if hb_genes:
        # === Identify hemoglobin genes (pattern HB but excluding pseudogenes HBxP) === #
        adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]"))
        print("#hb genes: ", adata.var["hb"].sum())
        qc_var_list.append("hb")
        qc_var_remove.extend(["total_counts_hb", "log1p_total_counts_hb"])

    # === Calculate QC metrics (adds pct_counts_mt etc.) === #
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=qc_var_list,
        percent_top=[20],
        log1p=True,
        inplace=True,
    )
    adata.obs = adata.obs[[x for x in adata.obs.columns if x not in qc_var_remove]]

    # === Store raw MT counts & remove MT genes from feature space (retain info in obsm) === #
    adata.obsm["MT"] = adata[:, adata.var["mt"].values].X.toarray()
    adata = adata[:, ~adata.var["mt"].values].copy()

    return adata


def mad_outlier(adata, metric, nmads, upper_only=False):
    """Return a boolean vector denoting outliers based on Median Absolute Deviation (MAD).

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    metric : str
        Column in `adata.obs` to test (e.g., 'log1p_total_counts').
    nmads : int or float
        Number of MADs from median to define an outlier boundary.
    upper_only : bool, default False
        If True, only flag very large values (right tail). If False, flag both tails.

    Returns
    -------
    pandas.Series[bool]
        True for outlier cells.
    """
    M = adata.obs[metric]

    if not upper_only:
        return (M < np.median(M) - nmads * mad(M)) | (M > np.median(M) + nmads * mad(M))

    return M > np.median(M) + nmads * mad(M)


def pp(
    adata,
    doublet_detection=False,
    solo=True,
    dd=True,
    scrublet=True,
    hb_genes=False,
    remove_definite_outliers=True,
):
    """Preprocess an AnnData object: optional doublet detection + QC + outlier annotation.

    Parameters
    ----------
    adata : AnnData
        Input AnnData with raw counts in `.X`.
    doublet_detection : bool, default False
        Whether to run any doublet detection methods.
    solo : bool, default True
        Use scVI + SOLO model for doublet classification (requires scvi-tools).
    dd : bool, default True
        Use the `doubletdetection` BoostClassifier method.
    scrublet : bool, default True
        Use Scanpy's wrapper for Scrublet.
    hb_genes : bool, default False
        Include hemoglobin genes in QC if True.
    remove_definite_outliers : bool, default True
        Apply coarse filters (e.g. min_genes) before computing QC metrics.

    Returns
    -------
    AnnData
        Processed AnnData with:
          - Doublet columns (if enabled) and consensus `doublet_intersected`.
          - QC metrics & outlier flag `Outlier`.
          - Mitochondrial genes removed from feature space but counts stored in `obsm['MT']`.
    """
    # ===================================================================================================
    # Step 1. Doublet Detection (optional): run selected algorithms and create consensus annotation
    # ===================================================================================================
    if doublet_detection:
        assert adata.X.max().is_integer(), "adata.X must be raw counts"
        if solo:
            adata_dd = deepcopy(adata)
            # === DD Method1: SCVI SOLO === #
            sc.pp.filter_genes(adata_dd, min_cells=3)
            sc.pp.highly_variable_genes(
                adata_dd, n_top_genes=2000, subset=True, flavor="seurat_v3"
            )
            scvi.model.SCVI.setup_anndata(adata_dd, layer="counts")
            vae = scvi.model.SCVI(adata_dd)
            vae.train(batch_size=512, early_stopping=True)
            solo = scvi.external.SOLO.from_scvi_model(vae)
            solo.train(batch_size=512)
            df = solo.predict()
            df["prediction"] = solo.predict(soft=False)
            # df.index = df.index.map(lambda x: x[:-2])
            df["dif"] = df.doublet - df.singlet
            doublets = df[(df.prediction == "doublet") & (df.dif > 0.2)]

            adata.obs["doublet_solo"] = adata.obs.index.isin(doublets.index)
        if dd:
            # === DD Method2: DoubletDetection (Boosted classifier on simulated artificial doublets) === #
            clf = doubletdetection.BoostClassifier(
                n_iters=10,
                clustering_algorithm="louvain",
                standard_scaling=True,
                pseudocount=0.1,
                n_jobs=-1,
            )
            doublets = clf.fit(adata.X).predict(p_thresh=1e-3, voter_thresh=0.5)
            doublet_score = clf.doublet_score()
            adata.obs["doublet_dd"] = doublets
            adata.obs["dd_score"] = doublet_score
        if scrublet:
            # === DD Method3 Scrublet (based on simulated doublet neighborhood densities) === #
            sc.pp.scrublet(adata, expected_doublet_rate=0.1)
            adata.obs.rename(
                columns={
                    "predicted_doublet": "doublet_scrublet",
                    "doublet_score": "scrublet_score",
                },
                inplace=True,
            )
        # === Final Doublet Annotation: if all methods agree, then it is a doublet === #
        doublet_cols = adata.obs.columns[adata.obs.columns.str.contains("doublet_")]
        adata.obs["doublet_intersected"] = adata.obs[doublet_cols].sum(axis=1) == len(
            doublet_cols
        )

    # ===================================================================================================
    # Step 2. Hard filters (optional): remove definite low-quality cells
    # ===================================================================================================
    if remove_definite_outliers:
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        # MT gene filter
        adata = adata[adata.obs.pct_counts_mt < 25].copy()

    # ===================================================================================================
    # Step 3. QC metric computation (mt, ribo, optional hb) + removal of MT genes from feature space
    # ===================================================================================================
    adata = add_qc(adata, hb_genes=hb_genes)

    # ===================================================================================================
    # Step 4. Annotate potential outliers (do NOT remove): MAD-based thresholding across metrics
    # ===================================================================================================
    # Sum booleans: any metric flagged => Outlier True
    bool_vector = (
        mad_outlier(adata, "log1p_total_counts", 5)
        + mad_outlier(adata, "log1p_n_genes_by_counts", 5)
        + mad_outlier(adata, "pct_counts_in_top_20_genes", 5)
        + mad_outlier(adata, "pct_counts_mt", 3, upper_only=True)
    )
    adata.obs["Outlier"] = bool_vector  # annotating but not removing !!

    return adata


# ===================================================================================================
# Example Usage
# ===================================================================================================
if __name__ == "__main__":
    ################################################################################################
    ## Example: Load a raw count AnnData, run preprocessing with triple doublet detection methods ##
    ################################################################################################
    # Adjust these paths / parameters for your dataset.
    example_path = "./Data/example_raw_counts.h5ad"

    print(f"[INFO] Loading AnnData from: {example_path}")
    adata_example = sc.read_h5ad(example_path)

    # Ensure raw counts accessible for scVI (store counts in a layer if not present)
    if "counts" not in adata_example.layers:
        adata_example.layers["counts"] = adata_example.X.copy()

    # Run full preprocessing with all doublet methods enabled
    adata_example = pp(
        adata_example,
        doublet_detection=True,
        solo=True,
        dd=True,
        scrublet=True,
        hb_genes=False,
        remove_definite_outliers=True,
    )

    # Summaries
    doublet_cols = [c for c in adata_example.obs.columns if c.startswith("doublet_")]
    print("\n[RESULT] Doublet method columns:", doublet_cols)
    for col in doublet_cols:
        print(
            f"  {col}: {adata_example.obs[col].sum()} / {adata_example.n_obs} flagged"
        )
    if "doublet_intersected" in adata_example.obs:
        print(
            f"  doublet_intersected: {adata_example.obs['doublet_intersected'].sum()} / {adata_example.n_obs} consensus"
        )

    # Optionally write results
    out_annot_path = "Results/doublet_qc_annotations.csv"
    adata_example.obs.to_csv(out_annot_path)
    print(f"[INFO] Saved annotations to {out_annot_path}")
