import os
import warnings
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import collections
from copy import deepcopy
from scipy import io
from scipy.stats import median_abs_deviation as mad

# ====== Ingest (Integrate / Label Transfer) ====== #
def ingest_label(adata, adata_ref, transfer_column, batch_key=None):
    """
    Ingests batches of an AnnData object into a reference AnnData object using a specified transfer column and batch key.

    Parameters:
    adata (AnnData): AnnData containing the batches to be ingested (normalized : run cluster_adata() first).
    adata_ref (AnnData): The reference AnnData (normalized + PCA + neighbor + UMAP : run cluster_adata() first).
    transfer_column (str): The column to transfer the labels from the reference
    batch_key (str): The column indicating different batches.

    Returns:
    AnnData: The ingested AnnData object with all batches integrated into the reference.
    """
    if batch_key is None:
        batch_key = "_batch"
        adata.obs[batch_key] = "query"
    adata_ref.obs[batch_key] = "reference"

    # save the original barcode to avoid unwanted suffixes
    adata.obs["_barcode"] = adata.obs.index
    adata_ref.obs["_barcode"] = adata_ref.obs.index

    var_names = adata_ref.var_names.intersection(adata.var_names)
    adata_ref = adata_ref[:, var_names]
    adata = adata[:, var_names]

    adatas = []
    batch_list = list(adata.obs[batch_key].unique())
    for batch in batch_list:
        ad = adata[adata.obs[batch_key] == batch].copy()
        print(f"... integrating batch {batch} into reference")
        ad.obs[f"{transfer_column}_orig"] = ad.obs[
            transfer_column
        ].copy()  # save the original cell type
        sc.tl.ingest(ad, adata_ref, obs=transfer_column)
        adatas.append(ad)

    adata_ingested = adata_ref.concatenate(
        adatas, batch_categories=["Ref"] + batch_list
    )

    # === Post-processing === #
    adata_ingested.obs[transfer_column] = adata_ingested.obs[transfer_column].astype(
        "category"
    )
    # fix category ordering
    adata_ingested.obs[transfer_column] = adata_ingested.obs[
        transfer_column
    ].cat.reorder_categories(adata_ref.obs[transfer_column].cat.categories)
    # fix category coloring
    adata_ingested.uns[f"{transfer_column}_colors"] = adata_ref.uns[
        f"{transfer_column}_colors"
    ]

    if batch_key == "_batch":
        adata_ingested.obs.drop(columns=[batch_key], inplace=True)

    # reset the original barcode
    adata_ingested.obs.set_index("_barcode", inplace=True)
    # adata_ingested.obs.drop(columns=['_barcode'], inplace=True)

    return adata_ingested