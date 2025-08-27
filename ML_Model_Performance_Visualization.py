###############################################################################################
##        Performance of across different models and metrics                                 ##
##                                                                                           ##
##  This script generates heatmap of model performance across different models               ##
##  For each model (and target), get performance metrics                                     ##
##  (precision, recall, specificity, etc.)                                                   ##
##                                                                                           ##
###############################################################################################

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import (
    confusion_matrix,
    classification_report,
    f1_score,
    precision_score,
    accuracy_score,
    recall_score,
)

datadir = "Data"
outdir = "Results"
os.makedirs(outdir, exist_ok=True)

pred_df = pd.read_csv(
    os.path.join(
        datadir,
        "CLL_RT_GenotypingProfiles_clean.csv",
    ),
)


# * ====== Performance metrics ====== * #
def evaluate_prediction(
    pred_df, pred_col="pred", y_col="Y", labels=[0, 1], pos_label=1, report=True
):
    # Validate the labels
    if (
        not pd.Series(pred_df[y_col].unique()).isin(labels).all()
        or not pd.Series(pred_df[pred_col].unique()).isin(labels).all()
        # not pd.Series(labels).isin(pred_df[y_col].unique()).all()
        # or not pd.Series(labels).isin(pred_df[pred_col].unique()).all()
    ):
        print("Labels mismatch between the prediction and the target.")
        print(f"\tLabels Specified: {labels}")
        print(f"\tPrediction labels: {pred_df[pred_col].unique()}")
        print(f"\tTarget labels: {pred_df[y_col].unique()}")
        # acc = prec = recall = f1 = 0
        # return acc, prec, recall, f1, None

    # === Report === #
    # Confusion matrix
    conf_matrix = confusion_matrix(pred_df[y_col], pred_df[pred_col], labels=labels)

    # Classification report
    if report:
        print(classification_report(pred_df[y_col], pred_df[pred_col]))

    # === Metrics scores === #
    # TN, FP, FN, TP
    tn, fp, fn, tp = confusion_matrix(
        pred_df[y_col], pred_df[pred_col], labels=labels
    ).ravel()
    if report:
        print(f"TN: {tn}, FP: {fp}, FN: {fn}, TP: {tp}")

    # Accuracy
    acc_man = (tp + tn) / (tp + tn + fp + fn) if (tp + tn + fp + fn) != 0 else np.nan
    acc = accuracy_score(pred_df[y_col], pred_df[pred_col])
    acc = acc if not np.isnan(acc) else acc_man
    acc_man = acc_man if not np.isnan(acc_man) else acc

    # Precision
    prec_man = tp / (tp + fp) if (tp + fp) != 0 else np.nan
    prec = precision_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    prec = prec if not np.isnan(prec) else prec_man
    prec_man = prec_man if not np.isnan(prec_man) else prec

    # Recall
    recall_man = tp / (tp + fn) if (tp + fn) != 0 else np.nan
    recall = recall_score(
        pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label
    )
    recall = recall if not np.isnan(recall) else recall_man
    recall_man = recall_man if not np.isnan(recall_man) else recall

    # F1 score
    f1_man = (
        2 * (prec_man * recall_man) / (prec_man + recall_man)
        if (prec_man + recall_man) != 0
        and not (np.isnan(prec_man) or np.isnan(recall_man))
        else np.nan
    )
    f1 = f1_score(pred_df[y_col], pred_df[pred_col], labels=labels, pos_label=pos_label)
    f1 = f1 if not np.isnan(f1) else f1_man
    f1_man = f1_man if not np.isnan(f1_man) else f1

    if report:
        print(f"F1: {f1}, F1_man: {f1_man}")
        print(f"Precision: {prec}, Recall: {recall}")

    return acc, prec, recall, f1, conf_matrix


def custom_performance_metric(y_true, y_pred, y_value_of_interest, alpha=0.9):
    """
    Custom performance metric that combines overall F1 score and recall - false positive rate (FPR) for a specific class.

    Parameters:
    - y_true (array-like): True labels.
    - y_pred (array-like): Predicted labels.
    - y_value_of_interest (int or str): The class of interest.

    Returns:
    - float: The composite performance metric value.
    """

    # * Weighted combination of precision and recall (more emphasis on precision) of the y_value_of_interest
    report = classification_report(y_true, y_pred, output_dict=True, zero_division=0)
    precision = report[str(y_value_of_interest)]["precision"]
    recall = report[str(y_value_of_interest)]["recall"]
    performance = alpha * precision + (1 - alpha) * recall

    return performance


# TODO: For each target and each model (wanna get output df with each row specific to a model and a target)
performance_df = []
targets = pred_df["target"].unique()
# targets = ["PLCG2_Mut2", "BTK", "JUNB"]  # ! For testing
for target in targets:
    # print(f"\n# === {target} === #")
    target_df = pred_df.query(f"target == '{target}'")

    # * === Gather all predictions and probabilities === #
    target_df.rename(
        columns={"ensemble_hard_confidence": "ensemble_hard_prob_MUT"}, inplace=True
    )
    pred_cols = target_df.columns[target_df.columns.str.contains("_pred$")].tolist()
    prob_cols = target_df.columns[target_df.columns.str.contains("_prob_")].tolist()

    # TODO: For each of model, adjust threshold of prediction probability
    # TODO: for the given threshold of MUT prediction probability, calculate performance
    models = [c.replace("_pred", "") for c in pred_cols]
    for model_name in models:
        model = model_name

        pred_col = f"{model}_pred"
        prob_col = f"{model}_prob_MUT"

        # ! Replace FalsePositive with WT
        target_df[pred_col] = target_df[pred_col].str.replace("FalsePositive", "WT")

        # Calculate false positive rate
        false_positive_rate = (
            target_df.query("expected == 0")[pred_col]
            .value_counts(normalize=True)
            .get("MUT")
        )
        # Calculate true positive rate
        true_positive_rate = (
            target_df.query("expected == 1")[pred_col]
            .value_counts(normalize=True)
            .get("MUT")
        )

        acc, prec, recall, f1, _conf_matrix = evaluate_prediction(
            target_df,
            pred_col=pred_col,
            y_col="gt",
            labels=["WT", "MUT"],
            pos_label="MUT",
            report=False,
        )

        combined_performance = custom_performance_metric(
            target_df["gt"],
            target_df[pred_col],
            y_value_of_interest="MUT",
        )

        new_row = {
            "target": target,
            "expected_sample": target_df["expected_sample"].unique().item(),
            "model": model_name,
            # "acc": acc,
            "prec": prec,
            "recall": recall,
            "f1": f1,
            "combined_performance": combined_performance,
            "mut_freq_expected_sample": true_positive_rate,
            "mut_freq_non_expected_sample": false_positive_rate,
        }
        performance_df.append(new_row)
performance_df = pd.DataFrame(performance_df)


y = "target"
x = "model"
df = performance_df.copy()

col_order = [
    "logistic_regression",
    "random_forest",
    "knn",
    "naive_bayes",
    "xgboost",
    # "gradient_boosting",
    "hist_gradient_boosting",
    # "adaboost",
    # "mlp",
    "ensemble_hard",
    "ensemble_soft",
]
row_order = [
    "B2M",
    "PLCG2_Mut1",
    "PLCG2_Mut3",
    "PLCG2_Mut2",
    "EMD",
    "IRF8_Mut1",
    "MCL1",
    "POU2F2_Mut1",
    "SRRM2",
    "BTK",
    "JUNB",
    "SF3B1",
    "PRKDC",
    "IKZF3",
    "JAK1",
]

fontsize = 14
values = ["prec", "recall", "combined_performance"]
for value in values:
    df_pivot = df.pivot(
        index=y,
        columns=x,
        values=value,
    )
    df_pivot = df_pivot.reindex(row_order)
    df_pivot = df_pivot[col_order]

    # Create the clustermap with the reordered data
    ax = sns.clustermap(
        df_pivot,  # Use pre-reordered data
        cmap="Blues",
        cbar_pos=(1.21, 0.45, 0.03, 0.35),
        figsize=(0.55 * df_pivot.shape[1], 0.55 * df_pivot.shape[0]),
        row_cluster=False,
        col_cluster=False,
        dendrogram_ratio=(0.0, 0.00001),
        col_linkage=None,
        annot=False,  # data_raw_reordered.values,
    )

    # Set the Colorbar label font size
    cbar = ax.ax_cbar
    cbar.set_ylabel(
        cbar.get_ylabel(), fontsize=fontsize + 2, fontweight="normal", labelpad=20
    )  # Set font size and weight
    cbar.set_yticklabels(
        cbar.get_yticklabels(), fontsize=fontsize, weight="normal"
    )  # Set the colorbar tick labels font sizes

    # Remove x and y labels
    ax.ax_heatmap.set_xlabel("", fontsize=fontsize + 2, weight="bold", labelpad=10)
    ax.ax_heatmap.set_ylabel("", fontsize=fontsize + 2, weight="bold", labelpad=10)

    # Format x-axis tick labels to show only 2 decimal places if numeric
    # Use the actual column names from df_pivot instead of tick positions
    new_x_ticks = []

    for col_name in df_pivot.columns:
        try:
            # If the column name is numeric, format it to 2 decimal places
            new_text = f"{float(col_name):.1f}"
            new_x_ticks.append(new_text)
        except ValueError:
            # If it's not numeric, keep the original text
            new_x_ticks.append(col_name)

    # Set the axes tick labels font sizes
    ytick_texts = ax.ax_heatmap.set_yticklabels(
        ax.ax_heatmap.get_yticklabels(),
        fontsize=fontsize,
        weight="normal",
        rotation=0,
        ha="left",
    )
    xtick_texts = ax.ax_heatmap.set_xticklabels(
        new_x_ticks, fontsize=fontsize, weight="normal"
    )

    # Title
    title_text = ax.ax_heatmap.set_title(
        f"{value}",
        fontsize=fontsize + 2,
        weight="bold",
        pad=20,
        color="0.2",
    )

    fig = ax.figure
    fig.savefig(
        os.path.join(
            outdir,
            f"Heatmap_performance_per_model_{value}.pdf",
        ),
        bbox_inches="tight",
        dpi=300,
    )
    plt.show()
