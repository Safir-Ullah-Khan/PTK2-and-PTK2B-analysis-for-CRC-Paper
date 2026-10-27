#!/usr/bin/env python3
import argparse
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

GENE_ORDER = [
    "Epcam", "Cdh1", "Cdh6", "Cdh17",
    "Krt8", "Krt14", "Krt15", "Krt18", "Krt19", "Krt80",
    "Pcdhb19", "Pcdhga3", "Pcdhgb6", "Pcdh11x",
]

GROUP_MAP = {
    "HGS2":  ["HGS2A", "HGS2B", "HGS2C"],
    "KO":    ["M38KOA", "M38KOB", "M38KOC"],
    "WT":    ["M38WTA", "M38WTB", "M38WTC"],
    "K457A": ["M38457A", "M38457B", "M38457C"],
}
COLUMN_ORDER = ["HGS2", "KO", "WT", "K457A"]


def set_style():
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    for key in ("text.color", "axes.labelcolor", "xtick.color", "ytick.color"):
        mpl.rcParams[key] = "black"


def build_tables(fpkm_path):
    fpkm = pd.read_csv(fpkm_path, sep="\t")
    expr = fpkm[fpkm["gene_name"].isin(GENE_ORDER)].set_index("gene_name")
    expr = expr.reindex(GENE_ORDER)
    missing = expr[expr.isna().all(axis=1)].index.tolist()
    if missing:
        print("WARNING - genes not found in FPKM matrix:", missing)
    expr_avg = pd.DataFrame({g: expr[cols].mean(axis=1) for g, cols in GROUP_MAP.items()})
    expr_avg = expr_avg[COLUMN_ORDER]
    log_expr = np.log2(expr_avg.astype(float) + 1)
    v = log_expr.values
    z = (v - v.mean(axis=1, keepdims=True)) / v.std(axis=1, keepdims=True)
    z_df = pd.DataFrame(z, index=log_expr.index, columns=log_expr.columns)
    return expr_avg, z_df


def plot_heatmap(z_df, out_dir):
    n = len(z_df)
    fig, ax = plt.subplots(figsize=(4.2, n * 0.35))
    sns.heatmap(
        z_df, cmap="RdBu_r", center=0, linewidths=0.5, linecolor="white",
        square=True,
        cbar_kws={"label": "z-score\n(log2 FPKM+1)", "shrink": 0.5, "aspect": 15, "pad": 0.03},
        yticklabels=True, xticklabels=True, ax=ax,
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9, fontweight="bold")
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=8, fontweight="bold", style="italic", rotation=0)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(f"PYK2-dependent epithelial genes (n={n})", fontsize=10, fontweight="bold", pad=10)
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("z-score (log2 FPKM+1)", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "Heatmap_14gene_log2_zscore.png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.savefig(os.path.join(out_dir, "Heatmap_14gene_log2_zscore.pdf"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True)
    parser.add_argument("--out", default="outputs")
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    set_style()
    fpkm_path = os.path.join(args.base, "3.Quant", "1.Count", "gene_fpkm.xls")
    expr_avg, z_df = build_tables(fpkm_path)
    expr_avg.round(2).to_excel(os.path.join(args.out, "Heatmap_14gene_averaged_FPKM.xlsx"))
    z_df.round(3).to_excel(os.path.join(args.out, "Heatmap_14gene_zscores.xlsx"))
    plot_heatmap(z_df, args.out)
    print(f"Saved to: {args.out}")


if __name__ == "__main__":
    main()
