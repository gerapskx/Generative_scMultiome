#!/usr/bin/env python3
"""UMAP panels coloured by the expression of individual genes.

Reads the raw gene x cell counts written by the ingest step, normalises to
log1p(CP10k), and overlays the n24 UMAP coordinates. Non-expressing cells stay
as a light grey background so the expressing population is easy to see.

    python src/plot_gene_umap.py --step n24 --genes XIST UTY TMEM156 TSIX

Writes one PNG per gene plus a combined grid to the Windows Desktop.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"
DESKTOP = Path("/mnt/c/Users/gerar/OneDrive/Desktop")


def extract_counts(mtx: Path, features: Path, genes: list[str]) -> pd.DataFrame:
    """Pull only the requested gene rows out of the full Matrix Market file."""
    names = pd.read_csv(features, header=None, names=["gene"])["gene"]
    index = {name: i + 1 for i, name in enumerate(names)}  # 1-based rows
    missing = [g for g in genes if g not in index]
    if missing:
        raise SystemExit(f"genes not in {features}: {', '.join(missing)}")
    wanted = {index[g]: g for g in genes}

    cache = mtx.parent.parent / "features" / "umap_gene_counts.tsv"
    if cache.exists():
        cached = pd.read_csv(cache, sep="\t")
        if set(genes) <= set(cached["gene"]):
            print(f"reusing {cache}")
            return cached[cached["gene"].isin(genes)].copy()

    opener = gzip.open if mtx.suffix == ".gz" else open
    rows = []
    with opener(mtx, "rt") as fh:
        next(fh)
        next(fh)
        for line in fh:
            gene, cell, value = line.split()
            name = wanted.get(int(gene))
            if name is not None:
                rows.append((name, int(cell), int(value)))
    frame = pd.DataFrame(rows, columns=["gene", "cell", "count"])

    try:
        frame.to_csv(cache, sep="\t", index=False)
        print(f"cached {cache}")
    except OSError:
        pass
    return frame


def plot_gene(ax, coords, expr, gene, point_size=1.6, cmap="viridis"):
    background = coords.loc[~coords.index.isin(expr.index)]
    ax.scatter(background["umap_1"], background["umap_2"], s=point_size * 0.6,
               color="#E8E8E8", linewidths=0, rasterized=True, zorder=1)
    order = expr["value"].sort_values().index
    sca = ax.scatter(expr.loc[order, "umap_1"], expr.loc[order, "umap_2"],
                     s=point_size, c=expr.loc[order, "value"], cmap=cmap,
                     vmin=0, vmax=float(np.percentile(expr["value"], 99)),
                     linewidths=0, rasterized=True, zorder=2)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1", fontsize=8)
    ax.set_ylabel("UMAP2", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_aspect("equal", adjustable="box")
    ax.margins(0.02)
    for spine in ax.spines.values():
        spine.set_color("#BBBBBB")

    frac = 100 * len(expr) / len(coords)
    ax.set_title(f"{gene}  ({frac:.1f}% of cells)", fontsize=11,
                 fontweight="bold")
    cbar = plt.colorbar(sca, ax=ax, shrink=0.62, pad=0.02, aspect=18)
    cbar.set_label("log1p(CP10k)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    return sca


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", default="n24")
    ap.add_argument("--genes", nargs="+",
                    default=["XIST", "UTY", "TMEM156", "TSIX"])
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--outdir", default=str(DESKTOP))
    ap.add_argument("--point-size", type=float, default=1.6)
    args = ap.parse_args()

    ingest = SCALE / args.step / "ingest"
    mtx = ingest / "rna_counts.mtx.gz"
    if not mtx.exists():
        mtx = ingest / "rna_counts.mtx"
    obs = pd.read_csv(ingest / "obs.csv",
                      usecols=["cell", "nCount_RNA", "group", "donor"])
    coords = pd.read_csv(SCALE / args.step / "cluster" / "obs_annotated.csv",
                         usecols=["cell", "umap_1", "umap_2", "celltype"])
    frame = obs.merge(coords, on="cell", how="inner")
    frame = frame.dropna(subset=["umap_1", "umap_2"]).set_index("cell")

    counts = extract_counts(mtx, ingest / "rna_features.tsv", args.genes)
    # The Matrix Market columns are 1-based positions in the same order as
    # obs.csv (verified: column sums equal nCount_RNA), not cell names.
    counts["cell"] = obs["cell"].to_numpy()[counts["cell"].to_numpy() - 1]
    counts["log1p_cp10k"] = np.log1p(
        1e4 * counts["count"] / frame.loc[counts["cell"], "nCount_RNA"].to_numpy()
    )

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    summary = []
    panels = []
    for gene in args.genes:
        sub = counts[counts["gene"] == gene].set_index("cell")
        expr = pd.DataFrame({
            "umap_1": frame.loc[sub.index, "umap_1"],
            "umap_2": frame.loc[sub.index, "umap_2"],
            "value": sub["log1p_cp10k"],
        })
        panels.append((gene, expr))

        for group, part in frame.groupby("group"):
            in_group = part.index.intersection(sub.index)
            summary.append({
                "gene": gene, "group": group, "cells": len(part),
                "pct_expressing": round(100 * len(in_group) / len(part), 2),
                "mean_log1p_cp10k": round(
                    float(sub.loc[in_group, "log1p_cp10k"].mean()), 3),
            })

    ncol = 2 if len(panels) > 1 else 1
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.4 * ncol, 4.8 * nrow),
                             squeeze=False)
    for ax, (gene, expr) in zip(axes.ravel(), panels):
        plot_gene(ax, frame, expr, gene, point_size=args.point_size)
    for ax in axes.ravel()[len(panels):]:
        ax.axis("off")
    fig.tight_layout()

    grid = outdir / f"GSE193688_UMAP_genes_{args.step}.png"
    if grid.exists():
        grid.unlink()
    fig.savefig(grid, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {grid}")

    for gene, expr in panels:
        fig, ax = plt.subplots(figsize=(5.4, 4.8))
        plot_gene(ax, frame, expr, gene)
        fig.tight_layout()
        out = outdir / f"GSE193688_UMAP_{gene}_{args.step}.png"
        if out.exists():
            out.unlink()
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)
        print(f"wrote {out}")

    table = pd.DataFrame(summary)
    table.to_csv(outdir / f"gene_expression_summary_{args.step}.csv", index=False)
    print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
