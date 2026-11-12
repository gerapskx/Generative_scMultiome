#!/usr/bin/env python3
"""Automatic cluster annotation from a canonical marker panel.

Scaling steps re-cluster at every N, so the manual annotation gate cannot be used.
This script scores each joint cluster on a fixed panel of canonical glial / neuronal
markers, z-scores the scores across clusters, and assigns the argmax cell type.

The score is a relative ranking, not a probability. `margin` (best minus runner-up)
and `best_score` are written so low-confidence clusters are visible.

Usage:
    python src/auto_annotate.py --ingest results/scale/n08/ingest \\
                               --clusters results/scale/n08/cluster/clusters.csv \\
                               --markers results/scale/n08/cluster/cluster_markers.csv \\
                               --out results/scale/n08/cluster

Outputs (out/):
    celltype_map.csv               cluster, celltype, score, runner_up, ...
    cluster_annotation_scores.csv   cluster x celltype z-scores
    obs_annotated.csv              clusters.csv + celltype
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PANEL = {
    "Oligodendrocyte": ["MBP", "PLP1", "MOG", "MAG", "CNP", "CLDN11", "MOBP"],
    "Astrocyte": ["GFAP", "AQP4", "SLC1A2", "SLC1A3", "GJA1", "ALDH1L1"],
    "Microglia": ["CSF1R", "P2RY12", "CX3CR1", "AIF1", "TMEM119", "C1QA", "C1QB"],
    "OPC": ["PDGFRA", "OLIG1", "OLIG2", "CSPG4", "VCAN", "LHFPL3"],
    "Neuron": ["SYT1", "SNAP25", "RBFOX3", "TH", "SLC6A3", "SLC17A7", "GAD1", "GAD2"],
    "Endothelial": ["CLDN5", "FLT1", "PECAM1", "VWF", "KDR"],
    "Pericyte": ["RGS5", "PDGFRB", "ACTA2", "NOTCH3"],
    "Tcell": ["CD3D", "CD3E", "CD3G", "PTPRC", "CD2"],
    "Ependymal": ["FOXJ1", "PIFO", "DNAH9", "CFAP43"],
}


def panel_expression(ingest: Path, panel: dict[str, list[str]]) -> tuple[np.ndarray, list[str], list[str]]:
    """Return log1p CP10K expression for the panel genes only, plus names."""
    genes = [l for l in (ingest / "rna_features.tsv").read_text().splitlines() if l]
    gene_pos = {g: i for i, g in enumerate(genes)}
    wanted = sorted({g for genes_ in panel.values() for g in genes_ if g in gene_pos})
    if not wanted:
        raise SystemExit("none of the panel genes are present in rna_features.tsv")
    row_pos = {gene_pos[g] + 1: j for j, g in enumerate(wanted)}

    obs = pd.read_csv(ingest / "obs.csv", dtype={"cell": str})
    n_cells = len(obs)
    mat = np.zeros((len(wanted), n_cells), dtype=np.float64)
    lib = np.zeros(n_cells, dtype=np.float64)

    mtx = ingest / "rna_counts.mtx.gz"
    if not mtx.exists():
        mtx = ingest / "rna_counts.mtx"
    # MatrixMarket: line 1 header, line 2 dimensions, then row col value.
    for chunk in pd.read_csv(
        mtx, sep=r"\s+", skiprows=2, header=None, names=["row", "col", "val"],
        chunksize=2_000_000, engine="c", dtype={"row": np.int64, "col": np.int64,
                                                 "val": np.float64},
    ):
        cols = chunk["col"].to_numpy() - 1
        vals = chunk["val"].to_numpy()
        lib += np.bincount(cols, weights=vals, minlength=n_cells)
        sel = chunk[chunk["row"].isin(row_pos)]
        if not sel.empty:
            rows = sel["row"].map(row_pos).to_numpy()
            cols = sel["col"].to_numpy() - 1
            np.add.at(mat, (rows, cols), sel["val"].to_numpy())

    lib[lib == 0] = 1.0
    return np.log1p(mat / lib[None, :] * 1e4).T, wanted, obs["cell"].tolist()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--clusters", default="results/cluster/clusters.csv")
    ap.add_argument("--markers", default="results/cluster/cluster_markers.csv")
    ap.add_argument("--out", default="results/cluster")
    ap.add_argument("--min-score", type=float, default=0.5)
    args = ap.parse_args()

    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    clusters = pd.read_csv(args.clusters, dtype={"cell": str, "cluster": str})
    if "cluster" not in clusters.columns:
        raise SystemExit(f"{args.clusters} needs a `cluster` column")

    x, genes, cells = panel_expression(ingest, PANEL)
    if len(cells) != x.shape[0]:
        raise SystemExit("obs row count does not match the count matrix")
    cell_index = {c: i for i, c in enumerate(cells)}
    missing = [c for c in clusters["cell"] if c not in cell_index]
    if missing:
        raise SystemExit(f"{len(missing)} cluster cells are absent from the ingest")

    cl = clusters.set_index("cell").loc[cells, "cluster"].to_numpy()
    cluster_levels = sorted(pd.unique(clusters["cluster"].astype(str)))
    means = np.vstack([x[cl == c].mean(axis=0) for c in cluster_levels])
    global_mean = x.mean(axis=0, keepdims=True)

    # Relative enrichment, not a z-score across clusters: in a cohort dominated
    # by one cell type, z-scoring would penalise the markers that are high in
    # many clusters (the oligodendrocyte genes here). log2(cluster/global) keeps
    # the absolute enrichment and is comparable across cell types.
    eps = 1e-6
    z = np.log2((means + eps) / (global_mean + eps))

    gene_index = {g: i for i, g in enumerate(genes)}
    scores = {}
    for celltype, panel_genes in PANEL.items():
        idx = [gene_index[g] for g in panel_genes if g in gene_index]
        scores[celltype] = z[:, idx].mean(axis=1) if idx else np.full(len(cluster_levels), -np.inf)
    score_df = pd.DataFrame(scores, index=cluster_levels)
    score_df.index.name = "cluster"

    best = score_df.idxmax(axis=1)
    best_score = score_df.max(axis=1)
    runner = score_df.apply(lambda r: r.drop(r.idxmax()).max(), axis=1)
    runner_name = score_df.apply(lambda r: r.drop(r.idxmax()).idxmax(), axis=1)

    cellmap = pd.DataFrame({
        "cluster": cluster_levels,
        "celltype": best.to_numpy(),
        "score": best_score.to_numpy(),
        "runner_up": runner_name.to_numpy(),
        "runner_up_score": runner.to_numpy(),
        "margin": (best_score - runner).to_numpy(),
        "low_confidence": (best_score < args.min_score).to_numpy(),
    })
    cellmap.to_csv(out / "celltype_map.csv", index=False)
    score_df.reset_index().to_csv(out / "cluster_annotation_scores.csv", index=False)

    lookup = dict(zip(cellmap["cluster"], cellmap["celltype"]))
    annotated = clusters.copy()
    annotated["celltype"] = annotated["cluster"].map(lookup)
    annotated.to_csv(out / "obs_annotated.csv", index=False)

    counts = annotated["celltype"].value_counts()
    summary = {
        "panel_genes_present": genes,
        "n_clusters": len(cluster_levels),
        "n_cells": len(annotated),
        "min_score": args.min_score,
        "low_confidence_clusters": cellmap.loc[cellmap["low_confidence"], "cluster"].tolist(),
        "celltypes": counts.to_dict(),
    }
    (out / "auto_annotation_summary.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote {out}/celltype_map.csv, cluster_annotation_scores.csv, "
          f"obs_annotated.csv ({len(annotated)} cells)")
    print(counts.to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
