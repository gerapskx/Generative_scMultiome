#!/usr/bin/env python3
"""Stage 1 - assemble the exported MTX files into AnnData.

Reads results/ingest/ and writes two h5ad files with identical obs and cell order.
RNA and ATAC have different feature spaces, so they are separate AnnData objects
rather than layers of one.

Usage: python src/assemble.py [--ingest results/ingest] [--out data/processed]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def find(path_stem: Path) -> Path:
    """Prefer the gzipped variant if it exists."""
    gz = Path(str(path_stem) + ".gz")
    if gz.exists():
        return gz
    if path_stem.exists():
        return path_stem
    raise FileNotFoundError(f"neither {path_stem} nor {gz} exists")


def load_modality(ingest: Path, name: str, obs: pd.DataFrame) -> ad.AnnData:
    mtx = find(ingest / f"{name}_counts.mtx")
    features = (ingest / f"{name}_features.tsv").read_text().splitlines()

    mat = scipy.io.mmread(mtx)          # features x cells
    if sp.issparse(mat):
        mat = mat.tocsr()
    else:
        mat = sp.csr_matrix(mat)

    if mat.shape[0] != len(features):
        raise ValueError(f"{name}: {mat.shape[0]} rows but {len(features)} features")
    if mat.shape[1] != len(obs):
        raise ValueError(f"{name}: {mat.shape[1]} columns but {len(obs)} cells")

    if mat.dtype.kind == "f" and not np.allclose(mat.data, np.round(mat.data)):
        raise ValueError(f"{name}: counts are not integers")

    adata = ad.AnnData(
        X=mat.T.tocsr().astype(np.int32),   # cells x features
        obs=obs.copy(),
        var=pd.DataFrame(index=pd.Index(features, name="feature")),
    )
    adata.obs_names = obs.index
    return adata


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--out", default="data/processed")
    args = ap.parse_args()

    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(ingest / "obs.csv", index_col="cell", dtype={"cell": str})
    if obs.index.has_duplicates:
        raise ValueError("duplicated cell names in obs.csv")
    for col in ("sample", "donor", "group"):
        if col not in obs.columns:
            raise ValueError(f"obs.csv missing column {col!r}")

    rna = load_modality(ingest, "rna", obs)
    atac = load_modality(ingest, "atac", obs)

    if not (rna.obs_names.equals(atac.obs_names)):
        raise ValueError("RNA and ATAC cell order diverged during assembly")

    # canonical group order: Aged is the reference everywhere downstream
    obs["group"] = pd.Categorical(obs["group"], categories=["Aged", "PD", "Young"])
    rna.obs["group"] = obs["group"].values
    atac.obs["group"] = obs["group"].values

    rna.write_h5ad(out / "rna.h5ad", compression="gzip")
    atac.write_h5ad(out / "atac.h5ad", compression="gzip")

    print(f"rna  : {rna.n_obs} cells x {rna.n_vars} genes  -> {out / 'rna.h5ad'}")
    print(f"atac : {atac.n_obs} cells x {atac.n_vars} peaks -> {out / 'atac.h5ad'}")
    print(f"groups: {obs['group'].value_counts().to_dict()}")
    print(f"donors: {obs['donor'].nunique()}  samples: {obs['sample'].nunique()}")

    if obs["donor"].nunique() < 2 and obs["group"].nunique() > 1:
        print(
            "WARNING: one donor per group - condition is confounded with donor",
            file=sys.stderr,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
