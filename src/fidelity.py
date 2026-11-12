#!/usr/bin/env python3
"""Stage 7 - ML fidelity of the scDesign3 simulations.

Two independent questions, both answered per modality (RNA, ATAC):

1. Discriminator: can a classifier tell real holdout cells from synthetic cells?
   A faithful simulator gives AUC ~0.5. Controls bracket the scale:
     identity        real vs real (two halves)     -> AUC ~0.5 by construction
     positive        Aged vs PD real cells          -> AUC ~1.0 (the test has power)
     naive_poisson   real vs Poisson(marginal means) -> how much structure is lost
                                                      by a marginal-only model
   Report `auc` and `auc_dist = |auc - 0.5|`: a bad simulator can be separable in
   either direction, so distance from chance is the quantity to compare.

2. Label transfer: train on one domain, predict `celltype` on the other.
   Synthetic -> real is the direction that matters: a simulator that preserves
   biology should produce cells a real-trained problem can still classify, and a
   synthetic-trained classifier should still recognise real cells.

Every comparison is repeated per simulation seed, so the reported mean/sd is
across seeds rather than across folds only.

Usage:
    python src/fidelity.py [--ingest results/ingest] [--features results/features]
                           [--split results/split/split.csv]
                           [--obs results/cluster/obs_annotated.csv]
                           [--sim results/simulation] [--mode alt]
                           [--out results/fidelity] [--seeds 101,202]

Outputs (results/fidelity/):
    discriminator.csv      modality, comparison, model, auc_mean, auc_sd, ...
    label_transfer.csv     modality, direction, model, accuracy, balanced_accuracy, ...
    per_class.csv          per-celltype recall for synthetic -> real
    summary.json           parameters + headline numbers
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RNG = np.random.default_rng(20261110)


# --------------------------------------------------------------------- input --
def find(path_stem: Path) -> Path:
    """Prefer the gzipped variant if it exists."""
    gz = Path(str(path_stem) + ".gz")
    if gz.exists():
        return gz
    if path_stem.exists():
        return path_stem
    raise FileNotFoundError(f"neither {path_stem} nor {gz} exists")


def read_mtx(path: Path) -> sp.csr_matrix:
    with path.open("rb") as fh:
        if path.suffix == ".gz":
            import gzip

            with gzip.open(fh, "rb") as gz:
                return sp.csr_matrix(scipy.io.mmread(gz))
        return sp.csr_matrix(scipy.io.mmread(fh))


def load_features(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


def subset_counts(
    ingest: Path, modality: str, features: list[str], cells: list[str], all_cells: list[str]
) -> sp.csr_matrix:
    """Return features x cells for the requested feature/cell subset."""
    mat = read_mtx(find(ingest / f"{modality}_counts.mtx"))  # features x cells
    feat_all = load_features(ingest / f"{modality}_features.tsv")
    if mat.shape[0] != len(feat_all):
        raise SystemExit(f"{modality}: matrix rows != feature list")

    feat_pos = {name: i for i, name in enumerate(feat_all)}
    missing = [f for f in features if f not in feat_pos]
    if missing:
        raise SystemExit(f"{modality}: {len(missing)} selected features absent, e.g. {missing[:3]}")
    rows = np.array([feat_pos[f] for f in features], dtype=np.int64)

    cell_pos = {name: i for i, name in enumerate(all_cells)}
    cols = np.array([cell_pos[c] for c in cells], dtype=np.int64)
    return mat[rows][:, cols].tocsr()


# ---------------------------------------------------------------- transform --
def normalise(mat: sp.csr_matrix, modality: str) -> np.ndarray:
    """Per-cell normalisation, shared across real and synthetic cells.

    RNA: log1p CPM.  ATAC: TF-IDF, the same transform used for feature selection.
    """
    mat = mat.astype(np.float64).T.tocsr()  # cells x features
    if modality == "rna":
        lib = np.asarray(mat.sum(axis=1)).ravel()
        lib[lib == 0] = 1.0
        out = mat.multiply(1e6 / lib[:, None]).tocsr()
        return np.log1p(out.toarray())
    # ATAC TF-IDF
    lib = np.asarray(mat.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    tf = mat.multiply(1.0 / lib[:, None]).tocsr()
    df = np.asarray((mat > 0).sum(axis=0)).ravel()
    idf = np.log1p(mat.shape[0] / np.maximum(df, 1))
    return (tf.multiply(idf[None, :])).toarray()


def prepare(real: np.ndarray, syn: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Scale jointly so real and synthetic share one feature space."""
    both = np.vstack([real, syn])
    scaler = StandardScaler().fit(both)
    return scaler.transform(real), scaler.transform(syn)


# ------------------------------------------------------------------- models --
def classifier(name: str):
    if name == "logistic_regression":
        return make_pipeline(
            LogisticRegression(max_iter=3000, C=1.0, class_weight="balanced")
        )
    if name == "random_forest":
        return RandomForestClassifier(
            n_estimators=300, min_samples_leaf=2, n_jobs=-1, random_state=0,
            class_weight="balanced_subsample",
        )
    raise SystemExit(f"unknown model: {name}")


def cv_auc(x: np.ndarray, y: np.ndarray, model: str, folds: int) -> float:
    y = np.asarray(y).astype(int)
    n_min = int(np.bincount(y).min())
    folds = max(2, min(folds, n_min))
    cv = StratifiedKFold(n_splits=folds, shuffle=True, random_state=0)
    proba = cross_val_predict(classifier(model), x, y, cv=cv, method="predict_proba")
    return float(roc_auc_score(y, proba[:, 1]))


def transfer(
    x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, y_test: np.ndarray,
    model: str,
) -> dict:
    clf = classifier(model).fit(x_train, y_train)
    pred = clf.predict(x_test)
    labels = sorted(set(y_test))
    cm = confusion_matrix(y_test, pred, labels=labels)
    recall = {
        lab: float(cm[i, i] / cm[i].sum()) if cm[i].sum() else float("nan")
        for i, lab in enumerate(labels)
    }
    return {
        "accuracy": float((pred == y_test).mean()),
        "balanced_accuracy": float(balanced_accuracy_score(y_test, pred)),
        "macro_f1": float(f1_score(y_test, pred, average="macro")),
        "per_class_recall": recall,
    }


# -------------------------------------------------------------------- naive ---
def naive_poisson(real_counts: sp.csr_matrix, groups: np.ndarray, modality: str) -> np.ndarray:
    """Poisson counts with the per-group marginal mean, library size preserved.

    This is the "simulator with no biology" control: it keeps each gene's mean but
    destroys co-expression and overdispersion, so a good discriminator should
    reject it.
    """
    counts = real_counts.T.tocsr()  # cells x features
    out = np.zeros(counts.shape, dtype=np.float64)
    for g in np.unique(groups):
        sel = np.flatnonzero(groups == g)
        block = counts[sel]
        lib = np.asarray(block.sum(axis=1)).ravel()
        mean_profile = np.asarray(block.sum(axis=0)).ravel()
        mean_profile = mean_profile / max(mean_profile.sum(), 1.0)
        for k, i in enumerate(sel):
            mu = mean_profile * max(lib[k], 1.0)
            out[i] = RNG.poisson(mu)
    if modality == "rna":
        lib = out.sum(axis=1)
        lib[lib == 0] = 1.0
        return np.log1p(out * (1e6 / lib[:, None]))
    lib = out.sum(axis=1)
    lib[lib == 0] = 1.0
    tf = out / lib[:, None]
    df = (out > 0).sum(axis=0)
    idf = np.log1p(out.shape[0] / np.maximum(df, 1))
    return tf * idf[None, :]


# --------------------------------------------------------------------- main --
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--split", default="results/split/split.csv")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--sim", default="results/simulation")
    ap.add_argument("--mode", default="alt")
    ap.add_argument("--out", default="results/fidelity")
    ap.add_argument("--seeds", default="", help="comma list; default = every seed in manifest")
    ap.add_argument("--models", default="logistic_regression,random_forest")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-cells", type=int, default=0,
                    help="subsample each domain for quick runs (0 = all)")
    args = ap.parse_args()

    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(args.obs, dtype={"cell": str})
    split = pd.read_csv(args.split, dtype={"cell": str})
    manifest = pd.read_csv(Path(args.sim) / "manifest.csv", keep_default_na=False)
    manifest["seed"] = manifest["seed"].astype(int)

    holdout = split.loc[split["split"] == "holdout", "cell"].tolist()
    fit = split.loc[split["split"] == "fit", "cell"].tolist()
    meta = obs.set_index("cell")

    # every cell of the cohort, in the order the matrices were written
    all_cells = obs["cell"].tolist()

    want = [s for s in args.seeds.split(",") if s.strip()]
    if want:
        manifest = manifest[manifest["seed"].astype(str).isin(want)]
    manifest = manifest[manifest["mode"] == args.mode]
    if manifest.empty:
        raise SystemExit(f"no simulations for mode={args.mode}")
    seeds = sorted(manifest["seed"].unique())

    models = [m for m in args.models.split(",") if m]
    disc_rows: list[dict] = []
    transfer_rows: list[dict] = []
    per_class_rows: list[dict] = []

    for modality in ("rna", "atac"):
        feats = load_features(Path(args.features) / f"{modality}_features.tsv")
        if len(feats) < 10:
            raise SystemExit(f"{modality}: too few selected features")
        print(f"[{modality}] {len(feats)} features", flush=True)

        real_counts = subset_counts(ingest, modality, feats, holdout, all_cells)
        real_x = normalise(real_counts, modality)
        y_ct = meta.loc[holdout, "celltype"].astype(str).to_numpy()
        y_cond = meta.loc[holdout, "group"].astype(str).to_numpy()

        if args.max_cells and len(holdout) > args.max_cells:
            keep = np.sort(RNG.choice(len(holdout), args.max_cells, replace=False))
            real_x, y_ct, y_cond = real_x[keep], y_ct[keep], y_cond[keep]

        # ---- identity control: real vs real -------------------------------
        half = len(real_x) // 2
        perm = RNG.permutation(len(real_x))
        x_id = np.vstack([real_x[perm[:half]], real_x[perm[half:]]])
        y_id = np.r_[np.zeros(half), np.ones(len(real_x) - half)]
        x_id = StandardScaler().fit_transform(x_id)
        for model in models:
            disc_rows.append({
                "modality": modality, "comparison": "identity", "model": model,
                "seed": np.nan, "auc": cv_auc(x_id, y_id, model, args.folds),
                "n_real": half, "n_syn": len(real_x) - half,
            })

        # ---- positive control: two real, biologically distinct groups -----
        # Within-cell feature shuffling is too weak on sparse counts, so use a
        # contrast the pipeline is actually meant to detect: condition.  If the
        # discriminator cannot separate Aged from PD cells, it has no power to
        # judge the simulation.  (In the 2-donor pilot this is also a donor
        # contrast; that is fine, it only calibrates the test.)
        conds = np.unique(y_cond)
        if len(conds) >= 2:
            a, b = conds[0], conds[1]
            xa, xb = real_x[y_cond == a], real_x[y_cond == b]
            n_pos = min(len(xa), len(xb))
            if n_pos >= 10:
                x_pos = np.vstack([xa[:n_pos], xb[:n_pos]])
                y_pos = np.r_[np.zeros(n_pos), np.ones(n_pos)]
                x_pos = StandardScaler().fit_transform(x_pos)
                for model in models:
                    disc_rows.append({
                        "modality": modality, "comparison": "positive_control",
                        "model": model, "seed": np.nan,
                        "auc": cv_auc(x_pos, y_pos, model, args.folds),
                        "n_real": n_pos, "n_syn": n_pos,
                        "detail": f"{a}_vs_{b}",
                    })

        # ---- naive Poisson control ---------------------------------------
        x_nv = np.vstack([real_x, naive_poisson(real_counts, y_cond, modality)])
        y_nv = np.r_[np.zeros(len(real_x)), np.ones(len(real_x))]
        x_nv = StandardScaler().fit_transform(x_nv)
        for model in models:
            disc_rows.append({
                "modality": modality, "comparison": "naive_poisson", "model": model,
                "seed": np.nan, "auc": cv_auc(x_nv, y_nv, model, args.folds),
                "n_real": len(real_x), "n_syn": len(real_x),
            })

        # ---- per-seed real vs synthetic ----------------------------------
        for seed in seeds:
            row = manifest[(manifest["modality"] == modality) &
                           (manifest["mode"] == args.mode) &
                           (manifest["seed"] == seed)]
            if row.empty:
                continue
            sim_path = Path(row.iloc[0]["path"])
            syn = read_mtx(sim_path).tocsr()
            if syn.shape[0] != len(feats):
                raise SystemExit(f"{sim_path}: {syn.shape[0]} rows, expected {len(feats)}")
            sim_obs = pd.read_csv(sim_path.parent / "obs.csv", dtype={"cell": str})
            if syn.shape[1] != len(sim_obs):
                raise SystemExit(f"{sim_path}: column count != simulated obs")
            syn_x = normalise(syn, modality)

            n = min(len(real_x), len(syn_x))
            if args.max_cells:
                n = min(n, args.max_cells)
            r_idx = np.sort(RNG.choice(len(real_x), n, replace=False))
            s_idx = np.sort(RNG.choice(len(syn_x), n, replace=False))
            xr, xs = prepare(real_x[r_idx], syn_x[s_idx])

            x = np.vstack([xr, xs])
            y = np.r_[np.zeros(n), np.ones(n)]
            for model in models:
                disc_rows.append({
                    "modality": modality, "comparison": "real_vs_synthetic",
                    "model": model, "seed": int(seed),
                    "auc": cv_auc(x, y, model, args.folds), "n_real": n, "n_syn": n,
                })

            # ---- label transfer ------------------------------------------
            y_syn = sim_obs["celltype"].astype(str).to_numpy()[s_idx]
            y_real = y_ct[r_idx]
            shared = sorted(set(y_syn) & set(y_real))
            keep_syn = np.isin(y_syn, shared)
            keep_real = np.isin(y_real, shared)
            for model in models:
                if keep_syn.sum() >= 10 and keep_real.sum() >= 10:
                    res = transfer(xs[keep_syn], y_syn[keep_syn],
                                   xr[keep_real], y_real[keep_real], model)
                    transfer_rows.append({
                        "modality": modality, "direction": "synthetic_to_real",
                        "model": model, "seed": int(seed), **{k: res[k] for k in
                        ("accuracy", "balanced_accuracy", "macro_f1")},
                        "n_train": int(keep_syn.sum()), "n_test": int(keep_real.sum()),
                    })
                    for lab, rec in res["per_class_recall"].items():
                        per_class_rows.append({
                            "modality": modality, "direction": "synthetic_to_real",
                            "model": model, "seed": int(seed), "celltype": lab,
                            "recall": rec,
                        })
                if keep_syn.sum() >= 10 and keep_real.sum() >= 10:
                    res = transfer(xr[keep_real], y_real[keep_real],
                                   xs[keep_syn], y_syn[keep_syn], model)
                    transfer_rows.append({
                        "modality": modality, "direction": "real_to_synthetic",
                        "model": model, "seed": int(seed), **{k: res[k] for k in
                        ("accuracy", "balanced_accuracy", "macro_f1")},
                        "n_train": int(keep_real.sum()), "n_test": int(keep_syn.sum()),
                    })

            # ---- label-transfer control: permuted labels -----------------
            y_perm = RNG.permutation(y_syn)
            for model in models:
                res = transfer(xs, y_perm, xr, y_ct[r_idx], model)
                transfer_rows.append({
                    "modality": modality, "direction": "permuted_control",
                    "model": model, "seed": int(seed), **{k: res[k] for k in
                    ("accuracy", "balanced_accuracy", "macro_f1")},
                    "n_train": len(y_perm), "n_test": len(y_ct[r_idx]),
                })

    disc = pd.DataFrame(disc_rows)
    if "detail" not in disc.columns:
        disc["detail"] = ""
    disc["detail"] = disc["detail"].fillna("")
    # distance from chance: a faithful simulator sits near 0, either direction
    disc["auc_dist"] = (disc["auc"] - 0.5).abs()
    disc.to_csv(out / "discriminator.csv", index=False)
    pd.DataFrame(transfer_rows).to_csv(out / "label_transfer.csv", index=False)
    pd.DataFrame(per_class_rows).to_csv(out / "per_class.csv", index=False)

    summary = {
        "mode": args.mode,
        "seeds": [int(s) for s in seeds],
        "models": models,
        "n_holdout": len(holdout),
        "discriminator_auc": (
            disc.groupby(["modality", "comparison", "model"])
            .agg(auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                 auc_dist_mean=("auc_dist", "mean"), n=("auc", "size"))
            .reset_index().to_dict(orient="records")
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote {out}/discriminator.csv, label_transfer.csv, per_class.csv, summary.json")
    print(disc.groupby(["modality", "comparison", "model"])["auc"]
          .agg(["mean", "std"]).round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
