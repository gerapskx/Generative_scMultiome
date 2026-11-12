#!/usr/bin/env python3
"""Synthetic-to-real utility proof of concept.

The fidelity module asks whether synthetic cells look real. This module asks a
different question: are they useful for a task? A classifier is trained on synthetic
cells generated from training donors and tested on held-out real donors.

Baselines
    real_trained     real training donors -> held-out donors (upper bound)
    synthetic        synthetic training cells -> held-out donors (the PoC)
    naive_poisson    marginal-only synthetic -> held-out donors (lower bound)
    permuted_labels  synthetic with shuffled condition labels -> chance
    majority         majority class -> chance

The learning curve varies the number of synthetic training cells at a fixed donor
split. Donor-disjoint evaluation is only defined from n08 onward; at n02 the fold
table is marked `not_estimable`.

Usage:
    python src/synthetic_utility.py [--ingest results/ingest]
                                   [--features results/features]
                                   [--obs results/cluster/obs_annotated.csv]
                                   [--sim results/utility/sim]
                                   [--out results/utility]

Outputs (out/):
    utility.csv          one row per fold x modality x size x seed x method
    learning_curve.csv    mean/sd across seeds
    summary.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    roc_auc_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fidelity import (  # noqa: E402
    find,
    load_features,
    naive_poisson,
    normalise,
    read_mtx,
)


def metrics(y_true: np.ndarray, proba: np.ndarray, pred: np.ndarray) -> dict:
    out = {"auc": float("nan"), "balanced_accuracy": float("nan"),
           "macro_f1": float("nan"), "brier": float("nan")}
    if len(np.unique(y_true)) == 2:
        out["auc"] = float(roc_auc_score(y_true, proba))
        out["brier"] = float(brier_score_loss(y_true, proba))
    out["balanced_accuracy"] = float(balanced_accuracy_score(y_true, pred))
    out["macro_f1"] = float(f1_score(y_true, pred, average="macro"))
    return out


def fit_predict(
    x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray, y_test: np.ndarray
) -> dict:
    model = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=3000, C=1.0, class_weight="balanced"),
    )
    model.fit(x_train, y_train)
    proba = model.predict_proba(x_test)[:, 1]
    pred = model.predict(x_test)
    return metrics(y_test, proba, pred)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--sim", default="results/utility/sim")
    ap.add_argument("--out", default="results/utility")
    ap.add_argument("--label", default="condition", choices=["condition", "celltype"])
    ap.add_argument("--folds", type=int, default=3)
    ap.add_argument("--sizes", default="500,1000,2000,5000,10000")
    ap.add_argument("--seed", type=int, default=20261110)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    sizes = [int(v) for v in args.sizes.split(",") if v]

    folds_path = Path(args.sim) / "folds.csv"
    manifest_path = Path(args.sim) / "manifest.csv"
    if not folds_path.exists() or not manifest_path.exists():
        raise SystemExit(f"missing {folds_path} or {manifest_path}; run utility_sim first")
    folds = pd.read_csv(folds_path, keep_default_na=False)
    manifest = pd.read_csv(manifest_path, keep_default_na=False)

    if folds.empty or "status" not in folds.columns or not (folds["status"] == "ok").any():
        payload = {
            "label": args.label,
            "status": "not_estimable",
            "reason": "at least one condition has fewer than two donors",
            "rows": [],
        }
        (out / "summary.json").write_text(json.dumps(payload, indent=2))
        pd.DataFrame().to_csv(out / "utility.csv", index=False)
        pd.DataFrame().to_csv(out / "learning_curve.csv", index=False)
        print("utility: not estimable (fewer than two donors per condition)")
        return 0

    obs = pd.read_csv(args.obs, dtype={"cell": str})
    obs["condition"] = obs["group"].astype(str)
    meta = obs.set_index("cell")

    rows: list[dict] = []
    for _, fold_row in folds[folds["status"] == "ok"].iterrows():
        fold = int(fold_row["fold"])
        train_donors = set(str(fold_row["train_donors"]).split(";"))
        holdout_donors = set(str(fold_row["holdout_donors"]).split(";"))
        train_cells = obs.loc[obs["donor"].isin(train_donors), "cell"].tolist()
        holdout_cells = obs.loc[obs["donor"].isin(holdout_donors), "cell"].tolist()
        if not train_cells or not holdout_cells:
            continue

        fold_dir = Path(args.sim) / f"fold{fold}"
        sim_obs = pd.read_csv(fold_dir / "obs.csv", dtype={"cell": str})
        y_syn_all = (sim_obs["condition"].astype(str) == "PD").astype(int).to_numpy()
        y_real_train = (meta.loc[train_cells, "condition"] == "PD").astype(int).to_numpy()
        y_real_test = (meta.loc[holdout_cells, "condition"] == "PD").astype(int).to_numpy()
        if len(np.unique(y_real_train)) < 2 or len(np.unique(y_real_test)) < 2:
            continue

        for modality in ("rna", "atac"):
            feats = load_features(Path(args.features) / f"{modality}_features.tsv")
            feat_all = load_features(ingest / f"{modality}_features.tsv")
            fpos = {f: i for i, f in enumerate(feat_all)}
            rows_idx = np.array([fpos[f] for f in feats])
            cpos = {c: i for i, c in enumerate(obs["cell"])}
            train_idx = np.array([cpos[c] for c in train_cells])
            test_idx = np.array([cpos[c] for c in holdout_cells])

            counts = read_mtx(find(ingest / f"{modality}_counts.mtx"))
            x_train = normalise(counts[rows_idx][:, train_idx], modality)
            x_test = normalise(counts[rows_idx][:, test_idx], modality)

            # ---- real-trained upper bound --------------------------------
            res = fit_predict(x_train, y_real_train, x_test, y_real_test)
            rows.append({
                "fold": fold, "modality": modality, "label": args.label,
                "synthetic_size": 0, "seed": 0, "method": "real_trained",
                "n_train": len(y_real_train), "n_test": len(y_real_test),
                **res,
            })

            # ---- majority baseline --------------------------------------
            majority = int(np.bincount(y_syn_all).argmax())
            pred = np.full(len(y_real_test), majority)
            proba = np.full(len(y_real_test), float(majority))
            rows.append({
                "fold": fold, "modality": modality, "label": args.label,
                "synthetic_size": 0, "seed": 0, "method": "majority",
                "n_train": len(y_real_train), "n_test": len(y_real_test),
                **metrics(y_real_test, proba, pred),
            })

            # ---- naive Poisson control ----------------------------------
            x_naive = naive_poisson(counts[rows_idx][:, train_idx], y_real_train,
                                    modality)
            res = fit_predict(x_naive, y_real_train, x_test, y_real_test)
            rows.append({
                "fold": fold, "modality": modality, "label": args.label,
                "synthetic_size": len(y_real_train), "seed": 0,
                "method": "naive_poisson",
                "n_train": len(y_real_train), "n_test": len(y_real_test),
                **res,
            })

            # ---- synthetic simulations ---------------------------------
            for _, sim_row in manifest[(manifest["fold"] == fold) &
                                       (manifest["modality"] == modality)].iterrows():
                seed = int(sim_row["seed"])
                sim = read_mtx(Path(sim_row["path"]))
                if sim.shape[0] != len(feats):
                    raise SystemExit(f"{sim_row['path']}: {sim.shape[0]} rows")
                x_syn = normalise(sim, modality)
                if len(x_syn) != len(y_syn_all):
                    raise SystemExit(f"{sim_row['path']}: cell count != sim obs")
                for size in sizes:
                    k = min(size, len(x_syn))
                    take = np.sort(rng.choice(len(x_syn), k, replace=False))
                    xs, ys = x_syn[take], y_syn_all[take]
                    if len(np.unique(ys)) < 2:
                        continue

                    res = fit_predict(xs, ys, x_test, y_real_test)
                    rows.append({
                        "fold": fold, "modality": modality, "label": args.label,
                        "synthetic_size": size, "seed": seed,
                        "method": "synthetic",
                        "n_train": k, "n_test": len(y_real_test),
                        **res,
                    })

                    y_perm = rng.permutation(ys)
                    res = fit_predict(xs, y_perm, x_test, y_real_test)
                    rows.append({
                        "fold": fold, "modality": modality, "label": args.label,
                        "synthetic_size": size, "seed": seed,
                        "method": "permuted_labels",
                        "n_train": k, "n_test": len(y_real_test),
                        **res,
                    })

    table = pd.DataFrame(rows)
    table.to_csv(out / "utility.csv", index=False)

    if table.empty:
        curve = pd.DataFrame()
        summary = {"label": args.label, "status": "no_estimable_folds", "rows": []}
    else:
        curve = (table.groupby(["modality", "method", "synthetic_size"])
                 .agg(auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                      bal_mean=("balanced_accuracy", "mean"),
                      bal_sd=("balanced_accuracy", "std"),
                      n=("auc", "size")).reset_index())
        curve.to_csv(out / "learning_curve.csv", index=False)
        summary = {
            "label": args.label,
            "status": "ok",
            "folds": int(table["fold"].nunique()),
            "sizes": sizes,
            "headline": curve.round(4).to_dict(orient="records"),
        }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))

    print(f"wrote {out}/utility.csv, learning_curve.csv, summary.json")
    if not table.empty:
        print(curve.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
