#!/usr/bin/env python3
"""Rank cell types by how strongly the disease state affects them.

For every cell type and modality this script reports:

  de_burden        significant condition features / testable features
  mean_abs_lfc     mean |log2FC| among significant features (PD vs Aged)
  donor_auc        condition classifier on held-out donors (not cells)
  eta2_condition   mean share of feature variance explained by condition
  eta2_donor       the same for donor; the ratio shows confounding
  reversal_*       counterfactual magnitude from the model branch
  composition_*    PD vs Aged share of the cell type (descriptive)

The donor-held-out AUC is undefined when a condition has only one donor, which is
the case at n02. It is the metric that makes the scaling series interpretable: a
high cell-level AUC with one donor per condition is donor fingerprinting.

Usage:
    python src/celltype_impact.py [--ingest results/ingest]
                                 [--features results/features]
                                 [--split results/split/split.csv]
                                 [--obs results/cluster/obs_annotated.csv]
                                 [--reversal results/reversal]
                                 [--out results/impact]

Outputs (out/):
    celltype_impact.csv   one row per modality x celltype
    de_features.csv       per-feature condition test within each celltype
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
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fidelity import find, load_features, normalise, read_mtx  # noqa: E402
from hypothesis import bh, de_test  # noqa: E402


def eta_squared(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    """Between-group sum of squares / total sum of squares, per column."""
    grand = values.mean(axis=0, keepdims=True)
    total = ((values - grand) ** 2).sum(axis=0)
    between = np.zeros(values.shape[1])
    for lab in np.unique(labels):
        block = values[labels == lab]
        between += block.shape[0] * (block.mean(axis=0) - grand[0]) ** 2
    out = np.divide(between, total, out=np.zeros_like(between), where=total > 0)
    return out


def donor_auc(
    x: np.ndarray, y: np.ndarray, donors: np.ndarray, folds: int
) -> tuple[float, float, int]:
    """Leave-donor-out AUC / balanced accuracy for a condition classifier."""
    n_donors = len(np.unique(donors))
    if n_donors < 3 or min(np.bincount(y, minlength=2)) < 2:
        return float("nan"), float("nan"), 0
    splitter = GroupKFold(n_splits=min(folds, n_donors))
    aucs, bals, weights = [], [], []
    for train, test in splitter.split(x, y, groups=donors):
        if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
            continue
        model = make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced"),
        )
        model.fit(x[train], y[train])
        proba = model.predict_proba(x[test])[:, 1]
        aucs.append(roc_auc_score(y[test], proba))
        bals.append(balanced_accuracy_score(y[test], model.predict(x[test])))
        weights.append(len(test))
    if not aucs:
        return float("nan"), float("nan"), 0
    w = np.asarray(weights, dtype=float)
    return float(np.average(aucs, weights=w)), float(np.average(bals, weights=w)), len(aucs)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--split", default="results/split/split.csv")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--reversal", default="results/reversal")
    ap.add_argument("--out", default="results/impact")
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--has-reversal", default="true")
    ap.add_argument("--max-eta-features", type=int, default=500)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    obs = pd.read_csv(args.obs, dtype={"cell": str, "cluster": str})
    split = pd.read_csv(args.split, dtype={"cell": str})
    obs = obs.merge(split[["cell", "split"]], on="cell", how="left")
    obs["condition"] = obs["group"].astype(str)
    obs["is_pd"] = (obs["condition"] == "PD").astype(int)

    reversal = None
    rev_path = Path(args.reversal) / "model_delta.csv"
    if args.has_reversal.lower() in ("true", "yes", "1") and rev_path.exists():
        reversal = pd.read_csv(rev_path)

    impact_rows, de_rows = [], []
    for modality in ("rna", "atac"):
        feats = load_features(Path(args.features) / f"{modality}_features.tsv")
        counts = read_mtx(find(ingest / f"{modality}_counts.mtx"))
        feat_all = load_features(ingest / f"{modality}_features.tsv")
        fpos = {f: i for i, f in enumerate(feat_all)}
        missing = [f for f in feats if f not in fpos]
        if missing:
            raise SystemExit(f"{modality}: {len(missing)} selected features absent")
        rows_idx = np.array([fpos[f] for f in feats])
        cells = obs["cell"].tolist()
        cpos = {c: i for i, c in enumerate(cells)}
        cols_idx = np.array([cpos[c] for c in cells])
        x_all = normalise(counts[rows_idx][:, cols_idx], modality)
        print(f"[{modality}] {x_all.shape[0]} cells x {x_all.shape[1]} features",
              flush=True)

        for ct, sub in obs.groupby("celltype", sort=True):
            n_pd = int((sub["condition"] == "PD").sum())
            n_aged = int((sub["condition"] == "Aged").sum())
            if min(n_pd, n_aged) < args.min_cells:
                print(f"  skip {ct}: PD={n_pd} Aged={n_aged}", flush=True)
                continue

            idx = np.array([cpos[c] for c in sub["cell"]])
            x = x_all[idx]
            y = sub["is_pd"].to_numpy()
            donors = sub["donor"].astype(str).to_numpy()
            donor_levels = pd.unique(donors)

            # ---- condition DE within the cell type -------------------------
            res = de_test(x, y == 1, y == 0, "wilcox")
            testable = (x[y == 1].var(axis=0) > 0) | (x[y == 0].var(axis=0) > 0)
            sig = (res["q"] <= args.fdr) & testable
            n_sig = int(sig.sum())
            lfc = res["log2fc"]
            abs_sig = np.abs(lfc[sig])
            de_rows.append(pd.DataFrame({
                "modality": modality, "celltype": ct, "feature": feats,
                "log2fc": lfc, "p": res["p"], "q": res["q"],
                "testable": testable, "significant": sig,
                "mean_pd": x[y == 1].mean(axis=0),
                "mean_aged": x[y == 0].mean(axis=0),
            }))

            # ---- donor-held-out classifier --------------------------------
            auc, bal, n_folds = donor_auc(x, y, donors, args.folds)

            # ---- variance decomposition ----------------------------------
            var = x.var(axis=0)
            top = np.argsort(var)[::-1][: min(args.max_eta_features, x.shape[1])]
            eta_cond = eta_squared(x[:, top], y)
            eta_don = eta_squared(x[:, top], donors)
            eta_cond_mean = float(np.nanmean(eta_cond))
            eta_don_mean = float(np.nanmean(eta_don))
            denom = eta_cond_mean + eta_don_mean
            condition_share = eta_cond_mean / denom if denom > 0 else float("nan")

            # ---- counterfactual reversal ---------------------------------
            rev_mean = rev_frac = rev_median = float("nan")
            if reversal is not None:
                rsub = reversal[(reversal["modality"] == modality) &
                                 (reversal["celltype"] == ct)]
                if not rsub.empty:
                    rev_mean = float(rsub["log2fc"].abs().mean())
                    rev_median = float(rsub["log2fc"].abs().median())
                    q95 = float(rsub["log2fc_permuted"].abs().quantile(0.95))
                    if np.isfinite(q95):
                        rev_frac = float((rsub["log2fc"].abs() > q95).mean())

            impact_rows.append({
                "modality": modality, "celltype": ct,
                "n_cells": len(sub), "n_pd": n_pd, "n_aged": n_aged,
                "n_donors": len(donor_levels),
                "n_donors_pd": len(np.unique(donors[y == 1])),
                "n_donors_aged": len(np.unique(donors[y == 0])),
                "n_features": x.shape[1],
                "n_testable": int(testable.sum()),
                "n_de_sig": n_sig,
                "de_burden": n_sig / max(int(testable.sum()), 1),
                "mean_abs_lfc_sig": float(abs_sig.mean()) if abs_sig.size else float("nan"),
                "max_abs_lfc_sig": float(abs_sig.max()) if abs_sig.size else float("nan"),
                "donor_auc": auc,
                "donor_balanced_accuracy": bal,
                "donor_folds": n_folds,
                "eta2_condition": eta_cond_mean,
                "eta2_donor": eta_don_mean,
                "condition_variance_share": condition_share,
                "pd_fraction": n_pd / len(sub),
                "reversal_mean_abs_lfc": rev_mean,
                "reversal_median_abs_lfc": rev_median,
                "reversal_frac_above_perm_q95": rev_frac,
                "condition_equals_donor": bool(
                    len(donor_levels) == 2 and
                    len(np.unique(donors[y == 1])) == 1 and
                    len(np.unique(donors[y == 0])) == 1
                ),
            })
            print(f"  {ct}: n={len(sub)} sig={n_sig} "
                  f"auc={auc:.3f} eta_cond={eta_cond_mean:.3f}", flush=True)

    impact = pd.DataFrame(impact_rows)
    impact.to_csv(out / "celltype_impact.csv", index=False)
    pd.concat(de_rows, ignore_index=True).to_csv(out / "de_features.csv", index=False)

    summary = {
        "fdr": args.fdr,
        "min_cells_per_stratum": args.min_cells,
        "folds": args.folds,
        "modalities": ["rna", "atac"],
        "n_celltypes": int(impact["celltype"].nunique()) if not impact.empty else 0,
        "headline": (
            impact.sort_values(["modality", "de_burden"], ascending=[True, False])
            .groupby("modality")
            .head(1)[["modality", "celltype", "de_burden", "donor_auc"]]
            .to_dict(orient="records") if not impact.empty else []
        ),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {out}/celltype_impact.csv, de_features.csv, summary.json")
    if not impact.empty:
        print(impact.round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
