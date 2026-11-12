#!/usr/bin/env python3
"""Stage 9b - counterfactual reversal, classifier branch.

Companion to `R/04_reversal.R`. The R branch asks what the fitted scDesign3
marginal model predicts for a PD nucleus under the Aged condition. This branch
asks a different question with a different model class: train an Aged-vs-PD
classifier on the fit cells of one cell type, then score every selected feature
by how much it pushes a held-out PD cell toward the Aged class.

For a logistic regression with standardised features, the contribution of
feature j to the log-odds of Aged for cell i is exactly `coef_j * z_ij`. The
per-feature score is the mean over held-out PD cells, so a positive score means
"higher values move the cell toward Aged".

Controls
    reverse_direction   score held-out Aged cells toward PD; the ranking should
                        anticorrelate with the PD -> Aged ranking
    permuted_condition  retrain with condition labels shuffled; contributions
                        must collapse toward zero

Only cell types selected by the model branch (>= min_cells_per_stratum in both
conditions) are scored. In the 2-donor pilot that is Oligodendrocyte only.

Usage:
    python src/reversal.py [--ingest results/ingest]
                           [--features results/features]
                           [--split results/split/split.csv]
                           [--obs results/cluster/obs_annotated.csv]
                           [--reversal results/reversal]
                           [--out results/reversal]
                           [--model logistic_regression]
                           [--folds 5] [--seed 20261110]

Outputs (results/reversal/):
    classifier_contribution.csv  feature, mean contribution, reverse, permuted
    combined.csv                 model and classifier scores side by side
    top_features.csv             top features by each method
    summary.json                 agreement between the two methods + controls
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------------------------- input --
def find(path_stem: Path) -> Path:
    gz = Path(str(path_stem) + ".gz")
    if gz.exists():
        return gz
    if path_stem.exists():
        return path_stem
    raise FileNotFoundError(f"neither {path_stem} nor {gz} exists")


def read_mtx(path: Path) -> sp.csr_matrix:
    if path.suffix == ".gz":
        import gzip

        with gzip.open(path, "rb") as fh:
            return sp.csr_matrix(scipy.io.mmread(fh))
    with path.open("rb") as fh:
        return sp.csr_matrix(scipy.io.mmread(fh))


def load_features(path: Path) -> list[str]:
    return [line for line in path.read_text().splitlines() if line]


def normalise(mat: sp.csr_matrix, modality: str) -> np.ndarray:
    """cells x features. RNA: log1p CPM. ATAC: TF-IDF (same as fidelity.py)."""
    mat = mat.astype(np.float64).T.tocsr()
    lib = np.asarray(mat.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    if modality == "rna":
        return np.log1p(mat.multiply(1e6 / lib[:, None]).toarray())
    tf = mat.multiply(1.0 / lib[:, None]).tocsr()
    df = np.asarray((mat > 0).sum(axis=0)).ravel()
    idf = np.log1p(mat.shape[0] / np.maximum(df, 1))
    return tf.multiply(idf[None, :]).toarray()


def make_model(name: str):
    if name == "logistic_regression":
        return make_pipeline(
            StandardScaler(),
            LogisticRegression(max_iter=5000, C=1.0, class_weight="balanced"),
        )
    raise SystemExit(f"unknown model: {name} (contributions need a linear model)")


def contributions(model, x: np.ndarray) -> np.ndarray:
    """Mean coef * z over rows of x; positive = toward the positive class."""
    scaler = model.named_steps["standardscaler"]
    clf = model.named_steps["logisticregression"]
    z = scaler.transform(x)
    return (z * clf.coef_.ravel()[None, :]).mean(axis=0)


# ---------------------------------------------------------------------- main --
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--split", default="results/split/split.csv")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--reversal", default="results/reversal",
                    help="directory holding model_delta.csv from R/04_reversal.R")
    ap.add_argument("--out", default="results/reversal")
    ap.add_argument("--model", default="logistic_regression")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--target", default="PD")
    ap.add_argument("--reference", default="Aged")
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--min-target-cells", type=int, default=10)
    ap.add_argument("--top-n", default="50,100,200")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    top_ns = [int(v) for v in args.top_n.split(",") if v]

    obs = pd.read_csv(args.obs, dtype={"cell": str})
    split = pd.read_csv(args.split, dtype={"cell": str})
    meta = obs.set_index("cell")

    model_delta = pd.read_csv(Path(args.reversal) / "model_delta.csv")
    if model_delta.empty:
        raise SystemExit("model_delta.csv is empty; run R/04_reversal.R first")

    fit_cells = set(split.loc[split["split"] == "fit", "cell"])
    holdout_cells = set(split.loc[split["split"] == "holdout", "cell"])
    all_cells = obs["cell"].tolist()

    rows, combined_rows, summary = [], [], []

    for (modality, celltype), sub in model_delta.groupby(["modality", "celltype"]):
        feats = sub["feature"].tolist()
        print(f"[{modality}/{celltype}] {len(feats)} features", flush=True)

        counts = read_mtx(find(ingest / f"{modality}_counts.mtx"))
        feat_all = load_features(ingest / f"{modality}_features.tsv")
        fpos = {f: i for i, f in enumerate(feat_all)}
        cpos = {c: i for i, c in enumerate(all_cells)}
        rows_idx = np.array([fpos[f] for f in feats])
        cols_idx = np.array([cpos[c] for c in all_cells])
        x = normalise(counts[rows_idx][:, cols_idx], modality)
        x = pd.DataFrame(x, index=all_cells, columns=feats)

        in_ct = meta.index[meta["celltype"] == celltype]
        fit_ct = [c for c in in_ct if c in fit_cells]
        hold_ct = [c for c in in_ct if c in holdout_cells]
        y_fit = (meta.loc[fit_ct, "group"] == args.reference).astype(int).to_numpy()

        if len(np.unique(y_fit)) < 2:
            print(f"  skip {modality}/{celltype}: one condition among fit cells")
            continue

        # ---- observed model: PD -> Aged -----------------------------------
        model = make_model(args.model).fit(x.loc[fit_ct].to_numpy(), y_fit)
        cv = StratifiedKFold(n_splits=max(2, min(args.folds, np.bincount(y_fit).min())),
                             shuffle=True, random_state=0)
        pred = cross_val_predict(make_model(args.model), x.loc[fit_ct].to_numpy(),
                                 y_fit, cv=cv)
        cv_bal_acc = float(balanced_accuracy_score(y_fit, pred))

        pd_cells = [c for c in hold_ct
                    if meta.loc[c, "group"] == args.target]
        aged_cells = [c for c in hold_ct
                      if meta.loc[c, "group"] == args.reference]
        if len(pd_cells) < args.min_target_cells:
            print(f"  skip {modality}/{celltype}: {len(pd_cells)} holdout PD cells")
            continue

        contrib = contributions(model, x.loc[pd_cells].to_numpy())

        # ---- control: reverse direction (Aged -> PD) ----------------------
        # scoring Aged cells toward PD is the negative of the same linear score
        contrib_reverse = -contributions(model, x.loc[aged_cells].to_numpy())

        # ---- control: permuted condition ----------------------------------
        y_perm = rng.permutation(y_fit)
        model_perm = make_model(args.model).fit(x.loc[fit_ct].to_numpy(), y_perm)
        contrib_perm = contributions(model_perm, x.loc[pd_cells].to_numpy())

        clf_rows = pd.DataFrame({
            "modality": modality, "celltype": celltype, "feature": feats,
            "contribution": contrib,
            "contribution_reverse": contrib_reverse,
            "contribution_permuted": contrib_perm,
            "n_fit": len(fit_ct), "n_pd_holdout": len(pd_cells),
            "n_aged_holdout": len(aged_cells),
            "cv_balanced_accuracy": cv_bal_acc,
            "model": args.model,
        })
        rows.append(clf_rows)

        # ---- agreement with the model branch ------------------------------
        merged = sub[["feature", "log2fc", "log2fc_permuted"]].merge(
            clf_rows[["feature", "contribution", "contribution_reverse",
                      "contribution_permuted"]],
            on="feature", how="inner")
        merged["modality"] = modality
        merged["celltype"] = celltype
        combined_rows.append(merged)

        rho = merged[["log2fc", "contribution"]].corr(method="spearman").iloc[0, 1]
        rho_rev = merged[["log2fc", "contribution_reverse"]].corr(
            method="spearman").iloc[0, 1]
        rho_perm = merged[["log2fc_permuted", "contribution_permuted"]].corr(
            method="spearman").iloc[0, 1]

        overlaps = {}
        for n in top_ns:
            k = min(n, len(merged))
            a = set(merged.reindex(merged["log2fc"].abs().sort_values(
                ascending=False).index)["feature"].head(k))
            b = set(merged.reindex(merged["contribution"].abs().sort_values(
                ascending=False).index)["feature"].head(k))
            overlaps[f"top{n}_overlap"] = len(a & b) / max(k, 1)

        summary.append({
            "modality": modality, "celltype": celltype,
            "n_features": len(merged),
            "n_fit": len(fit_ct), "n_pd_holdout": len(pd_cells),
            "cv_balanced_accuracy": cv_bal_acc,
            "spearman_model_vs_classifier": float(rho),
            "spearman_reverse_direction": float(rho_rev),
            "spearman_permuted": float(rho_perm),
            "classifier_permuted_max_abs": float(np.max(np.abs(contrib_perm))),
            "model_permuted_max_abs": float(np.max(np.abs(merged["log2fc_permuted"]))),
            **overlaps,
        })
        print(f"  cv bal acc {cv_bal_acc:.3f}; rho {rho:.3f}; "
              f"top100 overlap {overlaps.get('top100_overlap', float('nan')):.2f}",
              flush=True)

    clf = pd.concat(rows, ignore_index=True)
    combined = pd.concat(combined_rows, ignore_index=True)
    clf.to_csv(out / "classifier_contribution.csv", index=False)
    combined.to_csv(out / "combined.csv", index=False)

    top = []
    for (modality, celltype), sub in combined.groupby(["modality", "celltype"]):
        for method, col in (("model", "log2fc"), ("classifier", "contribution")):
            t = sub.reindex(sub[col].abs().sort_values(ascending=False).index).head(200)
            top.append(t.assign(method=method, rank=np.arange(1, len(t) + 1)))
    pd.DataFrame(pd.concat(top, ignore_index=True)).to_csv(
        out / "top_features.csv", index=False)

    payload = {
        "target": args.target, "reference": args.reference,
        "model": args.model, "top_n": top_ns,
        "strata": summary,
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2))
    print(f"wrote {out}/classifier_contribution.csv, combined.csv, "
          f"top_features.csv, summary.json")
    if summary:
        print(pd.DataFrame(summary).round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
