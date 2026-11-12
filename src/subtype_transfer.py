#!/usr/bin/env python3
"""Synthetic-to-real subtype transfer.

For each target cell type and modality this script builds a 4 x 4 transfer matrix
over four domains:

    real          real cells from the donor-disjoint split
    syn_alt       scDesign3 `subtype + condition` simulation
    syn_null      scDesign3 null refit (`condition` only)
    syn_null_zero scDesign3 alt fit with subtype coefficients zeroed

Training on synthetic and testing on real is the direction that matters: a simulator
that preserves the subtype axis should produce cells a real-trained problem can still
classify, and a synthetic-trained classifier should still recognise real cells.
The reverse direction (real -> synthetic) is reported alongside it so the asymmetry
is visible.

Evaluation is donor-disjoint for `real` and seed-paired for synthetic domains. Every
pair is reported overall, by condition and by donor (donors only exist for real test
cells). Controls:

    condition_only   a one-feature condition classifier -> subtype (confound check)
    permuted_labels  train labels shuffled -> chance
    majority         majority training class -> chance
    donor_only       donor-mean profiles, real -> real (donor-level ceiling)

Usage:
    python src/subtype_transfer.py --root results/scale/n24/subtype \
        --ingest results/scale/n24/ingest --out results/scale/n24/subtype/transfer
"""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (balanced_accuracy_score, brier_score_loss,
                             f1_score, roc_auc_score)
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

DOMAINS = ["real", "syn_alt", "syn_null", "syn_null_zero"]
MODES = {"syn_alt": "alt", "syn_null": "null", "syn_null_zero": "null_zero"}


def binary(subtype) -> np.ndarray:
    """Map s1/s2 (or 0/1) subtype labels to 0/1 with s2 == 1."""
    arr = np.asarray(subtype)
    if arr.dtype.kind in "iu" or arr.dtype.kind == "f":
        return arr.astype(int)
    return (arr.astype(str) == "s2").astype(int)


def read_mtx(stem: Path) -> scipy.sparse.csr_matrix:
    gz = Path(str(stem) + ".gz")
    path = gz if gz.exists() else stem
    with (gzip.open(path, "rb") if str(path).endswith(".gz") else open(path, "rb")) as fh:
        return scipy.io.mmread(fh).tocsr()


def load_real(root: Path, ingest: Path, celltype: str, modality: str):
    labels = pd.read_csv(root / celltype / modality / "labels.csv",
                         dtype={"cell": str, "donor": str, "group": str})
    obs = pd.read_csv(ingest / "obs.csv", dtype={"cell": str})
    counts = read_mtx(ingest / f"{modality}_counts_sel.mtx")
    if counts.shape[1] != len(obs):
        raise SystemExit(f"{modality}: matrix columns != ingest obs rows")
    col = obs.set_index("cell").index.get_indexer(labels["cell"])
    if (col < 0).any():
        raise SystemExit(f"{modality}: labels contain cells absent from ingest")
    x = counts[:, col].T.tocsr()
    return {
        "x": x,
        "y": binary(labels["subtype"]),
        "condition": labels["group"].to_numpy(),
        "donor": labels["donor"].to_numpy(),
        "split": labels["split"].to_numpy(),
    }


def load_syn(root: Path, celltype: str, modality: str, mode: str, seed: int):
    base = root / celltype / modality / "simulation" / mode
    obs = pd.read_csv(base / "obs.csv", dtype={"cell": str})
    counts = read_mtx(base / f"seed{seed}_counts.mtx")
    if counts.shape[1] != len(obs):
        raise SystemExit(f"{mode}/{seed}: matrix columns != obs rows")
    return {
        "x": counts.T.tocsr(),
        "y": binary(obs["subtype"]),
        "condition": obs["condition"].astype(str).to_numpy(),
        "donor": np.array(["synthetic"] * len(obs)),
        "split": np.array(["synthetic"] * len(obs)),
    }


def normalise(train: np.ndarray, test: np.ndarray, modality: str):
    train = np.asarray(train, dtype=np.float64)
    test = np.asarray(test, dtype=np.float64)
    lib_train = train.sum(axis=1)
    lib_test = test.sum(axis=1)
    lib_train[lib_train <= 0] = 1.0
    lib_test[lib_test <= 0] = 1.0
    if modality == "rna":
        train = np.log1p(train / lib_train[:, None] * 1e4)
        test = np.log1p(test / lib_test[:, None] * 1e4)
    else:
        tf_train = train / lib_train[:, None]
        tf_test = test / lib_test[:, None]
        # IDF is fit on the training domain only.
        idf = np.log1p(train.shape[0] / np.maximum((train > 0).sum(axis=0), 1))
        train = tf_train * idf
        test = tf_test * idf
    return train, test


def subsample(x, y, condition, donor, n, seed, include_donor=False):
    if n <= 0 or len(y) <= n:
        return x, y, condition, donor
    strata = np.char.add(np.char.add(y.astype(str), "|"), condition.astype(str))
    if include_donor:
        strata = np.char.add(np.char.add(strata, "|"), donor.astype(str))
    rng = np.random.default_rng(seed)
    keep = []
    for s in np.unique(strata):
        idx = np.flatnonzero(strata == s)
        take = max(1, int(round(n * len(idx) / len(y))))
        take = min(take, len(idx))
        keep.extend(rng.choice(idx, size=take, replace=False).tolist())
    keep = np.array(sorted(set(keep)))
    if len(keep) > n:
        keep = np.sort(rng.choice(keep, size=n, replace=False))
    return x[keep], y[keep], condition[keep], donor[keep]


def score(y_true, prob, pred, n_train, n_test):
    out = {}
    try:
        out["auc"] = float(roc_auc_score(y_true, prob))
    except ValueError:
        out["auc"] = np.nan
    out["balanced_accuracy"] = float(balanced_accuracy_score(y_true, pred))
    out["macro_f1"] = float(f1_score(y_true, pred, average="macro"))
    try:
        out["brier"] = float(brier_score_loss(y_true, prob))
    except ValueError:
        out["brier"] = np.nan
    out["n_train"] = int(n_train)
    out["n_test"] = int(n_test)
    return out


def evaluate(y_train, x_train, y_test, x_test, condition_test, donor_test,
              model: str, seed: int):
    y_train = np.asarray(y_train).astype(int)
    y_test = np.asarray(y_test).astype(int)
    if len(np.unique(y_train)) < 2:
        return None
    if x_train.ndim == 1:
        x_train = x_train[:, None]
        x_test = x_test[:, None]
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train)
    x_test = scaler.transform(x_test)
    clf = make_classifier(model, seed)
    clf.fit(x_train, y_train)
    prob = clf.predict_proba(x_test)[:, 1]
    pred = (prob >= 0.5).astype(int)
    metrics = score(y_test, prob, pred, len(y_train), len(y_test))
    return metrics, prob, pred


def make_classifier(model: str, seed: int):
    if model == "random_forest":
        return RandomForestClassifier(
            n_estimators=100, max_features="sqrt",
            class_weight="balanced_subsample", n_jobs=4, random_state=seed,
        )
    return LogisticRegression(
        C=0.1, class_weight="balanced", max_iter=2000,
        solver="liblinear", random_state=seed,
    )


def stratified_take(idx, y, condition, donor, n, seed, include_donor):
    if n <= 0 or len(idx) <= n:
        return idx
    strata = np.char.add(np.char.add(y.astype(str), "|"), condition.astype(str))
    if include_donor:
        strata = np.char.add(np.char.add(strata, "|"), donor.astype(str))
    rng = np.random.default_rng(seed)
    keep = []
    for s in np.unique(strata[idx]):
        sub = idx[strata[idx] == s]
        take = max(1, int(round(n * len(sub) / len(idx))))
        take = min(take, len(sub))
        keep.extend(rng.choice(sub, size=take, replace=False).tolist())
    keep = np.array(sorted(set(keep)))
    if len(keep) > n:
        keep = np.sort(rng.choice(keep, size=n, replace=False))
    return keep


def add_rows(rows, celltype, modality, model, train_domain, test_domain, seed,
             res, y_test, condition_test, donor_test, min_stratum):
    if res is None:
        return
    metrics, prob, pred = res
    for metric, value in metrics.items():
        if metric in ("n_train", "n_test"):
            continue
        rows.append({
            "celltype": celltype, "modality": modality, "model": model,
            "train_domain": train_domain, "test_domain": test_domain,
            "seed": seed, "metric": metric, "stratum": "overall",
            "value": value, "n_train": metrics["n_train"],
            "n_test": metrics["n_test"],
        })
    if prob is None:
        return
    for cond in np.unique(condition_test):
        mask = condition_test == cond
        if mask.sum() < min_stratum or len(np.unique(y_test[mask])) < 2:
            continue
        m = score(y_test[mask], prob[mask], pred[mask],
                  metrics["n_train"], int(mask.sum()))
        for metric, value in m.items():
            if metric in ("n_train", "n_test"):
                continue
            rows.append({
                "celltype": celltype, "modality": modality, "model": model,
                "train_domain": train_domain, "test_domain": test_domain,
                "seed": seed, "metric": metric,
                "stratum": f"condition:{cond}", "value": value,
                "n_train": metrics["n_train"], "n_test": int(mask.sum()),
            })
    for donor in np.unique(donor_test):
        if donor == "synthetic":
            continue
        mask = donor_test == donor
        if mask.sum() < min_stratum or len(np.unique(y_test[mask])) < 2:
            continue
        m = score(y_test[mask], prob[mask], pred[mask],
                  metrics["n_train"], int(mask.sum()))
        for metric, value in m.items():
            if metric in ("n_train", "n_test"):
                continue
            rows.append({
                "celltype": celltype, "modality": modality, "model": model,
                "train_domain": train_domain, "test_domain": test_domain,
                "seed": seed, "metric": metric,
                "stratum": f"donor:{donor}", "value": value,
                "n_train": metrics["n_train"], "n_test": int(mask.sum()),
            })


def condition_control(cond_train, cond_test, y_train, y_test, c_test, d_test,
                      seed):
    xtr = (cond_train == "PD").astype(float)[:, None]
    xte = (cond_test == "PD").astype(float)[:, None]
    res = evaluate(y_train, xtr, y_test, xte, c_test, d_test,
                    "logistic_regression", seed)
    if res is None:
        return None
    return res[0], None, None


def permuted_control(y_train, x_train, y_test, x_test, c_test, d_test, seed):
    rng = np.random.default_rng(seed)
    y_perm = rng.permutation(y_train)
    res = evaluate(y_perm, x_train, y_test, x_test, c_test, d_test,
                    "logistic_regression", seed)
    if res is None:
        return None
    return res[0], None, None


def majority_control(y_train, y_test, c_test, d_test):
    vals, counts = np.unique(y_train, return_counts=True)
    pred = np.full(len(y_test), vals[np.argmax(counts)])
    prob = pred.astype(float)
    metrics = score(y_test, prob, pred, len(y_train), len(y_test))
    return metrics, None, None


def donor_control(x, y, condition, donor, train_mask, test_mask, modality,
                  seed, celltype):
    x_tr, x_te = normalise(x[train_mask], x[test_mask], modality)
    donors_tr = np.unique(donor[train_mask])
    donors_te = np.unique(donor[test_mask])
    if len(donors_tr) < 2 or len(donors_te) < 2:
        return []
    xtr = np.vstack([x_tr[donor[train_mask] == d].mean(axis=0) for d in donors_tr])
    xte = np.vstack([x_te[donor[test_mask] == d].mean(axis=0) for d in donors_te])
    ytr = np.array([
        np.bincount(np.asarray(y[train_mask])[donor[train_mask] == d]).argmax()
        for d in donors_tr])
    yte = np.array([
        np.bincount(np.asarray(y[test_mask])[donor[test_mask] == d]).argmax()
        for d in donors_te])
    res = evaluate(ytr, xtr, yte, xte, None, None, "logistic_regression", seed)
    if res is None:
        return []
    metrics = res[0]
    return [{
        "celltype": celltype, "modality": modality, "model": "donor_only",
        "train_domain": "real", "test_domain": "real", "seed": seed,
        "metric": metric, "stratum": "overall", "value": value,
        "n_train": len(donors_tr), "n_test": len(donors_te),
    } for metric, value in metrics.items() if metric not in ("n_train", "n_test")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/subtype")
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--out", default="results/subtype/transfer")
    ap.add_argument("--celltypes", default="Oligodendrocyte,Microglia")
    ap.add_argument("--modalities", default="rna,atac")
    ap.add_argument("--seeds", default="101,202,303,404,505")
    ap.add_argument("--models", default="logistic_regression,random_forest")
    ap.add_argument("--rf-pairs", default="real->real,syn_alt->real,real->syn_alt,syn_null->real,syn_null_zero->real")
    ap.add_argument("--max-train", type=int, default=10000)
    ap.add_argument("--max-test", type=int, default=10000)
    ap.add_argument("--min-stratum", type=int, default=20)
    args = ap.parse_args()

    root = Path(args.root)
    ingest = Path(args.ingest)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    seeds = [int(s) for s in args.seeds.split(",") if s]
    celltypes = [c.strip() for c in args.celltypes.split(",") if c.strip()]
    modalities = [m.strip() for m in args.modalities.split(",") if m.strip()]
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    rf_pairs = {p.strip() for p in args.rf_pairs.split(",") if p.strip()}

    rows = []
    jaccard_rows = []

    for celltype in celltypes:
        for modality in modalities:
            print(f"=== {celltype}/{modality} ===", flush=True)
            real = load_real(root, ingest, celltype, modality)
            real_train = real["split"] == "train"
            real_test = real["split"] == "test"
            x_real = real["x"].toarray()
            y_real = real["y"]
            c_real = real["condition"]
            d_real = real["donor"]

            # feature overlap with the fixed n24 feature set
            var = pd.read_csv(root / celltype / modality / "features.csv")
            fixed = pd.read_csv(ingest / f"{modality}_features_sel.tsv", header=None)[0]
            vset = set(var["feature"].astype(str))
            fset = set(fixed.astype(str))
            jaccard_rows.append({
                "celltype": celltype, "modality": modality,
                "variable_features": len(vset), "fixed_features": len(fset),
                "intersection": len(vset & fset), "union": len(vset | fset),
                "jaccard": len(vset & fset) / max(len(vset | fset), 1),
            })

            # ---- real -> real, fixed split, one row per model ----------------
            tr = np.flatnonzero(real_train)
            te = np.flatnonzero(real_test)
            tr = stratified_take(tr, y_real, c_real, d_real, args.max_train, 0, True)
            te = stratified_take(te, y_real, c_real, d_real, args.max_test, 1, True)
            xr_tr, xr_te = normalise(x_real[tr], x_real[te], modality)
            for model in models:
                res = evaluate(y_real[tr], xr_tr, y_real[te], xr_te,
                                c_real[te], d_real[te], model, 0)
                if res is None:
                    continue
                add_rows(rows, celltype, modality, model, "real", "real", 0,
                          res, y_real[te], c_real[te], d_real[te],
                          args.min_stratum)
            # Controls on the seed-free real -> real ceiling.
            cond = condition_control(
                c_real[tr], c_real[te], y_real[tr], y_real[te],
                c_real[te], d_real[te], 0)
            add_rows(rows, celltype, modality, "condition_only", "real",
                      "real", 0, cond, y_real[te], c_real[te], d_real[te],
                      args.min_stratum)
            perm = permuted_control(
                y_real[tr], xr_tr, y_real[te], xr_te, c_real[te],
                d_real[te], 0)
            add_rows(rows, celltype, modality, "permuted_labels", "real",
                      "real", 0, perm, y_real[te], c_real[te], d_real[te],
                      args.min_stratum)
            maj = majority_control(y_real[tr], y_real[te], c_real[te],
                                    d_real[te])
            add_rows(rows, celltype, modality, "majority", "real", "real", 0,
                      maj, y_real[te], c_real[te], d_real[te],
                      args.min_stratum)
            # donor-only ceiling
            donor_rows = donor_control(x_real, y_real, c_real, d_real,
                                       real_train, real_test, modality, 0,
                                       celltype)
            rows.extend(donor_rows)

            for seed in seeds:
                syn = {name: load_syn(root, celltype, modality, mode, seed)
                       for name, mode in MODES.items()}
                for name, data in syn.items():
                    if len(np.unique(data["y"])) < 2:
                        print(f"  {name}/{seed}: one subtype in synthetic data", flush=True)

                # main matrix: real and the three synthetic domains
                pairs = [(a, b) for a in DOMAINS for b in DOMAINS]
                for train_name, test_name in pairs:
                    if train_name == "real" and test_name == "real":
                        continue  # handled above, seed-free
                    pair = f"{train_name}->{test_name}"
                    pair_models = models
                    if "random_forest" in pair_models and pair not in rf_pairs:
                        pair_models = [m for m in pair_models if m != "random_forest"]
                    train = real if train_name == "real" else syn[train_name]
                    test = real if test_name == "real" else syn[test_name]

                    if train_name == "real":
                        tr = np.flatnonzero(train["split"] == "train")
                    elif train_name == test_name:
                        # split synthetic cells deterministically for syn -> syn
                        idx = np.arange(len(train["y"]))
                        tr, _ = train_test_split(
                            idx, test_size=0.5, random_state=seed,
                            stratify=np.char.add(
                                np.char.add(train["y"].astype(str), "|"),
                                train["condition"].astype(str)),
                        )
                    else:
                        tr = np.arange(len(train["y"]))
                    tr = stratified_take(
                        tr, train["y"], train["condition"], train["donor"],
                        args.max_train, seed, train_name != "real")

                    if test_name == "real":
                        te = np.flatnonzero(test["split"] == "test")
                    elif train_name == test_name:
                        idx = np.arange(len(test["y"]))
                        _, te = train_test_split(
                            idx, test_size=0.5, random_state=seed,
                            stratify=np.char.add(
                                np.char.add(test["y"].astype(str), "|"),
                                test["condition"].astype(str)),
                        )
                    else:
                        te = np.arange(len(test["y"]))
                    te = stratified_take(
                        te, test["y"], test["condition"], test["donor"],
                        args.max_test, seed + 1, test_name != "real")

                    xtr, xte = normalise(
                        train["x"][tr].toarray(), test["x"][te].toarray(), modality)
                    ytr = train["y"][tr]
                    yte = test["y"][te]
                    cte = test["condition"][te]
                    dte = test["donor"][te]

                    for model in pair_models:
                        res = evaluate(ytr, xtr, yte, xte, cte, dte,
                                        model, seed)
                        if res is None:
                            continue
                        add_rows(rows, celltype, modality, model, train_name,
                                  test_name, seed, res, yte, cte, dte,
                                  args.min_stratum)

                    # controls, logistic only
                    cond = condition_control(
                        train["condition"][tr], test["condition"][te], ytr,
                        yte, cte, dte, seed)
                    add_rows(rows, celltype, modality, "condition_only",
                              train_name, test_name, seed, cond, yte, cte, dte,
                              args.min_stratum)
                    perm = permuted_control(
                        ytr, xtr, yte, xte, cte, dte, seed)
                    add_rows(rows, celltype, modality, "permuted_labels",
                              train_name, test_name, seed, perm, yte, cte, dte,
                              args.min_stratum)
                    maj = majority_control(ytr, yte, cte, dte)
                    add_rows(rows, celltype, modality, "majority",
                              train_name, test_name, seed, maj, yte, cte, dte,
                              args.min_stratum)
                print(f"  seed {seed} done", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(out / "transfer.csv", index=False)
    pd.DataFrame(jaccard_rows).to_csv(out / "jaccard.csv", index=False)

    if table.empty:
        summary_json = {"status": "no_rows"}
    else:
        group = ["celltype", "modality", "model", "train_domain",
                 "test_domain", "metric", "stratum"]
        summary = (table.groupby(group, dropna=False)
                   .agg(value_mean=("value", "mean"), value_sd=("value", "std"),
                        n=("value", "size")).reset_index())
        summary.to_csv(out / "transfer_summary.csv", index=False)
        headline = (table[table["stratum"] == "overall"]
                    .groupby(["celltype", "modality", "model",
                              "train_domain", "test_domain", "metric"])
                    .agg(mean=("value", "mean"), sd=("value", "std"),
                         n=("value", "size")).reset_index())
        summary_json = {
            "status": "ok",
            "rows": int(len(table)),
            "headline": headline.round(4).to_dict(orient="records"),
        }
    (out / "summary.json").write_text(json.dumps(summary_json, indent=2))
    print(f"wrote {out}/transfer.csv ({len(table)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
