#!/usr/bin/env python3
"""Stage 9c - reversal proof of concept.

Validates the counterfactual reversal analysis on data where the truth is known,
and stress-tests it against the donor confound. Three tests, all cheap (effect
sizes + logistic regression, no scDesign3 refits):

1. Known-effect recovery
   The `alt` simulation carries the fitted condition effect; the real holdout cells
   are an independent draw. Rank features by the standardised PD - Aged effect
   (Cohen's d) on real holdout cells and on each `alt` seed. A reversal that
   recovers real biology should rank the same features (Spearman, top-k Jaccard).

2. Null specificity
   The `null` simulation has every condition coefficient zeroed. The same effect
   statistic should collapse toward zero: small SD, few features above the real
   |d| q95, and no correlation with the real ranking.

3. Donor-confound stress test
   Condition and donor are confounded when there is one donor per condition (n02)
   and only weakly separated at larger N. Measure:
     - donor separability (multiclass balanced accuracy, cell-level CV)
     - condition separability with donor-grouped CV (honest: test on unseen donors)
     - the same after removing the donor mean (donor-centred features)
     - donor-level label permutation: reassign whole donors to conditions and
       recompute condition separability; empirical p and confound index
   `identifiable=false` when a condition has a single donor.

Usage:
    python src/reversal_poc.py --ingest results/scale/n24/ingest \\
        --features results/scale/n24/features \\
        --split results/scale/n24/split/split.csv \\
        --obs results/scale/n24/cluster/obs_annotated.csv \\
        --sim results/scale/n24/simulation \\
        --out results/scale/n24/reversal/poc --label n24

Outputs (out/):
    recovery.csv        per modality x celltype x seed effect recovery
    donor_confound.csv  per modality x celltype confound metrics
    summary.json        parameters + aggregated headline numbers
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
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

EPS = 1e-9


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


def log_cpm(mat: sp.csr_matrix) -> sp.csr_matrix:
    """features x cells -> log1p CPM, still sparse."""
    lib = np.asarray(mat.sum(axis=0)).ravel().astype(np.float64)
    lib[lib == 0] = 1.0
    scaled = mat.astype(np.float64).multiply(1e6 / lib[None, :]).tocsr()
    return scaled.log1p()


def load_real(ingest: Path, modality: str, feats: list[str]) -> tuple[sp.csr_matrix, list[str]]:
    """Full real matrix for the selected features: features x cells."""
    sel = ingest / f"{modality}_counts_sel.mtx.gz"
    if sel.exists():
        mat = read_mtx(sel)
        names = load_features(ingest / f"{modality}_features_sel.tsv")
    else:
        mat = read_mtx(find(ingest / f"{modality}_counts.mtx"))
        names = load_features(ingest / f"{modality}_features.tsv")
    pos = {f: i for i, f in enumerate(names)}
    missing = [f for f in feats if f not in pos]
    if missing:
        raise SystemExit(f"{modality}: {len(missing)} features missing (first {missing[0]})")
    rows = np.array([pos[f] for f in feats])
    obs = pd.read_csv(ingest / "obs.csv", dtype={"cell": str})
    return log_cpm(mat[rows]), obs["cell"].tolist()


def load_sim(sim: Path, modality: str, mode: str, seed: int,
             feats: list[str]) -> tuple[sp.csr_matrix, pd.DataFrame]:
    stem = sim / modality / mode / f"seed{seed}_counts.mtx"
    mat = read_mtx(find(stem))
    if mat.shape[0] != len(feats):
        raise SystemExit(
            f"{stem}: {mat.shape[0]} rows but {len(feats)} selected features")
    obs = pd.read_csv(sim / modality / mode / "obs.csv", dtype={"cell": str})
    return log_cpm(mat), obs


# ------------------------------------------------------------------ effects --
def cohens_d(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """x, y: features x cells. Standardised mean difference x - y."""
    mx, my = x.mean(axis=1), y.mean(axis=1)
    vx, vy = x.var(axis=1, ddof=1), y.var(axis=1, ddof=1)
    pooled = np.sqrt((vx + vy) / 2.0)
    d = np.full_like(mx, np.nan)
    ok = pooled >= EPS
    d[ok] = (mx[ok] - my[ok]) / pooled[ok]
    return d


def effect_on(mat: sp.csr_matrix, celltypes: np.ndarray, groups: np.ndarray,
              celltype: str, target: str, reference: str) -> tuple[np.ndarray, int, int]:
    in_ct = celltypes == celltype
    tgt = in_ct & (groups == target)
    ref = in_ct & (groups == reference)
    if tgt.sum() < 2 or ref.sum() < 2:
        return np.full(mat.shape[0], np.nan), int(tgt.sum()), int(ref.sum())
    dense = mat[:, tgt | ref].toarray()
    tgt_mask = tgt[tgt | ref]
    return cohens_d(dense[:, tgt_mask], dense[:, ~tgt_mask]), int(tgt.sum()), int(ref.sum())


def rank_agreement(a: np.ndarray, b: np.ndarray, top_n: int) -> tuple[float, float]:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan"), float("nan")
    a, b = a[ok], b[ok]
    rho = float(pd.Series(a).corr(pd.Series(b), method="spearman"))
    k = min(top_n, len(a))
    ia = set(np.argsort(-np.abs(a))[:k])
    ib = set(np.argsort(-np.abs(b))[:k])
    return rho, len(ia & ib) / k


# ------------------------------------------------------------ donor confound --
def make_clf():
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=5000, C=1.0, class_weight="balanced"),
    )


def cv_balanced_accuracy(x: np.ndarray, y: np.ndarray, groups: np.ndarray | None,
                         seed: int = 0) -> float:
    """Balanced accuracy. Grouped CV when `groups` is given (unseen donors).

    With grouped CV the split count must be bounded by the number of donors per
    class, otherwise a training fold can contain a single class.
    """
    if groups is None:
        counts = np.bincount(y)
        n_splits = int(min(5, counts[counts > 0].min()))
        if n_splits < 2:
            return float("nan")
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        try:
            pred = cross_val_predict(make_clf(), x, y, cv=cv)
        except ValueError:
            return float("nan")
    else:
        per_class = [len(np.unique(groups[y == c])) for c in np.unique(y)]
        n_splits = int(min(5, min(per_class)))
        if n_splits < 2:
            return float("nan")
        cv = GroupKFold(n_splits=n_splits)
        try:
            pred = cross_val_predict(make_clf(), x, y, cv=cv, groups=groups)
        except ValueError:
            return float("nan")
    return float(balanced_accuracy_score(y, pred))


def donor_permutation_test(x: np.ndarray, y: np.ndarray, donors: np.ndarray,
                            n_perm: int, seed: int) -> np.ndarray:
    """Reassign whole donors to conditions, then recompute condition separability."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(donors)
    donor_label = {d: y[donors == d][0] for d in uniq}
    labels = np.array([donor_label[d] for d in uniq])
    null = []
    for _ in range(n_perm):
        perm = rng.permutation(labels)
        mapping = dict(zip(uniq, perm))
        y_perm = np.array([mapping[d] for d in donors])
        if len(np.unique(y_perm)) < 2:
            continue
        null.append(cv_balanced_accuracy(x, y_perm, donors, seed=seed))
    return np.array([v for v in null if np.isfinite(v)])


# ---------------------------------------------------------------------- main --
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--split", default="results/split/split.csv")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--sim", default="results/simulation")
    ap.add_argument("--out", default="results/reversal/poc")
    ap.add_argument("--label", default="")
    ap.add_argument("--seeds", default="101,202,303,404,505")
    ap.add_argument("--target", default="PD")
    ap.add_argument("--reference", default="Aged")
    ap.add_argument("--min-cells", type=int, default=30)
    ap.add_argument("--top-n", type=int, default=100)
    ap.add_argument("--permutations", type=int, default=10)
    ap.add_argument("--seed", type=int, default=20261110)
    args = ap.parse_args()

    ingest = Path(args.ingest)
    sim = Path(args.sim)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    seeds = [int(s) for s in args.seeds.split(",") if s]

    obs = pd.read_csv(args.obs, dtype={"cell": str}).set_index("cell")
    split = pd.read_csv(args.split, dtype={"cell": str})
    fit = set(split.loc[split["split"] == "fit", "cell"])
    cell_order = obs.index.tolist()
    cpos = {c: i for i, c in enumerate(cell_order)}

    recovery_rows, confound_rows = [], []

    for modality in ("rna", "atac"):
        feats = load_features(Path(args.features) / f"{modality}_features.tsv")
        if not feats:
            continue
        real_all, real_cells = load_real(ingest, modality, feats)
        if real_cells != cell_order:
            raise SystemExit(f"{modality}: ingest cell order differs from obs")

        first_obs = pd.read_csv(sim / modality / "alt" / "obs.csv", dtype={"cell": str})
        sim_cells = first_obs["cell"].tolist()
        hold_cols = np.array([cpos[c] for c in sim_cells])
        real_hold = real_all[:, hold_cols]
        hold_ct = obs.loc[sim_cells, "celltype"].to_numpy()
        hold_grp = obs.loc[sim_cells, "group"].to_numpy()

        alt_blocks = {s: load_sim(sim, modality, "alt", s, feats) for s in seeds}
        null_blocks = {s: load_sim(sim, modality, "null", s, feats) for s in seeds}

        for celltype in pd.unique(hold_ct):
            in_ct = hold_ct == celltype
            n_tgt = int((in_ct & (hold_grp == args.target)).sum())
            n_ref = int((in_ct & (hold_grp == args.reference)).sum())
            if min(n_tgt, n_ref) < args.min_cells:
                continue

            d_real, _, _ = effect_on(real_hold, hold_ct, hold_grp, celltype,
                                     args.target, args.reference)
            if not np.isfinite(d_real).any():
                continue
            q95 = float(np.nanquantile(np.abs(d_real), 0.95))

            for s in seeds:
                alt_mat, alt_obs = alt_blocks[s]
                null_mat, null_obs = null_blocks[s]
                d_alt, n_alt_t, n_alt_r = effect_on(
                    alt_mat, alt_obs["celltype"].to_numpy(),
                    alt_obs["condition"].to_numpy(), celltype,
                    args.target, args.reference)
                d_null, _, _ = effect_on(
                    null_mat, null_obs["celltype"].to_numpy(),
                    null_obs["condition"].to_numpy(), celltype,
                    args.target, args.reference)

                rho_alt, jac_alt = rank_agreement(d_real, d_alt, args.top_n)
                rho_null, jac_null = rank_agreement(d_real, d_null, args.top_n)
                sd_alt = float(np.nanstd(d_alt))
                sd_null = float(np.nanstd(d_null))

                recovery_rows.append({
                    "label": args.label, "modality": modality,
                    "celltype": celltype, "seed": s,
                    "n_real_target": n_tgt, "n_real_reference": n_ref,
                    "n_sim_target": n_alt_t, "n_sim_reference": n_alt_r,
                    "spearman_real_alt": rho_alt,
                    "topn_jaccard_alt": jac_alt,
                    "spearman_real_null": rho_null,
                    "topn_jaccard_null": jac_null,
                    "sd_real": float(np.nanstd(d_real)),
                    "sd_alt": sd_alt, "sd_null": sd_null,
                    "sd_ratio_alt_null": (sd_alt / sd_null if sd_null > EPS else float("nan")),
                    "real_abs_q95": q95,
                    "frac_null_above_real_q95": float(np.mean(np.abs(d_null) > q95)),
                })

            # ---- donor confound on the real fit cells ------------------------
            fit_cells = [c for c in cell_order
                         if c in fit and obs.at[c, "celltype"] == celltype]
            fit_cols = np.array([cpos[c] for c in fit_cells])
            x_fit = real_all[:, fit_cols].toarray().T  # cells x features
            y_fit = obs.loc[fit_cells, "group"].to_numpy()
            donor = obs.loc[fit_cells, "donor"].to_numpy()

            n_donors = len(np.unique(donor))
            n_tgt_donors = len(np.unique(donor[y_fit == args.target]))
            n_ref_donors = len(np.unique(donor[y_fit == args.reference]))
            identifiable = min(n_tgt_donors, n_ref_donors) >= 2

            donor_bal = cv_balanced_accuracy(x_fit, pd.factorize(donor)[0], None, args.seed)
            cond_bal = cv_balanced_accuracy(x_fit, (y_fit == args.target).astype(int),
                                            donor, args.seed)

            uniq_d = np.unique(donor)
            means = {d: x_fit[donor == d].mean(axis=0) for d in uniq_d}
            centred = x_fit - np.vstack([means[d] for d in donor])
            cond_bal_centred = cv_balanced_accuracy(
                centred, (y_fit == args.target).astype(int), donor, args.seed)

            if identifiable:
                null = donor_permutation_test(
                    x_fit, (y_fit == args.target).astype(int), donor,
                    args.permutations, args.seed)
                perm_mean = float(null.mean()) if len(null) else float("nan")
                perm_q95 = float(np.quantile(null, 0.95)) if len(null) else float("nan")
                p = float((1 + (null >= cond_bal).sum()) / (len(null) + 1)) if len(null) else float("nan")
                confound_index = (float((cond_bal - perm_mean) / (1 - perm_mean))
                                   if len(null) and perm_mean < 1 else float("nan"))
            else:
                perm_mean = perm_q95 = p = confound_index = float("nan")

            confound_rows.append({
                "label": args.label, "modality": modality, "celltype": celltype,
                "n_cells": int(len(fit_cells)), "n_donors": n_donors,
                "n_target_donors": n_tgt_donors,
                "n_reference_donors": n_ref_donors,
                "donor_balanced_accuracy": donor_bal,
                "condition_balanced_accuracy": cond_bal,
                "condition_balanced_accuracy_donor_centred": cond_bal_centred,
                "perm_mean": perm_mean, "perm_q95": perm_q95,
                "empirical_p": p, "confound_index": confound_index,
                "identifiable": identifiable,
            })

    recovery = pd.DataFrame(recovery_rows)
    confound = pd.DataFrame(confound_rows)
    recovery.to_csv(out / "recovery.csv", index=False)
    confound.to_csv(out / "donor_confound.csv", index=False)

    def mean(df, col):
        return float(df[col].mean()) if len(df) else float("nan")

    payload = {
        "label": args.label,
        "target": args.target, "reference": args.reference,
        "seeds": seeds, "top_n": args.top_n,
        "permutations": args.permutations,
        "recovery": {
            "spearman_real_alt": mean(recovery, "spearman_real_alt"),
            "topn_jaccard_alt": mean(recovery, "topn_jaccard_alt"),
            "spearman_real_null": mean(recovery, "spearman_real_null"),
            "frac_null_above_real_q95": mean(recovery, "frac_null_above_real_q95"),
            "sd_ratio_alt_null": mean(recovery, "sd_ratio_alt_null"),
        },
        "donor_confound": {
            "donor_balanced_accuracy": mean(confound, "donor_balanced_accuracy"),
            "condition_balanced_accuracy": mean(confound, "condition_balanced_accuracy"),
            "confound_index": mean(confound, "confound_index"),
            "identifiable_strata": int(confound["identifiable"].sum()) if len(confound) else 0,
            "strata": int(len(confound)),
        },
        "strata": recovery_rows,
        "donor_confound_strata": confound_rows,
    }
    (out / "summary.json").write_text(json.dumps(payload, indent=2, default=float))
    print(f"wrote {out}/recovery.csv, donor_confound.csv, summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
