#!/usr/bin/env python3
"""Stage 8 - null / alternative hypothesis benchmarking.

The point of the simulator is not to look pretty but to be usable as a test bed:
a method run on `null` data must not call anything, and on `alt` data it must
recover the effects that are really there.

Two contrasts are benchmarked:

  condition  within one cell type, Aged vs PD
  celltype   one cell type vs the rest, within one condition

`null` zeroes only the `condition` coefficients of the scDesign3 fit, so the
cell-type structure survives. That gives a built-in asymmetry worth checking:
condition contrasts must be null under `null`, while celltype contrasts must
still have power under `null`.

Ground truth is the same test run on the real fit cells. A feature is "true
positive" for a contrast when it is significant in the real data at
`--truth-fdr` and has |log2FC| >= `--truth-lfc`.

Usage:
    python src/hypothesis.py [--features results/features]
                             [--split results/split/split.csv]
                             [--obs results/cluster/obs_annotated.csv]
                             [--sim results/simulation]
                             [--out results/hypothesis]

Outputs (results/hypothesis/):
    tests.csv        one row per modality x mode x seed x contrast
    calibration.csv  null calibration aggregated over seeds
    power.csv        alt power aggregated over seeds
    truth.csv        the real-data reference calls per contrast
    summary.json
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
from scipy.stats import mannwhitneyu, ttest_ind

MIN_CELLS = 20


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


def lognorm(mat: sp.csr_matrix) -> np.ndarray:
    """cells x features, log1p CPM."""
    x = mat.T.tocsr().astype(np.float64)
    lib = np.asarray(x.sum(axis=1)).ravel()
    lib[lib == 0] = 1.0
    return np.log1p(x.multiply(1e6 / lib[:, None]).toarray())


def bh(p: np.ndarray) -> np.ndarray:
    n = p.size
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def de_test(x: np.ndarray, mask_a: np.ndarray, mask_b: np.ndarray, test: str) -> dict:
    """x is cells x features; returns log2FC (a - b), p, q."""
    a, b = x[mask_a], x[mask_b]
    if test == "wilcox":
        # mannwhitneyu is vectorised over the feature axis
        res = mannwhitneyu(a, b, axis=0, alternative="two-sided")
        p = res.pvalue
    elif test == "welch":
        p = ttest_ind(a, b, axis=0, equal_var=False).pvalue
    else:
        raise SystemExit(f"unknown test: {test}")
    p = np.nan_to_num(p, nan=1.0)
    # log2 fold change of the mean normalised expression, with a small offset so
    # zero means do not produce -inf
    fc = np.log2((a.mean(axis=0) + 1e-6) / (b.mean(axis=0) + 1e-6))
    return {"log2fc": fc, "p": p, "q": bh(p)}


# ---------------------------------------------------------------- contrasts --
def contrasts(obs: pd.DataFrame, min_cells: int = MIN_CELLS) -> list[dict]:
    """Build the contrast list from the celltype x condition table."""
    out = []
    tab = pd.crosstab(obs["celltype"], obs["condition"])
    for ct in tab.index:
        row = tab.loc[ct]
        if len(row) >= 2 and row.min() >= min_cells:
            a, b = row.index[0], row.index[1]
            out.append({
                "type": "condition", "name": f"{ct}: {a} vs {b}",
                "celltype": ct, "a": a, "b": b,
                "mask_a": (obs["celltype"] == ct) & (obs["condition"] == a),
                "mask_b": (obs["celltype"] == ct) & (obs["condition"] == b),
            })
    for cond in tab.columns:
        col = tab[cond]
        present = col[col >= min_cells].index
        for ct in present:
            out.append({
                "type": "celltype", "name": f"{cond}: {ct} vs rest",
                "celltype": ct, "a": ct, "b": "rest",
                "mask_a": (obs["celltype"] == ct) & (obs["condition"] == cond),
                "mask_b": (obs["celltype"] != ct) & (obs["condition"] == cond),
            })
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", default="results/ingest")
    ap.add_argument("--features", default="results/features")
    ap.add_argument("--split", default="results/split/split.csv")
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--sim", default="results/simulation")
    ap.add_argument("--out", default="results/hypothesis")
    ap.add_argument("--modality", default="rna,atac")
    ap.add_argument("--seeds", default="")
    ap.add_argument("--test", default="wilcox", choices=["wilcox", "welch"])
    ap.add_argument("--truth-fdr", type=float, default=0.05)
    ap.add_argument("--truth-lfc", type=float, default=0.5)
    ap.add_argument("--fdr-grid", default="0.01,0.05,0.1,0.2")
    ap.add_argument("--min-cells", type=int, default=MIN_CELLS)
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    fdr_grid = [float(v) for v in args.fdr_grid.split(",")]

    obs = pd.read_csv(args.obs, dtype={"cell": str})
    # the real obs calls it `group`; simulated obs.csv calls it `condition`
    if "condition" not in obs.columns and "group" in obs.columns:
        obs = obs.rename(columns={"group": "condition"})
    split = pd.read_csv(args.split, dtype={"cell": str})
    manifest = pd.read_csv(Path(args.sim) / "manifest.csv", keep_default_na=False)
    manifest["seed"] = manifest["seed"].astype(int)
    if args.seeds:
        want = set(args.seeds.split(","))
        manifest = manifest[manifest["seed"].astype(str).isin(want)]

    modalities = [m for m in args.modality.split(",") if m]
    rows, truth_rows = [], []

    for modality in modalities:
        feats = [l for l in (Path(args.features) / f"{modality}_features.tsv")
                 .read_text().splitlines() if l]
        print(f"[{modality}] {len(feats)} features", flush=True)

        # ---- real reference calls on the fit cells ------------------------
        fit_cells = split.loc[split["split"] == "fit", "cell"].tolist()
        fit_obs = obs.set_index("cell").loc[fit_cells].reset_index()
        real = read_mtx(find(Path(args.ingest) / f"{modality}_counts.mtx"))
        all_cells = obs["cell"].tolist()
        feat_all = [l for l in (Path(args.ingest) / f"{modality}_features.tsv")
                    .read_text().splitlines() if l]
        fpos = {f: i for i, f in enumerate(feat_all)}
        cpos = {c: i for i, c in enumerate(all_cells)}
        rows_idx = np.array([fpos[f] for f in feats])
        cols_idx = np.array([cpos[c] for c in fit_cells])
        real_x = lognorm(real[rows_idx][:, cols_idx])

        truth: dict[tuple, set] = {}
        fc_ref: dict[tuple, np.ndarray] = {}
        for spec in contrasts(fit_obs, args.min_cells):
            res = de_test(real_x, spec["mask_a"].to_numpy(), spec["mask_b"].to_numpy(),
                          args.test)
            sig = (res["q"] <= args.truth_fdr) & (np.abs(res["log2fc"]) >= args.truth_lfc)
            key = (spec["type"], spec["name"])
            truth[key] = set(np.flatnonzero(sig))
            fc_ref[key] = res["log2fc"]
            truth_rows.append({
                "modality": modality, "contrast_type": spec["type"],
                "contrast": spec["name"], "n_truth": int(sig.sum()),
                "n_a": int(spec["mask_a"].sum()), "n_b": int(spec["mask_b"].sum()),
            })
        print(f"[{modality}] {len(truth)} contrasts, "
              f"{sum(len(v) for v in truth.values())} truth features", flush=True)

        # ---- simulated runs ----------------------------------------------
        for _, mrow in manifest[manifest["modality"] == modality].iterrows():
            mode, seed = mrow["mode"], int(mrow["seed"])
            sim = read_mtx(Path(mrow["path"]))
            if sim.shape[0] != len(feats):
                raise SystemExit(f"{mrow['path']}: {sim.shape[0]} rows != {len(feats)}")
            sim_obs = pd.read_csv(Path(mrow["path"]).parent / "obs.csv",
                                  dtype={"cell": str})
            sim_x = lognorm(sim)
            for spec in contrasts(sim_obs, args.min_cells):
                key = (spec["type"], spec["name"])
                truth_set = truth.get(key, set())
                res = de_test(sim_x, spec["mask_a"].to_numpy(),
                              spec["mask_b"].to_numpy(), args.test)
                rec = {
                    "modality": modality, "mode": mode, "seed": seed,
                    "contrast_type": spec["type"], "contrast": spec["name"],
                    "n_a": int(spec["mask_a"].sum()), "n_b": int(spec["mask_b"].sum()),
                    "n_truth": len(truth_set),
                }
                for a in fdr_grid:
                    sig = res["q"] <= a
                    rec[f"n_sig_{a}"] = int(sig.sum())
                    if truth_set:
                        tp = len(set(np.flatnonzero(sig)) & truth_set)
                        rec[f"power_{a}"] = tp / len(truth_set)
                        rec[f"fdp_{a}"] = 1 - tp / max(int(sig.sum()), 1)
                    else:
                        rec[f"power_{a}"] = np.nan
                        rec[f"fdp_{a}"] = np.nan
                # effect-size agreement with the real contrast
                fc_real = fc_ref.get(key)
                if fc_real is not None and fc_real.size == res["log2fc"].size:
                    rec["fc_corr"] = float(np.corrcoef(fc_real, res["log2fc"])[0, 1])
                else:
                    rec["fc_corr"] = np.nan
                rows.append(rec)
            print(f"  {mode} seed {seed}: {len(contrasts(sim_obs, args.min_cells))} contrasts", flush=True)

    tests = pd.DataFrame(rows)
    tests.to_csv(out / "tests.csv", index=False)
    pd.DataFrame(truth_rows).to_csv(out / "truth.csv", index=False)

    def agg(df: pd.DataFrame, label: str) -> pd.DataFrame:
        if df.empty:
            return pd.DataFrame()
        cols = [c for c in df.columns if c.startswith(("n_sig_", "power_", "fdp_"))]
        cols += ["fc_corr"]
        g = df.groupby(["modality", "contrast_type", "contrast"])[cols]
        res = g.mean().reset_index()
        res.insert(0, "set", label)
        return res

    calib = agg(tests[tests["mode"] == "null"], "null")
    power = agg(tests[tests["mode"] == "alt"], "alt")
    pd.concat([calib, power], ignore_index=True).to_csv(out / "calibration.csv",
                                                        index=False)

    summary = {
        "test": args.test,
        "truth_fdr": args.truth_fdr,
        "truth_lfc": args.truth_lfc,
        "fdr_grid": fdr_grid,
        "seeds": sorted(tests["seed"].unique().tolist()) if not tests.empty else [],
        "modalities": modalities,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {out}/tests.csv, calibration.csv, truth.csv, summary.json")
    if not tests.empty:
        print(tests.groupby(["modality", "mode", "contrast_type"])
              [["n_sig_0.05", "power_0.05"]].mean().round(3).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
