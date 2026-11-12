#!/usr/bin/env python3
"""Merge the scaling steps into one comparison table.

Every step writes its own `results/scale/nXX/` tree. This script reads those
trees and emits one long table plus a compact summary:

    comparison.csv   N, modality, analysis, metric, celltype, label,
                     synthetic_size, seed, value
    summary.json     headline curves and the stability-vs-n24 numbers

The long format is deliberate: sample size and simulation iteration are separate
axes, so a metric can be plotted against either without reshaping.

Usage:
    python src/compare_scaling.py [--root results/scale]
                                 [--steps 2,8,16,24]
                                 [--fixed-step 24]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(path, keep_default_na=False)
    except Exception:
        return None


def read_json(path: Path):
    if not path.exists() or path.stat().st_size == 0:
        return None
    try:
        return json.loads(path.read_text())
    except Exception:
        return None


def read_features(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [l for l in path.read_text().splitlines() if l]


def step_features(step: Path, modality: str) -> list[str]:
    """The feature list this step selected for itself.

    Steps with `features.fixed_from` copy the n24 list into `features/` for the
    simulation, but still write their own selection under `features/selected/`.
    Feature stability must compare the step's own selection, not the copy.
    """
    own = step / "features" / "selected" / f"{modality}_features.tsv"
    if own.exists():
        return read_features(own)
    return read_features(step / "features" / f"{modality}_features.tsv")


def jaccard(a: list[str], b: list[str]) -> float:
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return float("nan")
    return len(sa & sb) / len(sa | sb)


def add(rows: list[dict], n: int, modality: str, analysis: str, metric: str,
        value, celltype: str = "", label: str = "", synthetic_size="",
        seed="") -> None:
    if value is None:
        return
    try:
        v = float(value)
    except (TypeError, ValueError):
        return
    if not np.isfinite(v):
        return
    rows.append({
        "N": n, "modality": modality, "analysis": analysis, "metric": metric,
        "celltype": celltype, "label": label,
        "synthetic_size": synthetic_size, "seed": seed, "value": v,
    })


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/scale")
    ap.add_argument("--steps", default="2,8,16,24")
    ap.add_argument("--fixed-step", type=int, default=24)
    args = ap.parse_args()

    root = Path(args.root)
    steps = [int(s) for s in args.steps.split(",") if s]
    fixed_dir = root / f"n{args.fixed_step:02d}" / "features"
    fixed_rna = read_features(fixed_dir / "rna_features.tsv")
    fixed_atac = read_features(fixed_dir / "atac_features.tsv")

    rows: list[dict] = []
    per_step: dict[int, dict] = {}
    de_by_step: dict[int, pd.DataFrame] = {}

    for n in steps:
        step = root / f"n{n:02d}"
        if not step.exists():
            print(f"skip n{n:02d}: missing {step}")
            continue
        per_step[n] = {}

        obs = read_csv(step / "ingest" / "obs.csv")
        annotated = read_csv(step / "cluster" / "obs_annotated.csv")
        qc = read_csv(step / "ingest" / "qc_per_sample.csv")
        clusters = read_csv(step / "cluster" / "clusters.csv")

        n_cells = len(obs) if obs is not None else None
        n_donors = obs["donor"].nunique() if obs is not None else None
        n_celltypes = annotated["celltype"].nunique() if annotated is not None else None
        n_clusters = clusters["cluster"].nunique() if clusters is not None else None
        doublet_rate = None
        if qc is not None and "doublets" in qc.columns and "cells_in" in qc.columns:
            doublet_rate = float(qc["doublets"].sum() / max(qc["cells_in"].sum(), 1))

        for metric, value in (
            ("n_cells", n_cells), ("n_donors", n_donors),
            ("n_celltypes", n_celltypes), ("n_clusters", n_clusters),
            ("doublet_rate", doublet_rate),
        ):
            add(rows, n, "both", "qc", metric, value)
        per_step[n].update({"n_cells": n_cells, "n_donors": n_donors,
                            "n_celltypes": n_celltypes, "n_clusters": n_clusters})

        # ---- feature stability --------------------------------------------
        step_rna = step_features(step, "rna")
        step_atac = step_features(step, "atac")
        add(rows, n, "rna", "features", "jaccard_vs_fixed", jaccard(step_rna, fixed_rna))
        add(rows, n, "atac", "features", "jaccard_vs_fixed", jaccard(step_atac, fixed_atac))
        add(rows, n, "rna", "features", "n_selected", len(step_rna))
        add(rows, n, "atac", "features", "n_selected", len(step_atac))

        # ---- fidelity ------------------------------------------------------
        disc = read_csv(step / "fidelity" / "discriminator.csv")
        if disc is not None:
            sub = disc[disc["comparison"] == "real_vs_synthetic"]
            for modality, g in sub.groupby("modality"):
                for model, h in g.groupby("model"):
                    add(rows, n, modality, "fidelity", f"auc_{model}", h["auc"].mean())
                    add(rows, n, modality, "fidelity", f"auc_sd_{model}", h["auc"].std())
                add(rows, n, modality, "fidelity", "auc_mean",
                    g.groupby("seed")["auc"].mean().mean())
                # iteration axis: one AUC per simulation seed
                for seed, h in g.groupby("seed"):
                    add(rows, n, modality, "fidelity", "auc_by_seed",
                        h["auc"].mean(), seed=str(seed))

        # ---- hypothesis ----------------------------------------------------
        tests = read_csv(step / "hypothesis" / "tests.csv")
        if tests is not None:
            cond = tests[tests["contrast_type"] == "condition"]
            for _, r in cond.iterrows():
                for metric in ("power_0.05", "fdp_0.05", "fc_corr"):
                    add(rows, n, r["modality"], "hypothesis_iter", metric,
                        r.get(metric), celltype=str(r.get("contrast", "")),
                        seed=str(r.get("seed", "")))

        calib = read_csv(step / "hypothesis" / "calibration.csv")
        if calib is not None:
            cond = calib[calib["contrast_type"] == "condition"]
            for modality, g in cond.groupby("modality"):
                for set_name in ("alt", "null"):
                    h = g[g["set"] == set_name]
                    if h.empty:
                        continue

                    def col_mean(col: str) -> float:
                        # Empty cells are kept as "" by read_csv; coerce before averaging.
                        return float(pd.to_numeric(h[col], errors="coerce").mean())

                    for col in h.columns:
                        if col.startswith("power_"):
                            fdr = col.split("_")[1]
                            add(rows, n, modality, "hypothesis",
                                f"{set_name}_power_fdr{fdr}", col_mean(col))
                        if col.startswith("fdp_"):
                            fdr = col.split("_")[1]
                            add(rows, n, modality, "hypothesis",
                                f"{set_name}_fdp_fdr{fdr}", col_mean(col))
                        if col.startswith("n_sig_"):
                            fdr = col.split("_")[1]
                            add(rows, n, modality, "hypothesis",
                                f"{set_name}_n_sig_fdr{fdr}", col_mean(col))
                    add(rows, n, modality, "hypothesis", f"{set_name}_fc_corr",
                        col_mean("fc_corr"))

        # ---- reversal ------------------------------------------------------
        rev_model = read_json(step / "reversal" / "model_summary.json")
        if rev_model:
            for s in rev_model.get("strata", []):
                for metric in ("reverse_correlation", "observed_sd", "permuted_sd",
                               "permuted_q95", "frac_observed_above_perm_q95"):
                    add(rows, n, s.get("modality", ""), "reversal", metric,
                        s.get(metric), celltype=s.get("celltype", ""))
        rev = read_json(step / "reversal" / "summary.json")
        if rev:
            for s in rev.get("strata", []):
                for metric in ("spearman_model_vs_classifier",
                               "spearman_reverse_direction", "spearman_permuted",
                               "cv_balanced_accuracy", "top100_overlap"):
                    add(rows, n, s.get("modality", ""), "reversal", metric,
                        s.get(metric), celltype=s.get("celltype", ""))

        # ---- reversal PoC --------------------------------------------------
        poc = read_json(step / "reversal" / "poc" / "summary.json")
        if poc:
            rec = poc.get("recovery", {})
            for metric in ("spearman_real_alt", "topn_jaccard_alt",
                           "spearman_real_null", "frac_null_above_real_q95",
                           "sd_ratio_alt_null"):
                add(rows, n, "", "reversal_poc", metric, rec.get(metric))
            con = poc.get("donor_confound", {})
            for metric in ("donor_balanced_accuracy",
                           "condition_balanced_accuracy", "confound_index"):
                add(rows, n, "", "reversal_poc", metric, con.get(metric))
            add(rows, n, "", "reversal_poc", "identifiable_strata",
                con.get("identifiable_strata"))

        # ---- impact --------------------------------------------------------
        impact = read_csv(step / "impact" / "celltype_impact.csv")
        if impact is not None:
            per_step[n]["impact"] = impact
            for _, r in impact.iterrows():
                for metric in ("de_burden", "mean_abs_lfc_sig", "max_abs_lfc_sig",
                               "donor_auc", "donor_balanced_accuracy",
                               "eta2_condition", "eta2_donor",
                               "condition_variance_share", "pd_fraction",
                               "reversal_mean_abs_lfc",
                               "reversal_frac_above_perm_q95"):
                    add(rows, n, r["modality"], "impact", metric, r.get(metric),
                        celltype=r["celltype"])
            for modality, g in impact.groupby("modality"):
                w = pd.to_numeric(g["n_cells"], errors="coerce").to_numpy(dtype=float)
                if not np.isfinite(w).all() or w.sum() <= 0:
                    continue
                # read_csv keeps empty cells as "", so coerce before averaging.
                auc = pd.to_numeric(g["donor_auc"], errors="coerce")
                burden = pd.to_numeric(g["de_burden"], errors="coerce")
                ok_auc = auc.notna().to_numpy()
                if ok_auc.any():
                    add(rows, n, modality, "impact", "donor_auc_weighted",
                        np.average(auc.to_numpy(dtype=float)[ok_auc],
                                   weights=w[ok_auc]))
                ok = burden.notna().to_numpy()
                if ok.any():
                    add(rows, n, modality, "impact", "de_burden_weighted",
                        np.average(burden.to_numpy(dtype=float)[ok],
                                   weights=w[ok]))

        # ---- utility -------------------------------------------------------
        curve = read_csv(step / "utility" / "learning_curve.csv")
        if curve is not None and not curve.empty:
            per_step[n]["utility"] = curve
            for _, r in curve.iterrows():
                for metric in ("auc_mean", "auc_sd", "bal_mean", "bal_sd"):
                    add(rows, n, r["modality"], "utility", metric, r.get(metric),
                        label=str(r.get("method", "")),
                        synthetic_size=int(r["synthetic_size"]))
        # iteration axis: raw per-seed utility rows
        util = read_csv(step / "utility" / "utility.csv")
        if util is not None and not util.empty:
            syn = util[util["method"] == "synthetic"]
            for _, r in syn.iterrows():
                add(rows, n, r["modality"], "utility_iter", "auc", r.get("auc"),
                    label=str(r.get("method", "")),
                    synthetic_size=r.get("synthetic_size", ""),
                    seed=str(r.get("seed", "")))

        # ---- DE table for stability ---------------------------------------
        de = read_csv(step / "impact" / "de_features.csv")
        if de is not None:
            de_by_step[n] = de

    # ---- stability vs the fixed step --------------------------------------
    ref_de = de_by_step.get(args.fixed_step)
    if ref_de is not None:
        ref = ref_de[["modality", "celltype", "feature", "log2fc"]].rename(
            columns={"log2fc": "log2fc_ref"})
        for n, de in de_by_step.items():
            if n == args.fixed_step:
                continue
            merged = de.merge(ref, on=["modality", "celltype", "feature"])
            for (modality, celltype), g in merged.groupby(["modality", "celltype"]):
                if len(g) < 10:
                    continue
                rho = spearmanr(g["log2fc"], g["log2fc_ref"]).statistic
                add(rows, n, modality, "stability", "spearman_vs_fixed", rho,
                    celltype=celltype)
                k = min(100, len(g))
                a = set(g.reindex(g["log2fc"].abs().sort_values(
                    ascending=False).index)["feature"].head(k))
                b = set(g.reindex(g["log2fc_ref"].abs().sort_values(
                    ascending=False).index)["feature"].head(k))
                add(rows, n, modality, "stability", "jaccard_top100_vs_fixed",
                    len(a & b) / max(k, 1), celltype=celltype)

    table = pd.DataFrame(rows)
    table.to_csv(root / "comparison.csv", index=False)

    # ---- compact summary -------------------------------------------------
    summary = {"steps": steps, "fixed_step": args.fixed_step, "headline": {}}
    if not table.empty:
        for analysis in sorted(table["analysis"].unique()):
            sub = table[table["analysis"] == analysis]
            summary["headline"][analysis] = (
                sub.groupby(["N", "modality", "metric"])["value"]
                .mean().reset_index().round(4).to_dict(orient="records")
            )
    (root / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {root}/comparison.csv and summary.json "
          f"({len(table)} rows, {len(steps)} steps)")
    if not table.empty:
        print(table.groupby(["N", "analysis"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
