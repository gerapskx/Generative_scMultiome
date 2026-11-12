#!/usr/bin/env python3
"""Merge the subtype-transfer experiment into one comparison table.

The transfer stage writes `subtype/transfer/transfer.csv` with one row per
cell type x modality x model x train/test domain x seed x metric x stratum. This
script adds the pieces that are not part of that matrix:

    qc            silhouette, depth/doublet/condition/donor controls, subtype sizes
    concordance   RNA vs ATAC adjusted Rand index on shared cells
    jaccard       overlap of the subtype variable features with the fixed n24 set
    separability  synthetic alt/null CV AUC from the transfer matrix

Outputs `subtype/comparison.csv` (long format) and `subtype/summary.json`
(headline means with the same column names as the scaling comparison).

Usage:
    python src/subtype_compare.py --root results/scale/n24/subtype
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score


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


def add(rows: list[dict], celltype: str, modality: str, analysis: str,
        metric: str, value, model: str = "", train_domain: str = "",
        test_domain: str = "", stratum: str = "", seed="") -> None:
    if value is None or value == "":
        return
    try:
        v = float(value)
    except (TypeError, ValueError):
        return
    if not np.isfinite(v):
        return
    rows.append({
        "celltype": celltype, "modality": modality, "analysis": analysis,
        "metric": metric, "model": model, "train_domain": train_domain,
        "test_domain": test_domain, "stratum": stratum, "seed": seed,
        "value": v,
    })


def subtype_ari(root: Path, celltype: str) -> tuple[float, int]:
    rna = read_csv(root / celltype / "rna" / "labels.csv")
    atac = read_csv(root / celltype / "atac" / "labels.csv")
    if rna is None or atac is None:
        return float("nan"), 0
    rna = rna.set_index("cell")
    atac = atac.set_index("cell")
    common = rna.index.intersection(atac.index)
    if len(common) < 10:
        return float("nan"), len(common)
    return float(adjusted_rand_score(
        rna.loc[common, "subtype"], atac.loc[common, "subtype"])), len(common)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="results/subtype")
    args = ap.parse_args()
    root = Path(args.root)

    transfer = read_csv(root / "transfer" / "transfer.csv")
    jaccard = read_csv(root / "transfer" / "jaccard.csv")
    if transfer is not None and not transfer.empty:
        transfer = transfer.copy()
        transfer["value"] = pd.to_numeric(transfer["value"], errors="coerce")

    rows: list[dict] = []

    # ---- transfer matrix and controls ---------------------------------------
    if transfer is not None and not transfer.empty:
        for _, r in transfer.iterrows():
            add(rows, str(r["celltype"]), str(r["modality"]), "transfer",
                str(r["metric"]), r["value"], model=str(r["model"]),
                train_domain=str(r["train_domain"]),
                test_domain=str(r["test_domain"]),
                stratum=str(r["stratum"]), seed=str(r["seed"]))
        # separability of each synthetic domain, from the within-domain pair
        sep = transfer[(transfer["model"] == "logistic_regression") &
                        (transfer["metric"] == "auc") &
                        (transfer["stratum"] == "overall") &
                        (transfer["train_domain"] == transfer["test_domain"])]
        for (ct, mod, domain), g in sep.groupby(
                ["celltype", "modality", "train_domain"]):
            add(rows, ct, mod, "separability", f"auc_{domain}", g["value"].mean(),
                train_domain=domain, test_domain=domain)

    # ---- variable-feature overlap -------------------------------------------
    if jaccard is not None and not jaccard.empty:
        for _, r in jaccard.iterrows():
            add(rows, str(r["celltype"]), str(r["modality"]), "jaccard",
                "variable_vs_fixed", r["jaccard"])
            add(rows, str(r["celltype"]), str(r["modality"]), "jaccard",
                "n_variable", r["variable_features"])

    celltypes = sorted({p.parent.parent.name for p in
                        (root).glob("*/*/qc.json")})
    for celltype in celltypes:
        ari, n_shared = subtype_ari(root, celltype)
        add(rows, celltype, "both", "concordance", "rna_atac_ari", ari)
        add(rows, celltype, "both", "concordance", "shared_cells", n_shared)
        for modality in ("rna", "atac"):
            qc = read_json(root / celltype / modality / "qc.json")
            if not qc:
                continue
            for metric in ("silhouette_train", "depth_auc_train",
                           "depth_auc_test", "doublet_auc_train",
                           "doublet_auc_test", "condition_auc"):
                add(rows, celltype, modality, "qc", metric, qc.get(metric))
            add(rows, celltype, modality, "qc", "donor_auc", qc.get("donor_auc"))
            for split in ("train", "test"):
                sizes = qc.get(f"subtype_sizes_{split}", {})
                for subtype, n in sizes.items():
                    add(rows, celltype, modality, "qc",
                        f"n_{split}_{subtype}", n, stratum=split)

    table = pd.DataFrame(rows)
    table.to_csv(root / "comparison.csv", index=False)

    # ---- compact summary ---------------------------------------------------
    summary = {"status": "ok" if not table.empty else "no_rows",
               "rows": int(len(table)), "headline": []}
    if not table.empty:
        overall = table[table["stratum"].isin(["", "overall"])]
        summary["headline"] = (
            overall.groupby(["celltype", "modality", "analysis", "metric",
                            "model", "train_domain", "test_domain"],
                           dropna=False)["value"]
            .agg(["mean", "std", "size"]).reset_index()
            .round(4).to_dict(orient="records")
        )
        # direction asymmetry: synthetic -> real minus real -> synthetic
        asym = transfer[(transfer["model"] == "logistic_regression") &
                        (transfer["metric"] == "auc") &
                        (transfer["stratum"] == "overall")] \
            if transfer is not None and not transfer.empty else pd.DataFrame()
        if not asym.empty:
            fwd = asym[(asym["train_domain"] == "syn_alt") &
                        (asym["test_domain"] == "real")].groupby(
                            ["celltype", "modality"])["value"].mean()
            rev = asym[(asym["train_domain"] == "real") &
                        (asym["test_domain"] == "syn_alt")].groupby(
                            ["celltype", "modality"])["value"].mean()
            for key in fwd.index.intersection(rev.index):
                summary.setdefault("asymmetry", []).append({
                    "celltype": key[0], "modality": key[1],
                    "syn_to_real_auc": round(float(fwd[key]), 4),
                    "real_to_syn_auc": round(float(rev[key]), 4),
                    "difference": round(float(fwd[key] - rev[key]), 4),
                })

    (root / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {root}/comparison.csv and summary.json "
          f"({len(table)} rows)")
    if not table.empty:
        print(table.groupby(["celltype", "analysis"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
