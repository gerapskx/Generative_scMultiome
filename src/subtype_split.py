#!/usr/bin/env python3
"""Donor-disjoint split for the subtype transfer experiment.

The main pipeline split (`results/.../split/split.csv`) is cell-level: fit and
holdout cells come from the same donors. That is appropriate for the condition
analysis, but it leaks donor identity into any synthetic-to-real test. This stage
defines a separate donor-level split, per target cell type:

- donors with fewer than `--min-cells-per-donor` cells in the target type are
  dropped (recorded in `dropped.csv`);
- within each condition, `--fraction` of the remaining donors are held out
  (at least one test donor and at least two training donors per condition);
- the remaining donors are the training donors.

The subtype labels themselves are fit on training donors only by
`R/06_subcluster.R`, so no held-out donor contributes to the subtype definition.

Usage:
    python src/subtype_split.py --obs results/scale/n24/cluster/obs_annotated.csv \
        --out results/scale/n24/subtype \
        --celltypes Oligodendrocyte,Microglia --fraction 0.30 \
        --min-cells-per-donor 20 --seed 20261112

Outputs (out/):
    donors.csv         celltype, donor, group, n_cells, split
    dropped.csv        celltype, donor, group, n_cells, reason
    split_summary.csv  celltype, split, group, donors, cells
"""

from __future__ import annotations

import argparse
import zlib
from pathlib import Path

import numpy as np
import pandas as pd


def stable_seed(seed: int, celltype: str, group: str) -> int:
    return int(seed + zlib.crc32(f"{celltype}|{group}".encode()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--out", default="results/subtype")
    ap.add_argument("--celltypes", default="Oligodendrocyte,Microglia")
    ap.add_argument("--fraction", type=float, default=0.30)
    ap.add_argument("--min-cells-per-donor", type=int, default=20)
    ap.add_argument("--min-train-per-condition", type=int, default=2)
    ap.add_argument("--min-test-per-condition", type=int, default=1)
    ap.add_argument("--seed", type=int, default=20261112)
    args = ap.parse_args()

    if not 0.0 < args.fraction < 1.0:
        raise SystemExit("--fraction must be in (0, 1)")

    obs = pd.read_csv(args.obs, dtype={"cell": str, "donor": str, "group": str})
    for col in ("cell", "celltype", "donor", "group"):
        if col not in obs.columns:
            raise SystemExit(f"{args.obs} missing column {col!r}")
    if obs["cell"].duplicated().any():
        raise SystemExit("duplicated cell ids")

    celltypes = [c.strip() for c in args.celltypes.split(",") if c.strip()]
    unknown = sorted(set(celltypes) - set(obs["celltype"].unique()))
    if unknown:
        raise SystemExit(f"unknown cell types: {', '.join(unknown)}")

    counts = (obs.groupby(["celltype", "donor", "group"], sort=True)
              .size().rename("n_cells").reset_index())

    donor_rows: list[dict] = []
    dropped_rows: list[dict] = []
    summary_rows: list[dict] = []

    for ct in celltypes:
        sub = counts[counts["celltype"] == ct].copy()
        if sub.empty:
            raise SystemExit(f"{ct}: no cells in {args.obs}")

        keep = sub[sub["n_cells"] >= args.min_cells_per_donor].copy()
        drop = sub[sub["n_cells"] < args.min_cells_per_donor].copy()
        for _, r in drop.iterrows():
            dropped_rows.append({
                "celltype": ct, "donor": r["donor"], "group": r["group"],
                "n_cells": int(r["n_cells"]),
                "reason": f"fewer than {args.min_cells_per_donor} cells in {ct}",
            })

        if keep.empty:
            raise SystemExit(f"{ct}: no donor has >= {args.min_cells_per_donor} cells")

        for group, g in keep.groupby("group", sort=True):
            donors = sorted(g["donor"].tolist())
            n = len(donors)
            if n < args.min_train_per_condition + args.min_test_per_condition:
                raise SystemExit(
                    f"{ct}/{group}: {n} eligible donors, need at least "
                    f"{args.min_train_per_condition + args.min_test_per_condition}"
                )
            n_hold = int(np.floor(args.fraction * n + 0.5))
            n_hold = max(args.min_test_per_condition, n_hold)
            n_hold = min(n_hold, n - args.min_train_per_condition)
            rng = np.random.default_rng(stable_seed(args.seed, ct, group))
            hold = set(rng.choice(donors, size=n_hold, replace=False).tolist())
            for _, r in g.iterrows():
                donor_rows.append({
                    "celltype": ct, "donor": r["donor"], "group": r["group"],
                    "n_cells": int(r["n_cells"]),
                    "split": "test" if r["donor"] in hold else "train",
                })

        # one summary row per cell type/split/group
        ct_rows = [r for r in donor_rows if r["celltype"] == ct]
        for split in ("train", "test"):
            for group in sorted(keep["group"].unique()):
                rows = [r for r in ct_rows
                        if r["split"] == split and r["group"] == group]
                summary_rows.append({
                    "celltype": ct, "split": split, "group": group,
                    "donors": len(rows),
                    "cells": int(sum(r["n_cells"] for r in rows)),
                })

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    donors = pd.DataFrame(donor_rows, columns=[
        "celltype", "donor", "group", "n_cells", "split"])
    dropped = pd.DataFrame(dropped_rows, columns=[
        "celltype", "donor", "group", "n_cells", "reason"])
    summary = pd.DataFrame(summary_rows, columns=[
        "celltype", "split", "group", "donors", "cells"])
    donors.to_csv(out / "donors.csv", index=False)
    dropped.to_csv(out / "dropped.csv", index=False)
    summary.to_csv(out / "split_summary.csv", index=False)

    print(f"wrote {out}/donors.csv ({len(donors)} donors)")
    print(summary.to_string(index=False))
    if not dropped.empty:
        print("\ndropped donors:")
        print(dropped.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
