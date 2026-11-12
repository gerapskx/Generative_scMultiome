#!/usr/bin/env python3
"""Stage 4 - stratified fit/holdout split.

Stratifies by `celltype` x `donor` so every cell type appears in both halves and
no donor is held out entirely (with a 2-donor pilot a leave-one-donor-out split
would make the holdout unlearnable).

Feature selection is later fit on the fit cells only; the holdout is never used to
choose genes or peaks.

Usage:
    python src/split.py [--obs results/cluster/obs_annotated.csv]
                        [--out results/split] [--fraction 0.3] [--seed 20261110]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--out", default="results/split")
    ap.add_argument("--fraction", type=float, default=0.30)
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--min-fit-per-stratum", type=int, default=1)
    args = ap.parse_args()

    if not 0.0 < args.fraction < 1.0:
        raise SystemExit("--fraction must be in (0, 1)")

    obs = pd.read_csv(args.obs, dtype={"cell": str, "cluster": str})
    for col in ("cell", "celltype", "donor"):
        if col not in obs.columns:
            raise SystemExit(f"{args.obs} missing column {col!r}")
    if obs["cell"].duplicated().any():
        raise SystemExit("duplicated cell ids")

    rng = np.random.default_rng(args.seed)
    split = pd.Series("fit", index=obs.index, dtype=object)

    for (celltype, donor), idx in obs.groupby(["celltype", "donor"], sort=True).groups.items():
        idx = np.asarray(list(idx))
        n = len(idx)
        n_hold = int(round(n * args.fraction))
        # never empty a stratum: keep at least min-fit-per-stratum in fit
        n_hold = min(n_hold, max(0, n - args.min_fit_per_stratum))
        if n_hold == 0:
            continue
        chosen = rng.choice(idx, size=n_hold, replace=False)
        split.loc[chosen] = "holdout"

    out = pd.DataFrame({"cell": obs["cell"], "split": split.values})
    out = out.merge(
        obs[["cell", "celltype", "donor", "group", "cluster"]], on="cell", how="left"
    )

    root = Path(args.out)
    root.mkdir(parents=True, exist_ok=True)
    out.to_csv(root / "split.csv", index=False)
    (root / "fit_cells.txt").write_text(
        "\n".join(out.loc[out["split"] == "fit", "cell"]) + "\n"
    )
    (root / "holdout_cells.txt").write_text(
        "\n".join(out.loc[out["split"] == "holdout", "cell"]) + "\n"
    )

    summary = (
        out.groupby(["split", "celltype", "donor"]).size().rename("cells").reset_index()
    )
    summary.to_csv(root / "split_summary.csv", index=False)

    print(f"wrote {root}/split.csv ({len(out)} cells)")
    print(out["split"].value_counts().to_string())
    missing = set(obs["celltype"]) - set(out.loc[out["split"] == "holdout", "celltype"])
    if missing:
        print(f"WARNING: cell types absent from holdout: {', '.join(sorted(missing))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
