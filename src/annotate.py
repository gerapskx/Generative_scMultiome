#!/usr/bin/env python3
"""Stage 3 - apply the manual annotation gate.

`R/01_cluster.R` writes clusters and a template `config/celltype_map.csv`. A human
fills in `celltype` per cluster. This script enforces the gate: every cluster must
have a label before anything downstream can run.

Usage:
    python src/annotate.py [--clusters results/cluster/clusters.csv]
                           [--map config/celltype_map.csv]
                           [--out results/cluster/obs_annotated.csv]
                           [--allow-unassigned]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", default="results/cluster/clusters.csv")
    ap.add_argument("--map", default="config/celltype_map.csv")
    ap.add_argument("--out", default="results/cluster/obs_annotated.csv")
    ap.add_argument("--allow-unassigned", action="store_true")
    args = ap.parse_args()

    clusters = pd.read_csv(args.clusters, dtype={"cell": str, "cluster": str})
    mapping = pd.read_csv(args.map, dtype=str, keep_default_na=False)

    if "cluster" not in mapping.columns or "celltype" not in mapping.columns:
        raise SystemExit(f"{args.map} needs `cluster` and `celltype` columns")

    mapping["cluster"] = mapping["cluster"].str.strip()
    mapping["celltype"] = mapping["celltype"].str.strip()
    mapping = mapping.drop_duplicates("cluster")

    observed = set(clusters["cluster"].unique())
    mapped = set(mapping.loc[mapping["celltype"] != "", "cluster"])

    missing = sorted(observed - mapped)
    if missing:
        if not args.allow_unassigned:
            detail = mapping[mapping["cluster"].isin(missing)]
            print(
                f"annotation gate: {len(missing)} cluster(s) have no celltype: "
                f"{', '.join(missing)}",
                file=sys.stderr,
            )
            print(f"edit {args.map} and re-run, or pass --allow-unassigned",
                  file=sys.stderr)
            print(detail.to_string(index=False), file=sys.stderr)
            return 1
        print(
            f"WARNING: {len(missing)} cluster(s) left unassigned -> celltype=UN: "
            f"{', '.join(missing)}",
            file=sys.stderr,
        )

    extra = sorted(mapped - observed)
    if extra:
        print(f"WARNING: map has clusters not in the data: {', '.join(extra)}",
              file=sys.stderr)

    lookup = dict(zip(mapping["cluster"], mapping["celltype"]))
    clusters["celltype"] = [
        (lookup.get(c) or "UN") for c in clusters["cluster"]
    ]

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    clusters.to_csv(out, index=False)

    counts = clusters["celltype"].value_counts()
    print(f"wrote {out} ({len(clusters)} cells, {counts.size} cell types)")
    for name, n in counts.items():
        print(f"  {name:<20} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
