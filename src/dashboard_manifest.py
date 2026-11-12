#!/usr/bin/env python3
"""Stage 10 - dashboard artifact manifest.

The Streamlit app reads `results/` directly. This script records which artifacts
exist and a few headline numbers, so `snakemake` has a concrete output for the
dashboard stage and a CI job can assert the app has data to show.

Usage: python src/dashboard_manifest.py [--results results] [--out results/dashboard/manifest.json]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

ARTIFACTS = {
    "ingest": "ingest/obs.csv",
    "cluster": "cluster/clusters.csv",
    "annotation": "cluster/obs_annotated.csv",
    "split": "split/split.csv",
    "features_rna": "features/rna_features.tsv",
    "features_atac": "features/atac_features.tsv",
    "simulation": "simulation/manifest.csv",
    "fidelity": "fidelity/summary.json",
    "hypothesis": "hypothesis/summary.json",
    "reversal": "reversal/summary.json",
    "auto_annotation": "cluster/celltype_map.csv",
    "impact": "impact/summary.json",
    "utility": "utility/summary.json",
    "utility_curve": "utility/learning_curve.csv",
    "subtype_transfer": "subtype/transfer/transfer.csv",
    "subtype_comparison": "subtype/comparison.csv",
    "subtype_summary": "subtype/summary.json",
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="results/dashboard/manifest.json")
    args = ap.parse_args()

    root = Path(args.results)
    present = {
        name: {"path": str(root / rel), "exists": (root / rel).exists()}
        for name, rel in ARTIFACTS.items()
    }

    headline: dict = {}
    obs_path = root / "cluster" / "obs_annotated.csv"
    if obs_path.exists():
        obs = pd.read_csv(obs_path)
        headline["n_cells"] = int(len(obs))
        headline["celltypes"] = obs["celltype"].value_counts().to_dict()
        headline["conditions"] = obs["group"].value_counts().to_dict()

    sim_path = root / "simulation" / "manifest.csv"
    if sim_path.exists():
        sim = pd.read_csv(sim_path, keep_default_na=False)
        headline["n_simulations"] = int(len(sim))
        headline["sim_modes"] = sorted(sim["mode"].unique().tolist())
        headline["sim_seeds"] = sorted(int(s) for s in sim["seed"].unique())

    impact_path = root / "impact" / "celltype_impact.csv"
    if impact_path.exists():
        impact = pd.read_csv(impact_path)
        headline["impact_celltypes"] = int(impact["celltype"].nunique())
        headline["impact_top"] = (
            impact.sort_values(["modality", "de_burden"], ascending=[True, False])
            .groupby("modality").head(1)[["modality", "celltype", "de_burden"]]
            .to_dict(orient="records")
        )

    utility_path = root / "utility" / "summary.json"
    if utility_path.exists():
        try:
            headline["utility"] = json.loads(utility_path.read_text())
        except Exception:
            pass

    subtype_path = root / "subtype" / "summary.json"
    if subtype_path.exists():
        try:
            headline["subtype"] = json.loads(subtype_path.read_text())
        except Exception:
            pass

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"artifacts": present, "headline": headline}, indent=2))
    print(f"wrote {out} ({sum(v['exists'] for v in present.values())}/"
          f"{len(present)} artifacts present)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
