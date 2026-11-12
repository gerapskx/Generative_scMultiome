#!/usr/bin/env python3
"""Plan the nested Aged/PD sample-size series.

Writes one sample sheet and one full config per step:

    config/subsets/full.csv        all cohort samples, use=yes
    config/subsets/n02.csv         step 2
    config/subsets/n08.csv         step 8
    config/subsets/n16.csv         step 16
    config/subsets/n24.csv         step 24
    config/scaling/full.yaml       one-time ingest over the whole cohort
    config/scaling/nXX.yaml        one run per step

Donor sets are nested: the donors in n02 are a subset of n08, n08 of n16, and
n16 of n24. The two pilot donors are always included. Counts per condition are
proportional to the cohort (9 Aged, 15 PD) with a minimum of one donor each.

Usage:
    python src/plan_scaling.py [--config config/config.yaml]
                              [--sheet config/samples.csv]
                              [--steps 2,8,16,24]
                              [--root results/scale]
                              [--fixed-step 24]
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
from pathlib import Path

import yaml


def largest_remainder(total: int, weights: dict[str, float], caps: dict[str, int]) -> dict[str, int]:
    """Allocate `total` seats proportionally, at least one per group, never above cap."""
    groups = list(weights)
    if total < len(groups):
        raise SystemExit(f"step {total} is smaller than the {len(groups)} groups")
    base = {g: max(1, int(total * weights[g])) for g in groups}
    # respect caps, then repair
    for g in groups:
        base[g] = min(base[g], caps[g])
    while sum(base.values()) < total:
        # add to the group with the largest fractional part that still has room
        candidates = [g for g in groups if base[g] < caps[g]]
        if not candidates:
            raise SystemExit(f"cannot allocate {total} donors: caps too small")
        g = max(candidates, key=lambda x: total * weights[x] - base[x])
        base[g] += 1
    while sum(base.values()) > total:
        candidates = [g for g in groups if base[g] > 1]
        g = min(candidates, key=lambda x: total * weights[x] - base[x])
        base[g] -= 1
    return base


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--sheet", default="config/samples.csv")
    ap.add_argument("--steps", default="")
    ap.add_argument("--root", default="")
    ap.add_argument("--subsets", default="config/subsets")
    ap.add_argument("--out", default="config/scaling")
    ap.add_argument("--fixed-step", type=int, default=0)
    args = ap.parse_args()

    base = yaml.safe_load(Path(args.config).read_text())
    scaling = base.get("scaling", {})
    steps = [int(s) for s in (args.steps or ",".join(str(s) for s in scaling.get("steps", [2, 8, 16, 24]))).split(",") if s]
    root = args.root or scaling.get("root", "results/scale")
    fixed_step = args.fixed_step or int(scaling.get("fixed_feature_step", 24))
    cohort = scaling.get("cohort", ["Aged", "PD"])
    anchors = set(scaling.get("anchor_donors", []))

    rows = list(csv.DictReader(Path(args.sheet).open()))
    cohort_rows = [r for r in rows if r["group"] in cohort]
    donors = sorted({r["donor"] for r in cohort_rows})
    group_of = {r["donor"]: r["group"] for r in cohort_rows}
    by_group: dict[str, list[str]] = {g: [] for g in cohort}
    for d in donors:
        by_group[group_of[d]].append(d)
    for g in cohort:
        by_group[g].sort()
        missing = anchors & set(by_group[g])
        # anchors are just ordering hints; a missing anchor is ignored
        by_group[g] = sorted(missing) + [d for d in by_group[g] if d not in missing]

    caps = {g: len(v) for g, v in by_group.items()}
    weights = {g: caps[g] / max(sum(caps.values()), 1) for g in cohort}
    if sum(caps.values()) < max(steps):
        raise SystemExit(
            f"cohort has {sum(caps.values())} donors, cannot reach step {max(steps)}"
        )

    subsets = Path(args.subsets)
    out = Path(args.out)
    subsets.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    # ---- full cohort sheet and ingest config --------------------------------
    full_rows = [dict(r, use="yes" if r["donor"] in set(donors) else "no") for r in rows]
    full_sheet = subsets / "full.csv"
    with full_sheet.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(full_rows)

    full_cfg = copy.deepcopy(base)
    full_cfg["results_dir"] = f"{root}/full"
    full_cfg["processed_dir"] = "data/processed/scale/full"
    full_cfg.setdefault("samples", {})["sheet"] = str(full_sheet)
    full_cfg["samples"]["regenerate"] = False
    full_cfg.setdefault("cluster", {})["auto_annotate"] = True
    full_cfg["cluster"]["map"] = f"{root}/full/cluster/celltype_map.csv"
    full_cfg.setdefault("features", {})["fixed_from"] = ""
    full_cfg.setdefault("scaling", {}).update(
        {"enabled": True, "subset": False, "ingest_only": True, "step": 0,
         "root": root, "full_ingest": f"{root}/full/ingest"}
    )
    full_cfg.setdefault("impact", {})["enabled"] = False
    full_cfg.setdefault("utility", {})["enabled"] = False
    (out / "full.yaml").write_text(yaml.safe_dump(full_cfg, sort_keys=False))

    # ---- nested steps ------------------------------------------------------
    plan_rows = []
    for step in steps:
        alloc = largest_remainder(step, weights, caps)
        selected: dict[str, list[str]] = {}
        for g in cohort:
            if alloc[g] > caps[g]:
                raise SystemExit(f"step {step}: {g} needs {alloc[g]} donors, only {caps[g]}")
            selected[g] = by_group[g][: alloc[g]]
        chosen = set().union(*selected.values())

        sheet = subsets / f"n{step:02d}.csv"
        step_rows = [dict(r, use="yes" if r["donor"] in chosen else "no") for r in rows]
        with sheet.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(step_rows)

        cfg = copy.deepcopy(base)
        cfg["results_dir"] = f"{root}/n{step:02d}"
        cfg["processed_dir"] = f"data/processed/scale/n{step:02d}"
        cfg.setdefault("samples", {})["sheet"] = str(sheet)
        cfg["samples"]["regenerate"] = False
        cfg.setdefault("cluster", {})["auto_annotate"] = True
        cfg["cluster"]["map"] = f"{root}/n{step:02d}/cluster/celltype_map.csv"
        cfg.setdefault("features", {})["fixed_from"] = (
            "" if step == fixed_step else f"{root}/n{fixed_step:02d}/features"
        )
        cfg.setdefault("scaling", {}).update(
            {"enabled": True, "subset": True, "ingest_only": False, "step": step,
             "root": root, "full_ingest": f"{root}/full/ingest", "assemble": False}
        )
        cfg.setdefault("impact", {})["enabled"] = True
        cfg.setdefault("utility", {})["enabled"] = True
        # The subtype transfer experiment is only run on the fixed (largest)
        # step, where the n24 feature set and the full donor set exist.
        cfg.setdefault("subtype", {})["enabled"] = step == fixed_step
        # The permutation control refits the full marginal model once per cell
        # type per modality. That is affordable on the 2-donor pilot but not on
        # 24 donors x several strata, so scaling steps keep the cheap model
        # contrast and skip the refit control.
        cfg.setdefault("reversal", {})["permute"] = False
        (out / f"n{step:02d}.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))

        plan_rows.append({
            "step": step,
            "n_donors": len(chosen),
            **{f"n_{g}": len(selected[g]) for g in cohort},
            "donors": ";".join(sorted(chosen)),
            "sheet": str(sheet),
            "config": str(out / f"n{step:02d}.yaml"),
        })
        print(f"n{step:02d}: {len(chosen)} donors "
              f"({', '.join(f'{g}={len(selected[g])}' for g in cohort)})")

    with (out / "plan.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(plan_rows[0]))
        w.writeheader()
        w.writerows(plan_rows)
    (out / "plan.json").write_text(json.dumps(plan_rows, indent=2))
    print(f"wrote {out}/plan.csv and {len(steps)} step configs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
