#!/usr/bin/env python3
"""Subset the one-time full ingest to a scaling step's donors.

The full ingest builds one consensus peak space and recounts every fragment once.
Re-running the recount per step would cost hours, so each step only filters the
shared matrices. Column order and the RNA/ATAC cell order are preserved.

Usage:
    python src/subset_ingest.py --full results/scale/full \\
                                --samples config/subsets/n08.csv \\
                                --out results/scale/n08/ingest
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from pathlib import Path


def selected_donors(sheet: Path) -> tuple[set[str], list[str]]:
    rows = list(csv.DictReader(sheet.open()))
    use = {r["donor"] for r in rows if r.get("use", "").strip().lower() in
           ("yes", "true", "t", "1", "y")}
    samples = [r["sample"] for r in rows if r["donor"] in use]
    if not use:
        raise SystemExit(f"no use=yes rows in {sheet}")
    return use, samples


def subset_mtx(src: Path, dst: Path, keep: dict[int, int]) -> tuple[int, int, int]:
    """Stream a MatrixMarket file, keeping the columns in `keep` (old -> new)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    body = dst.with_suffix(dst.suffix + ".body")
    kept = 0
    n_rows = 0
    with gzip.open(src, "rt") as fh, body.open("w") as out:
        header = fh.readline()
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_rows, n_cols, _ = (int(v) for v in line.split())
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            row, col, val = int(parts[0]), int(parts[1]), parts[2]
            new_col = keep.get(col - 1)
            if new_col is None:
                continue
            out.write(f"{row} {new_col + 1} {val}\n")
            kept += 1
    with gzip.open(dst, "wt") as out:
        out.write(header)
        out.write(f"{n_rows} {len(keep)} {kept}\n")
        with body.open("r") as fh:
            shutil.copyfileobj(fh, out)
    body.unlink()
    return n_rows, len(keep), kept


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    full = Path(args.full)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    donors, samples = selected_donors(Path(args.samples))

    with (full / "obs.csv").open() as fh:
        obs = list(csv.DictReader(fh))
    if not obs:
        raise SystemExit(f"empty obs.csv in {full}")
    fields = list(obs[0])
    keep_rows = [i for i, r in enumerate(obs) if r["donor"] in donors]
    if not keep_rows:
        raise SystemExit(f"no cells for donors {sorted(donors)}")
    keep = {old: new for new, old in enumerate(keep_rows)}

    with (out / "obs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for i in keep_rows:
            w.writerow(obs[i])

    shapes = {}
    for modality in ("rna", "atac"):
        src = full / f"{modality}_counts.mtx.gz"
        dst = out / f"{modality}_counts.mtx.gz"
        shapes[modality] = subset_mtx(src, dst, keep)
        print(f"{modality}: {shapes[modality][0]} x {shapes[modality][1]} "
              f"({shapes[modality][2]} nonzeros)")

    for name in ("rna_features.tsv", "atac_features.tsv", "consensus_peaks.bed",
                  "sessionInfo.txt"):
        shutil.copyfile(full / name, out / name)

    qc = list(csv.DictReader((full / "qc_per_sample.csv").open()))
    qc_fields = list(qc[0]) if qc else []
    with (out / "qc_per_sample.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=qc_fields)
        w.writeheader()
        for r in qc:
            if r["sample"] in samples:
                w.writerow(r)

    summary = {
        "full": str(full),
        "samples": sorted(samples),
        "donors": sorted(donors),
        "cells": len(keep_rows),
        "rna": {"rows": shapes["rna"][0], "cols": shapes["rna"][1],
                "nonzeros": shapes["rna"][2]},
        "atac": {"rows": shapes["atac"][0], "cols": shapes["atac"][1],
                 "nonzeros": shapes["atac"][2]},
    }
    (out / "subset_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {out} ({len(keep_rows)} cells, {len(donors)} donors)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
