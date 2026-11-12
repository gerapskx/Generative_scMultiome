#!/usr/bin/env python3
"""Pre-extract the selected features from an ingest matrix.

`R/03_simulate.R` and `R/05_utility_sim.R` only need the ~2,000 selected RNA
genes and ~2,000 selected ATAC peaks. Reading the full ATAC matrix (280M
nonzeros for n24) with Matrix::readMM peaks at several GB and repeatedly
triggered the OOM killer on this VM. This script streams the MatrixMarket file
and writes a tiny feature-subset copy, so the R jobs never touch the full
matrix.

It writes, inside the ingest directory:
    <modality>_counts_sel.mtx.gz   rows renumbered 1..K in selected order
    <modality>_features_sel.tsv    feature names, same order as the rows

The selected order is taken from `--features/<modality>_features.tsv`, so the
row order matches what the R scripts expect.

Usage:
    python src/subset_features.py --ingest results/scale/n24/ingest \\
                                  --features results/scale/n24/features
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path


def read_features(path: Path) -> list[str]:
    return [line.strip() for line in path.open() if line.strip()]


def subset_rows_named(
    src: Path,
    dst: Path,
    full_features: list[str],
    keep: dict[str, int],
) -> tuple[int, int, int]:
    """Stream a MatrixMarket file, keeping the rows named in `keep` (name -> new index)."""
    dst.parent.mkdir(parents=True, exist_ok=True)
    body = dst.with_suffix(dst.suffix + ".body")

    kept = 0
    with gzip.open(src, "rt") as fh, body.open("w") as out:
        header = fh.readline()
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_rows, n_cols, _ = (int(v) for v in line.split())
        if n_rows != len(full_features):
            raise SystemExit(
                f"{src}: header has {n_rows} rows but feature file has "
                f"{len(full_features)} lines"
            )
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            row = int(parts[0])
            name = full_features[row - 1]
            new_row = keep.get(name)
            if new_row is None:
                continue
            out.write(f"{new_row + 1} {parts[1]} {parts[2]}\n")
            kept += 1

    with gzip.open(dst, "wt") as out:
        out.write(header)
        out.write(f"{len(keep)} {n_cols} {kept}\n")
        with body.open("r") as fh:
            shutil.copyfileobj(fh, out)
    body.unlink()
    return n_rows, len(keep), kept


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ingest", required=True)
    ap.add_argument("--features", required=True)
    ap.add_argument("--modalities", default="rna,atac")
    ap.add_argument("--force", action="store_true", help="overwrite existing outputs")
    args = ap.parse_args()

    ingest = Path(args.ingest)
    features = Path(args.features)

    for modality in args.modalities.split(","):
        modality = modality.strip()
        if not modality:
            continue
        full_feats = read_features(ingest / f"{modality}_features.tsv")
        selected = read_features(features / f"{modality}_features.tsv")
        if not selected:
            raise SystemExit(f"no selected features in {features}")

        present = set(full_feats)
        keep = {name: i for i, name in enumerate(selected) if name in present}
        missing = [n for n in selected if n not in present]
        if missing:
            print(
                f"{modality}: {len(missing)} selected features are absent from the "
                f"ingest and will be skipped (first: {missing[:3]})",
                file=sys.stderr,
            )
        if len(keep) < 10:
            raise SystemExit(f"{modality}: fewer than 10 selected features present")

        dst = ingest / f"{modality}_counts_sel.mtx.gz"
        feats_out = ingest / f"{modality}_features_sel.tsv"
        if dst.exists() and feats_out.exists() and not args.force:
            print(f"{modality}: {dst} already exists, skipping (use --force to redo)")
            continue

        src = ingest / f"{modality}_counts.mtx.gz"
        if not src.exists():
            src = ingest / f"{modality}_counts.mtx"
        print(f"{modality}: streaming {src} ...", flush=True)
        n_rows, n_kept, nnz = subset_rows_named(src, dst, full_feats, keep)
        feats_out.write_text("\n".join(selected) + "\n")
        print(
            f"{modality}: {n_rows} -> {n_kept} features, {nnz} nonzeros -> {dst}",
            flush=True,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
