#!/usr/bin/env python3
"""Reversal scatter plots with the most significant features labelled.

One panel per modality x cell type: model log2FC (PD -> Aged) against the
classifier contribution toward Aged. The features at the top of *either*
ranking are labelled, so both the counterfactual model and the classifier get
a voice on the plot.

    python src/plot_reversal.py --step n24
    python src/plot_reversal.py --step n24 --top 8 --dpi 300

Writes 300 dpi PNGs plus a grid overview and a CSV of what was labelled to
/mnt/c/Users/gerar/OneDrive/Desktop/Reversal/ (override with --outdir).
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"
DESKTOP = Path("/mnt/c/Users/gerar/OneDrive/Desktop")

MODEL_COLOR = "#C44E52"
CLASSIFIER_COLOR = "#2A7F62"
POINT_COLOR = "#6C8EBF"
MODALITY_ORDER = ["rna", "atac"]


def short_feature(name: str) -> str:
    """Compact ATAC peak names: chr16-46386132-46392150 -> chr16:46.39Mb."""
    match = re.fullmatch(r"(chr[\w.]+)[-:](\d+)[-:](\d+)", name)
    if not match:
        return name
    chrom, start, end = match.groups()
    mid = (int(start) + int(end)) / 2
    if mid >= 1_000_000:
        return f"{chrom}:{mid / 1_000_000:.2f}Mb"
    return f"{chrom}:{mid / 1_000:.0f}kb"


def load_summary(path: Path) -> dict[tuple[str, str], float]:
    if not path.exists():
        return {}
    payload = json.loads(path.read_text())
    out = {}
    for row in payload.get("strata", []):
        rho = row.get("spearman_model_vs_classifier")
        if rho is not None:
            out[(row["modality"], row["celltype"])] = float(rho)
    return out


def pick_labels(sub: pd.DataFrame, top: pd.DataFrame, top_n: int):
    """Union of the top features by |log2FC| and by |contribution|."""
    chosen: dict[str, str] = {}
    for method in ("model", "classifier"):
        part = top[top["method"] == method].nsmallest(top_n, "rank")
        for feature in part["feature"]:
            chosen.setdefault(feature, method)
    rows = []
    for feature, method in chosen.items():
        hit = sub[sub["feature"] == feature]
        if hit.empty:
            continue
        rows.append((feature, method, float(hit["log2fc"].iloc[0]),
                     float(hit["contribution"].iloc[0])))
    return rows


def draw_panel(ax, sub, labels, rho, fontsize, point_size, show_axes=True):
    ax.scatter(sub["log2fc"], sub["contribution"], s=point_size, alpha=0.45,
               color=POINT_COLOR, linewidths=0, rasterized=True, zorder=2)
    ax.axhline(0, color="#999999", lw=0.7, zorder=1)
    ax.axvline(0, color="#999999", lw=0.7, zorder=1)

    # A handful of very low-count genes have |log2FC| in the tens; a symlog axis
    # keeps them visible without collapsing the bulk of the distribution.
    if sub["log2fc"].abs().max() > 8:
        ax.set_xscale("symlog", linthresh=2)
        ax.set_xlim(-30, 30)

    texts = []
    for feature, method, x, y in labels:
        colour = MODEL_COLOR if method == "model" else CLASSIFIER_COLOR
        texts.append(
            ax.text(x, y, short_feature(feature), fontsize=fontsize,
                    fontweight="bold", color=colour, zorder=5,
                    path_effects=[pe.withStroke(linewidth=2.5,
                                                foreground="white")])
        )

    if texts:
        from adjustText import adjust_text
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle="-", color="#333333", lw=0.6,
                            shrinkA=2, shrinkB=2),
            expand=(1.15, 1.35),
            force_text=(0.1, 0.25),
            ensure_inside_axes=True,
        )

    if show_axes:
        ax.set_xlabel("model log2FC (PD \u2192 Aged)", fontsize=8)
        ax.set_ylabel("classifier contribution (toward Aged)", fontsize=8)
        ax.tick_params(labelsize=7)
    if rho is not None:
        ax.text(0.02, 0.98, f"$\\rho$ = {rho:.2f}", transform=ax.transAxes,
                ha="left", va="top", fontsize=7, color="#555555")
    ax.margins(x=0.14, y=0.14)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", default="n24")
    ap.add_argument("--top", type=int, default=5,
                    help="features labelled per ranking method (default 5)")
    ap.add_argument("--grid-top", type=int, default=3,
                    help="features labelled per method in the grid overview")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--point-size", type=float, default=8)
    ap.add_argument("--outdir", default=str(DESKTOP / "Reversal"))
    ap.add_argument("--no-grid", action="store_true")
    args = ap.parse_args()

    rev = SCALE / args.step / "reversal"
    combined = pd.read_csv(rev / "combined.csv")
    top = pd.read_csv(rev / "top_features.csv")
    rho_map = load_summary(rev / "summary.json")

    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    celltypes = sorted(combined["celltype"].unique())
    modalities = [m for m in MODALITY_ORDER if m in set(combined["modality"])]

    manifest = []
    for modality in modalities:
        for celltype in celltypes:
            sub = combined[(combined["modality"] == modality)
                           & (combined["celltype"] == celltype)]
            part = top[(top["modality"] == modality)
                       & (top["celltype"] == celltype)]
            if sub.empty or part.empty:
                continue
            labels = pick_labels(sub, part, args.top)
            rho = rho_map.get((modality, celltype))

            fig, ax = plt.subplots(figsize=(5.2, 4.2))
            draw_panel(ax, sub, labels, rho, fontsize=7,
                       point_size=args.point_size)
            ax.set_title(f"{modality.upper()} / {celltype}", fontsize=10)
            fig.tight_layout()

            out = outdir / f"reversal_{modality}_{celltype}_{args.step}.png"
            if out.exists():
                out.unlink()  # Windows viewers lock open files (Errno 22)
            fig.savefig(out, dpi=args.dpi, bbox_inches="tight",
                        facecolor="white")
            plt.close(fig)
            print(f"wrote {out}")

            for feature, method, x, y in labels:
                manifest.append({
                    "modality": modality, "celltype": celltype,
                    "feature": feature, "short": short_feature(feature),
                    "method": method, "log2fc": x, "contribution": y,
                })

    if not args.no_grid:
        fig, axes = plt.subplots(len(modalities), len(celltypes),
                                 figsize=(4.2 * len(celltypes),
                                          3.8 * len(modalities)),
                                 squeeze=False)
        for i, modality in enumerate(modalities):
            for j, celltype in enumerate(celltypes):
                ax = axes[i][j]
                sub = combined[(combined["modality"] == modality)
                               & (combined["celltype"] == celltype)]
                part = top[(top["modality"] == modality)
                           & (top["celltype"] == celltype)]
                labels = pick_labels(sub, part, args.grid_top)
                draw_panel(ax, sub, labels, rho_map.get((modality, celltype)),
                           fontsize=6, point_size=args.point_size * 0.7,
                           show_axes=(i == len(modalities) - 1))
                ax.set_title(f"{modality.upper()} / {celltype}", fontsize=9)
                if i == len(modalities) - 1:
                    ax.set_xlabel("model log2FC (PD \u2192 Aged)", fontsize=7)
                if j == 0:
                    ax.set_ylabel("classifier contribution", fontsize=7)
                ax.tick_params(labelsize=6)

        from matplotlib.lines import Line2D
        handles = [
            Line2D([], [], color=MODEL_COLOR, marker="o", ls="", ms=4,
                   label=f"top {args.grid_top} by model |log2FC|"),
            Line2D([], [], color=CLASSIFIER_COLOR, marker="o", ls="", ms=4,
                   label=f"top {args.grid_top} by classifier |contribution|"),
        ]
        fig.legend(handles=handles, loc="lower center", ncol=2,
                   frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.01))
        fig.suptitle(f"Counterfactual PD \u2192 Aged reversal, {args.step}",
                     fontsize=12)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        out = outdir / f"reversal_grid_{args.step}.png"
        if out.exists():
            out.unlink()
        fig.savefig(out, dpi=args.dpi, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)
        print(f"wrote {out}")

    table = pd.DataFrame(manifest)
    table.to_csv(outdir / "labelled_features.csv", index=False)
    print(f"wrote {outdir / 'labelled_features.csv'} "
          f"({len(table)} labels across "
          f"{table[['modality', 'celltype']].drop_duplicates().shape[0]} panels)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
