#!/usr/bin/env python3
"""Publication figures for the cell-type impact analysis.

Figures, all 300 dpi PNG:

    impact_dumbbell_<step>.png      DE burden and donor-held-out AUC per cell type
    impact_burden_vs_auc.png        burden vs separability, one panel per modality
    impact_heatmap.png              burden and AUC by cell type x donor count
    impact_de_features_<step>_<modality>.png   top DE features per cell type

Two numbers disagree throughout this analysis and the figures are built to make
that visible:

    de_burden   fraction of testable features passing FDR in a cell-level
                PD-vs-Aged test. Its power grows with the number of cells.
    donor_auc   a condition classifier trained on some donors and tested on
                held-out donors, within one cell type. This is the only
                condition test in the table; 0.5 is chance.

A cell type can therefore have a large burden and an AUC at chance (Microglia
ATAC at n24: 0.81 vs 0.53). The dumbbell and the scatter put both on one row or
one point instead of asking the reader to compare two bar panels.

`n02` has one donor per condition, so `donor_auc` is undefined and every
cell-level contrast is also a donor contrast. The heatmaps start at n08 and say
so in the footer.

Usage:
    python src/plot_impact.py
    python src/plot_impact.py --steps n08 n16 n24 --outdir ~/Desktop/Celltype_impact
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import Bbox

from plot_utility import DESKTOP, MODALITY_COLORS, save

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"
CHANCE = 0.5
BURDEN_THRESHOLD = 0.4
DIRECTION_UP = "#C44E52"    # higher in PD
DIRECTION_DOWN = "#4C72B0"  # lower in PD (higher in Aged)
NEUTRAL = "#666666"


# ------------------------------------------------------------------ helpers --
def numeric(frame: pd.DataFrame, columns) -> pd.DataFrame:
    frame = frame.copy()
    for col in columns:
        if col in frame.columns:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
    return frame


def _cells(value) -> str:
    try:
        n = int(float(value))
    except (TypeError, ValueError):
        return "?"
    return f"{n / 1000:.1f}k" if n >= 1000 else str(n)


def load_step(step: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = SCALE / step / "impact"
    impact = pd.read_csv(base / "celltype_impact.csv")
    de = pd.read_csv(base / "de_features.csv")
    return impact, de


def comparison_impact() -> pd.DataFrame:
    """Long `comparison.csv` rows for the impact analysis, numeric where useful."""
    path = SCALE / "comparison.csv"
    table = pd.read_csv(path)
    imp = table[table["analysis"] == "impact"].copy()
    imp["value"] = pd.to_numeric(imp["value"], errors="coerce")
    imp["N"] = pd.to_numeric(imp["N"], errors="coerce")
    return imp


def metric_at(imp: pd.DataFrame, n, modality: str, celltype: str,
              metric: str) -> float:
    part = imp[(imp["N"] == n) & (imp["modality"] == modality)
               & (imp["celltype"] == celltype) & (imp["metric"] == metric)]
    if part.empty:
        return float("nan")
    return float(part["value"].mean())


def _place_labels(ax, points, fig, fontsize: float = 7.0) -> None:
    """Annotate (x, y, label) points without letting the labels overlap.

    Twelve labels on a 10 x 5 inch panel are enough to collide, and the points
    are not spread evenly, so each label tries a short list of offsets and keeps
    the first one that misses every previously placed label and every marker.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    placed: list[Bbox] = []

    # Reserve the markers themselves: a label must not sit on top of a point.
    for x, y, _ in points:
        px, py = ax.transData.transform((x, y))
        placed.append(Bbox.from_extents(px - 5, py - 5, px + 5, py + 5))

    candidates = ((6, 4, "left"), (6, -10, "left"), (-6, 4, "right"),
                  (-6, -10, "right"), (0, 9, "center"), (0, -14, "center"))
    for x, y, label in points:
        chosen = None
        for dx, dy, ha in candidates:
            text = ax.annotate(label, (x, y), xytext=(dx, dy),
                               textcoords="offset points", fontsize=fontsize,
                               ha=ha, va="bottom", color="#222222")
            fig.canvas.draw()
            box = text.get_window_extent(renderer=renderer).expanded(1.04, 1.18)
            if not any(box.overlaps(other) for other in placed):
                chosen = text
                break
            text.remove()
        if chosen is None:
            chosen = ax.annotate(label, (x, y), xytext=(6, 4),
                                 textcoords="offset points", fontsize=fontsize,
                                 ha="left", va="bottom", color="#222222")
            fig.canvas.draw()
        placed.append(chosen.get_window_extent(renderer=renderer))


def _style(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(labelsize=7)


# ------------------------------------------------------------------ dumbbell --
def build_dumbbell(impact: pd.DataFrame, title: str | None = None):
    """One row per cell type: hollow = DE burden, filled = donor-held-out AUC.

    Both live on 0-1, so they can share an axis. Rows are sorted by AUC when it
    exists; the connector length is the disagreement between the two metrics.
    """
    impact = numeric(impact, ("de_burden", "donor_auc", "n_cells"))
    impact = impact.dropna(subset=["de_burden"])
    has_auc = bool(np.isfinite(impact["donor_auc"]).any())
    mods = [m for m in ("rna", "atac") if (impact["modality"] == m).any()]
    if not mods:
        raise ValueError("impact frame has no rna/atac rows")

    n_rows = max(len(impact[impact["modality"] == m]) for m in mods)
    fig, axes = plt.subplots(1, len(mods), squeeze=False,
                             figsize=(11.4, 0.48 * n_rows + 2.0))
    for ax, modality in zip(axes[0], mods):
        part = impact[impact["modality"] == modality].copy()
        sort_key = "donor_auc" if has_auc else "de_burden"
        part = part.sort_values(sort_key, ascending=True, na_position="first")
        y = np.arange(len(part))
        colour = MODALITY_COLORS[modality]

        if has_auc:
            ax.hlines(y, part["de_burden"], part["donor_auc"],
                      color="#CCCCCC", lw=1.4, zorder=1)
            ax.scatter(part["de_burden"], y, s=54, facecolors="white",
                       edgecolors=colour, linewidths=1.8, zorder=3,
                       label="DE burden")
            ax.scatter(part["donor_auc"], y, s=54, color=colour, zorder=4,
                       label="donor-held-out AUC")
            for yi, burden, auc in zip(y, part["de_burden"], part["donor_auc"]):
                if np.isfinite(burden) and np.isfinite(auc) \
                        and abs(auc - burden) >= 0.15:
                    ax.annotate(f"{auc - burden:+.2f}",
                                ((burden + auc) / 2, yi),
                                xytext=(0, 5), textcoords="offset points",
                                ha="center", va="bottom", fontsize=6,
                                color="#666666")
        else:
            ax.barh(y, part["de_burden"], color=colour, alpha=0.78, height=0.6)
            ax.set_title(f"{modality.upper()} \u2014 AUC not estimable",
                         fontsize=9, fontweight="bold")

        labels = [f"{ct} ({_cells(n)})" for ct, n in
                  zip(part["celltype"], part["n_cells"])]
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7.5)
        ax.set_ylim(-0.6, len(part) - 0.4)
        ax.set_xlim(0, 1.0)
        ax.axvline(CHANCE, color="#AAAAAA", lw=0.9, ls=(0, (1, 2)), zorder=0)
        ax.set_xlabel("value (0\u20131)", fontsize=8)
        if has_auc:
            ax.set_title(f"{modality.upper()}", fontsize=9, fontweight="bold")
        _style(ax)

    handles = [
        Line2D([], [], marker="o", ls="", ms=6, markerfacecolor="white",
               markeredgecolor="#333333", markeredgewidth=1.6,
               label="DE burden (cell-level, power grows with n)"),
        Line2D([], [], marker="o", ls="", ms=6, color="#333333",
               label="donor-held-out AUC (condition test)"),
        Line2D([], [], color="#AAAAAA", lw=0.9, ls=(0, (1, 2)), label="chance"),
    ]
    if not has_auc:
        handles = [handles[0], handles[2]]
    if not has_auc:
        fig.text(0.5, 0.135,
                 "One donor per condition: the AUC is not estimable and "
                 "cell-level DE is a donor contrast.",
                 ha="center", fontsize=7.2, color="#C44E52")
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=7.4, bbox_to_anchor=(0.5, 0.01))
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0.19 if not has_auc else 0.06, 1,
                           0.96 if title else 1.0))
    return fig


# ------------------------------------------------------------------- scatter --
def build_scatter(impact: pd.DataFrame, title: str | None = None):
    """Burden on x, donor-held-out AUC on y, one panel per modality.

    Marker area scales with the number of cells. High x / low y is the
    "broad but donor-driven" corner (Microglia ATAC at n24); high x / high y is
    a broad and generalisable effect.
    """
    impact = numeric(impact, ("de_burden", "donor_auc", "n_cells"))
    impact = impact.dropna(subset=["de_burden", "donor_auc"])
    if impact.empty:
        return None
    mods = [m for m in ("rna", "atac") if (impact["modality"] == m).any()]
    if not mods:
        return None

    lo_x = min(0.0, float(impact["de_burden"].min()) - 0.05)
    hi_x = float(impact["de_burden"].max()) + 0.10
    lo_y = min(CHANCE, float(impact["donor_auc"].min())) - 0.06
    hi_y = max(CHANCE, float(impact["donor_auc"].max())) + 0.06
    lo_x = max(lo_x, -0.05)

    fig, axes = plt.subplots(1, len(mods), squeeze=False, sharey=True,
                             figsize=(5.6 * len(mods) + 0.6, 4.9))
    label_points: dict[int, list[tuple[float, float, str]]] = {}
    for ax, modality in zip(axes[0], mods):
        part = impact[impact["modality"] == modality]
        colour = MODALITY_COLORS[modality]
        cells = part["n_cells"].to_numpy(dtype=float)
        sizes = 40 + 300 * np.sqrt(cells / np.nanmax(cells))
        ax.scatter(part["de_burden"], part["donor_auc"], s=sizes, color=colour,
                   alpha=0.8, edgecolors="white", linewidths=1.0, zorder=3)

        ax.set_xlim(lo_x, hi_x)
        ax.set_ylim(lo_y, hi_y)
        # Explicit ticks inside the view: the default locator also emits labels
        # just outside the limits, and the right-hand one of the left panel
        # collides with the left-hand one of the right panel.
        ax.set_xticks(np.arange(0.0, hi_x + 1e-9, 0.2))
        ax.set_yticks(np.arange(np.ceil(lo_y / 0.1) * 0.1, hi_y + 1e-9, 0.1))
        # Light shading of the two informative corners.
        ax.add_patch(Rectangle((BURDEN_THRESHOLD, CHANCE), hi_x - BURDEN_THRESHOLD,
                               hi_y - CHANCE, facecolor="#EAF2EA",
                               edgecolor="none", zorder=0))
        ax.add_patch(Rectangle((BURDEN_THRESHOLD, lo_y), hi_x - BURDEN_THRESHOLD,
                               CHANCE - lo_y, facecolor="#F7EDEA",
                               edgecolor="none", zorder=0))
        ax.axhline(CHANCE, color="#AAAAAA", lw=0.9, ls=(0, (1, 2)), zorder=1)
        ax.axvline(BURDEN_THRESHOLD, color="#DDDDDD", lw=0.8, zorder=1)
        ax.set_xlabel("DE burden (significant / testable)", fontsize=8)
        if modality == mods[0]:
            ax.set_ylabel("donor-held-out AUC", fontsize=8)
        ax.set_title(f"{modality.upper()}", fontsize=9, fontweight="bold")
        _style(ax)

        label_points[id(ax)] = [
            (float(r.de_burden), float(r.donor_auc), str(r.celltype))
            for r in part.itertuples()
        ]

    handles = [
        Line2D([], [], marker="o", ls="", ms=6,
               color=MODALITY_COLORS[m], label=m.upper()) for m in mods
    ]
    handles.append(Line2D([], [], color="#AAAAAA", lw=0.9, ls=(0, (1, 2)),
                          label="chance"))
    fig.legend(handles=handles, loc="lower center", ncol=len(handles),
               frameon=False, fontsize=7.6, bbox_to_anchor=(0.5, 0.0))
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0.07, 1, 0.95 if title else 1.0))
    # Labels are placed after the layout settles: tight_layout moves the axes,
    # so a collision check before it would validate a geometry that no longer
    # exists in the saved figure.
    for ax in axes[0]:
        _place_labels(ax, label_points[id(ax)], fig)
    return fig


# ------------------------------------------------------------------ heatmaps --
def _heatmap(ax, matrix, rows, cols, title, cmap, vmin, vmax, centre,
             fmt=".2f"):
    cm = matplotlib.colormaps[cmap].copy()
    cm.set_bad("#E9E9E9")
    im = ax.imshow(matrix, cmap=cm, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, fontsize=7.5)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows, fontsize=7.5)
    ax.set_title(title, fontsize=9, fontweight="bold", pad=7)
    span = vmax - vmin
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if not np.isfinite(value):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=6.5,
                        color="#999999")
                continue
            colour = ("white" if abs(value - centre) > span * 0.34
                      else "#222222")
            ax.text(j, i, format(value, fmt), ha="center", va="center",
                    fontsize=7.0, color=colour)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)
    return im


def build_heatmaps(imp: pd.DataFrame, title: str | None = None):
    """2 x 2: modality (rows) x metric (columns), cell type on y, donors on x."""
    if imp.empty:
        return None
    ns = sorted(int(n) for n in imp["N"].dropna().unique())
    rank = (imp[imp["metric"] == "donor_auc"]
            .groupby("celltype")["value"].mean()
            .sort_values(ascending=False))
    celltypes = [ct for ct in rank.index
                 if ct in set(imp["celltype"].dropna())]
    if not celltypes:
        celltypes = sorted(imp["celltype"].dropna().unique())

    fig, axes = plt.subplots(2, 2, squeeze=False,
                             figsize=(11.6, 1.05 * len(celltypes) + 2.9))
    cols = [f"n{n:02d}" for n in ns]
    specs = (
        ("de_burden", "DE burden", "Blues", 0.0, 0.9, 0.0,
         "fraction significant"),
        ("donor_auc", "donor-held-out AUC", "RdBu_r", 0.35, 0.75, CHANCE,
         "AUC (red > 0.5)"),
    )
    for i, modality in enumerate(("rna", "atac")):
        for j, (metric, label, cmap, vmin, vmax, centre, cbar_label) in \
                enumerate(specs):
            ax = axes[i][j]
            matrix = np.full((len(celltypes), len(ns)), np.nan)
            for r, ct in enumerate(celltypes):
                for k, n in enumerate(ns):
                    matrix[r, k] = metric_at(imp, n, modality, ct, metric)
            im = _heatmap(ax, matrix, celltypes, cols,
                          f"{modality.upper()} \u2014 {label}",
                          cmap, vmin, vmax, centre)
            cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
            cbar.ax.tick_params(labelsize=6.2)
            cbar.set_label(cbar_label, fontsize=6.8)
    fig.text(0.5, 0.012,
             "n02 is excluded: one donor per condition, so the AUC is not "
             "estimable and cell-level DE is a donor contrast.",
             ha="center", fontsize=7, color="#666666")
    fig.suptitle(title or "Cell-type impact across cohorts", fontsize=12)
    fig.tight_layout(rect=(0, 0.035, 1, 0.94))
    return fig


# --------------------------------------------------------------- DE features --
def _feature_panel(ax, part, modality, top_n):
    part = part.copy()
    part = part[(part["significant"]) & (part["testable"])]
    if part.empty:
        ax.axis("off")
        ax.text(0.5, 0.5, "no significant features\n(FDR > 0.05)",
                ha="center", va="center", fontsize=8, color="#888888",
                transform=ax.transAxes)
        return
    part["abs"] = part["log2fc"].abs()
    part = part.sort_values("abs", ascending=False).head(top_n)
    part = part.sort_values("log2fc")
    y = np.arange(len(part))
    colours = [DIRECTION_UP if v > 0 else DIRECTION_DOWN
               for v in part["log2fc"]]
    # Marker area encodes significance: smaller q -> larger dot.
    strength = np.clip(-np.log10(part["q"].clip(lower=1e-12)), 0, 12)
    sizes = 16 + 110 * strength / 12
    ax.hlines(y, 0, part["log2fc"], color="#CFCFCF", lw=1.0, zorder=1)
    ax.scatter(part["log2fc"], y, s=sizes, c=colours, zorder=3,
               edgecolors="white", linewidths=0.5)
    ax.set_yticks(y)
    ax.set_yticklabels(part["feature"], fontsize=6.4)
    ax.axvline(0, color=NEUTRAL, lw=0.8, zorder=0)
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.margins(x=0.12)
    ax.set_xlabel("log2 FC (PD \u2212 Aged)", fontsize=7)
    _style(ax)


def build_de_features(de: pd.DataFrame, celltypes, modality: str,
                      top_n: int = 15, title: str | None = None):
    """Lollipop of the top |log2FC| significant features, one panel per type."""
    de = numeric(de, ("log2fc", "p", "q"))
    de["significant"] = de["significant"].astype(str).str.lower().isin(
        ("true", "1"))
    de["testable"] = de["testable"].astype(str).str.lower().isin(
        ("true", "1"))
    celltypes = list(celltypes)
    ncols = 1 if len(celltypes) == 1 else 2
    nrows = max(1, math.ceil(len(celltypes) / ncols))
    panel_h = 0.30 * top_n + 1.5
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(5.6 * ncols + 0.4, panel_h * nrows))
    for idx, ax in enumerate(axes.ravel()):
        if idx >= len(celltypes):
            ax.axis("off")
            continue
        ct = celltypes[idx]
        part = de[(de["modality"] == modality) & (de["celltype"] == ct)]
        _feature_panel(ax, part, modality, top_n)
        n_sig = int(part["significant"].sum()) if not part.empty else 0
        ax.set_title(f"{ct} ({n_sig} significant)", fontsize=8.5,
                     fontweight="bold")

    handles = [
        Line2D([], [], marker="o", ls="", ms=6, color=DIRECTION_UP,
               label="higher in PD"),
        Line2D([], [], marker="o", ls="", ms=6, color=DIRECTION_DOWN,
               label="lower in PD (higher in Aged)"),
        Line2D([], [], marker="o", ls="", ms=3.5, color="#999999",
               label="marker area = significance"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=7.2, bbox_to_anchor=(0.5, 0.0))
    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95 if title else 1.0))
    return fig


# ----------------------------------------------------------------------- CLI --
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steps", nargs="+", default=["n08", "n16", "n24"])
    ap.add_argument("--outdir", default=str(DESKTOP / "Celltype_impact"))
    args = ap.parse_args()

    steps = [s if s.startswith("n") else f"n{int(s):02d}" for s in args.steps]
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    for step in steps:
        impact, de = load_step(step)
        save(build_dumbbell(impact,
                            title=f"Cell-type impact \u2014 {step}"),
             outdir / f"impact_dumbbell_{step}.png")
        for modality in ("rna", "atac"):
            celltypes = sorted(de.loc[de["modality"] == modality,
                                      "celltype"].unique())
            if not celltypes:
                continue
            save(build_de_features(
                    de, celltypes, modality, top_n=6,
                    title=f"Top DE features \u2014 {step} {modality.upper()}"),
                 outdir / f"impact_de_features_{step}_{modality}.png")

    imp = comparison_impact()
    if imp.empty:
        print("comparison.csv has no impact rows; skipping cross-cohort panels")
        return 0

    # The scatter needs a cohort where donor AUC exists; n24 is the reference.
    reference = steps[-1]
    impact, _ = load_step(reference)
    fig = build_scatter(impact, title=f"Burden vs separability \u2014 {reference}")
    if fig is not None:
        save(fig, outdir / "impact_burden_vs_auc.png")

    fig = build_heatmaps(imp)
    if fig is not None:
        save(fig, outdir / "impact_heatmap.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
