#!/usr/bin/env python3
"""Publication figures for the synthetic-to-real utility analysis.

Figures, all 300 dpi PNG:

    utility_learning_curve_<step>.png   per cohort, 2x2 (RNA/ATAC x AUC/bal acc)
    utility_by_cohort.png               AUC vs donors + gains over controls
    utility_heatmap.png                 AUC by cohort x synthetic size (+ real)
    utility_delta_heatmap.png           synthetic minus size-matched real
    utility_learning_curve_all_cohorts.png   supplement overview

The learning curve uses the configured synthetic sizes on a log2 axis. The
`real_trained` ceiling, the `naive_poisson` floor and chance are drawn as
reference lines, never as points on the size axis. `real_matched` is the same
number of real training cells as the synthetic curve, so
`synthetic - real_matched` isolates the value of the synthetic data from the
cost of training on fewer cells.

Usage:
    python src/plot_utility.py
    python src/plot_utility.py --steps n08 n16 n24 --outdir ~/Desktop/Utility
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedFormatter, FixedLocator

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"
DESKTOP = Path("/mnt/c/Users/gerar/OneDrive/Desktop")

METHOD_COLORS = {
    "synthetic": "#4C72B0",
    "real_matched": "#55A868",
    "permuted_labels": "#999999",
    "real_trained": "#C44E52",
    "naive_poisson": "#DD8452",
    "majority": "#BBBBBB",
}
MODALITY_COLORS = {"rna": "#3B6EA5", "atac": "#C47F2A"}
SIZE_TICKS = [500, 1000, 2000, 5000, 10000]
SIZE_LABELS = {500: "500", 1000: "1k", 2000: "2k", 5000: "5k", 10000: "10k"}
CHANCE = 0.5


def load(step: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = SCALE / step / "utility"
    curve = pd.read_csv(base / "learning_curve.csv")
    raw = pd.read_csv(base / "utility.csv")
    folds_path = base / "sim" / "folds.csv"
    folds = pd.read_csv(folds_path) if folds_path.exists() else pd.DataFrame()
    return curve, raw, folds


def method_curve(curve, modality, method, metric):
    part = curve[(curve["modality"] == modality) & (curve["method"] == method)]
    part = part[pd.to_numeric(part["synthetic_size"], errors="coerce") > 0]
    part = part.sort_values("synthetic_size")
    return (pd.to_numeric(part["synthetic_size"], errors="coerce").to_numpy(),
            pd.to_numeric(part[f"{metric}_mean"], errors="coerce").to_numpy(),
            pd.to_numeric(part[f"{metric}_sd"], errors="coerce").to_numpy())


def single(curve, modality, method, metric):
    part = curve[(curve["modality"] == modality) & (curve["method"] == method)]
    if part.empty:
        return float("nan")
    return float(pd.to_numeric(part[f"{metric}_mean"], errors="coerce").iloc[0])


def draw_curve(ax, curve, raw, modality, metric, compact=False):
    label = "AUC" if metric == "auc" else "balanced accuracy"
    sizes = None
    for method, lw, z in (
        ("synthetic", 2.2, 5),
        ("real_matched", 1.6, 4),
        ("permuted_labels", 1.1, 3),
    ):
        x, mean, sd = method_curve(curve, modality, method, metric)
        if not len(x):
            continue
        sizes = x
        colour = METHOD_COLORS[method]
        if np.isfinite(sd).any():
            ax.fill_between(x, mean - sd, mean + sd, color=colour, alpha=0.16,
                            linewidth=0)
        ax.plot(x, mean, color=colour, lw=lw, marker="o",
                ms=4 if not compact else 3, zorder=z)
        if method == "synthetic" and not compact and raw is not None and len(raw):
            seed_rows = raw[(raw["modality"] == modality)
                            & (raw["method"] == "synthetic")
                            & (pd.to_numeric(raw["synthetic_size"],
                                             errors="coerce") > 0)]
            jitter = np.random.default_rng(7)
            raw_col = "auc" if metric == "auc" else "balanced_accuracy"
            for size, part in seed_rows.groupby("synthetic_size"):
                xs = float(size) * np.exp(jitter.uniform(-0.045, 0.045, len(part)))
                ax.scatter(xs, pd.to_numeric(part[raw_col], errors="coerce"),
                           s=7, color=colour, alpha=0.45, linewidths=0, zorder=6)

    ceiling = single(curve, modality, "real_trained", metric)
    floor = single(curve, modality, "naive_poisson", metric)
    ax.axhline(ceiling, color=METHOD_COLORS["real_trained"], lw=1.3, ls="--",
               zorder=2)
    ax.axhline(floor, color=METHOD_COLORS["naive_poisson"], lw=1.3, ls=":",
               zorder=2)
    ax.axhline(CHANCE, color="#AAAAAA", lw=0.9, ls=(0, (1, 2)), zorder=1)

    ax.set_xscale("log", base=2)
    ax.set_xticks(SIZE_TICKS)
    ax.xaxis.set_major_formatter(FixedFormatter([SIZE_LABELS[s] for s in SIZE_TICKS]))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.set_xlim(380, 13000)
    ax.set_xlabel("synthetic training cells", fontsize=8)
    ax.set_ylabel(f"{label} on held-out donors", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.set_title(f"{modality.upper()} \u2014 {label}", fontsize=9,
                 fontweight="bold")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    vals = [ceiling, floor, CHANCE]
    for method in ("synthetic", "real_matched", "permuted_labels"):
        _, mean, sd = method_curve(curve, modality, method, metric)
        if len(mean):
            vals += list(mean - np.nan_to_num(sd))
            vals += list(mean + np.nan_to_num(sd))
    vals = [v for v in vals if np.isfinite(v)]
    ax.set_ylim(min(vals) - 0.035, max(vals) + 0.035)

    if compact or sizes is None or not len(sizes):
        return
    syn_curve = method_curve(curve, modality, "synthetic", metric)[1]
    matched_curve = method_curve(curve, modality, "real_matched", metric)[1]
    if not len(syn_curve):
        return
    syn = float(syn_curve[-1])
    matched = float(matched_curve[-1]) if len(matched_curve) else float("nan")
    lines = [
        f"ceiling \u2212 synthetic   {ceiling - syn:+.2f}",
        f"synthetic \u2212 matched   {syn - matched:+.2f}",
        f"synthetic \u2212 marginal  {syn - floor:+.2f}",
    ]
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes,
            va="top", ha="left", fontsize=6.4, color="#333333",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor="#DDDDDD", alpha=0.85))


def legend_handles(compact=False):
    handles = [
        Line2D([], [], color=METHOD_COLORS["synthetic"], lw=2.2, marker="o",
               ms=4, label="synthetic"),
        Line2D([], [], color=METHOD_COLORS["real_matched"], lw=1.6, marker="o",
               ms=4, label="matched real (same n)"),
        Line2D([], [], color=METHOD_COLORS["permuted_labels"], lw=1.1,
               label="shuffled labels"),
        Line2D([], [], color=METHOD_COLORS["real_trained"], lw=1.3, ls="--",
               label="real donors (ceiling)"),
        Line2D([], [], color=METHOD_COLORS["naive_poisson"], lw=1.3, ls=":",
               label="marginal only (floor)"),
        Line2D([], [], color="#AAAAAA", lw=0.9, ls=(0, (1, 2)), label="chance"),
    ]
    if compact:
        handles = [h for h in handles if h.get_label() != "shuffled labels"]
    return handles


def cohort_label(step: str, folds: pd.DataFrame) -> str:
    try:
        n = int(step.lstrip("n"))
        name = f"{n} donors"
    except ValueError:
        name = step
    if folds is None or folds.empty or "holdout_donors" not in folds:
        return name
    holdout = len(str(folds["holdout_donors"].iloc[0]).split(";"))
    return f"{name} ({holdout} held out)"


def build_learning_curve(step, curve, raw, folds, title=None):
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 7.4))
    for i, modality in enumerate(("rna", "atac")):
        for j, metric in enumerate(("auc", "bal")):
            draw_curve(axes[i][j], curve, raw, modality, metric)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=3,
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.005))
    fig.suptitle(title or
                 f"Synthetic-to-real utility \u2014 {cohort_label(step, folds)}",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    return fig


def comparison_frame() -> pd.DataFrame:
    path = SCALE / "comparison.csv"
    table = pd.read_csv(path)
    util = table[table["analysis"] == "utility"].copy()
    util["value"] = pd.to_numeric(util["value"], errors="coerce")
    util["N"] = pd.to_numeric(util["N"], errors="coerce")
    return util


def value_at(util, n, modality, label, metric, size=None):
    part = util[(util["N"] == n) & (util["modality"] == modality)
                & (util["label"] == label) & (util["metric"] == metric)]
    if size is not None:
        part = part[pd.to_numeric(part["synthetic_size"], errors="coerce") == size]
    if part.empty:
        return float("nan")
    return float(part["value"].mean())


def build_by_cohort(steps, util):
    ns = sorted(int(s.lstrip("n")) for s in steps)
    hold = {}
    for s in steps:
        _, _, folds = load(s)
        hold[int(s.lstrip("n"))] = (
            len(str(folds["holdout_donors"].iloc[0]).split(";"))
            if not folds.empty else 0)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.6))
    for j, modality in enumerate(("rna", "atac")):
        ax = axes[0][j]
        vals_all = []
        for label, size, marker, name in (
            ("synthetic", 10000, "o", "synthetic @ 10k"),
            ("real_matched", 10000, "s", "matched real @ 10k"),
        ):
            colour = METHOD_COLORS[label]
            vals = [value_at(util, n, modality, label, "auc_mean", size)
                    for n in ns]
            vals_all += [v for v in vals if np.isfinite(v)]
            ax.plot(ns, vals, color=colour, lw=2.0, marker=marker, ms=5,
                    label=name)
            if label != "synthetic":
                continue
            # Only the synthetic points are labelled: the size-matched control
            # can sit within a few thousandths of them, and the delta panel
            # below already reports that difference exactly.
            for n, v in zip(ns, vals):
                if np.isfinite(v):
                    ax.annotate(f"{v:.2f}", (n, v), textcoords="offset points",
                                xytext=(0, 7), ha="center", fontsize=6.5,
                                color=colour)

        for label, ls, name in (
            ("real_trained", "--", "real donors (ceiling)"),
            ("naive_poisson", ":", "marginal only (floor)"),
        ):
            vals = [value_at(util, n, modality, label, "auc_mean") for n in ns]
            vals_all += [v for v in vals if np.isfinite(v)]
            ax.plot(ns, vals, color=METHOD_COLORS[label], lw=1.3, ls=ls,
                    label=name)

        # Headroom for the point labels: without it the top annotation is
        # pushed out of the axes and collides with the title.
        lo = min(vals_all + [CHANCE])
        hi = max(vals_all + [CHANCE])
        pad = max(0.02, 0.12 * (hi - lo))
        ax.set_ylim(lo - pad, hi + pad)

        ax.axhline(CHANCE, color="#AAAAAA", lw=0.9, ls=(0, (1, 2)))
        ax.set_xticks(ns)
        ax.set_xticklabels([f"{n}\n({hold[n]} held out)" for n in ns], fontsize=7)
        ax.set_xlim(min(ns) - 1.2, max(ns) + 1.2)
        ax.set_xlabel("donors (N)", fontsize=8)
        ax.set_ylabel("AUC on held-out donors", fontsize=8)
        ax.tick_params(labelsize=7)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

        below = [n for n in ns
                 if value_at(util, n, modality, "real_trained", "auc_mean") < CHANCE]
        # The below-chance flag goes in the title gap: the panels have no free
        # corner, and every in-axes spot collides with a series or its label.
        ax.set_title(f"{modality.upper()} \u2014 AUC by cohort", fontsize=9,
                     fontweight="bold", pad=13 if below else 6)
        if below:
            ax.text(0.5, 1.012,
                    "below chance: " + ", ".join(f"n{n:02d}" for n in below),
                    transform=ax.transAxes, ha="center", va="bottom",
                    fontsize=6.8, color="#C44E52")

    # Shared legends: an in-axes legend would cover the n24 points, and the two
    # bottom panels need a modality key that is not per-panel clutter.
    fig.legend(
        handles=[
            Line2D([], [], color=METHOD_COLORS["synthetic"], lw=2, marker="o",
                   ms=4, label="synthetic @ 10k"),
            Line2D([], [], color=METHOD_COLORS["real_matched"], lw=2,
                   marker="s", ms=4, label="matched real @ 10k"),
            Line2D([], [], color=METHOD_COLORS["real_trained"], lw=1.3, ls="--",
                   label="real donors (ceiling)"),
            Line2D([], [], color=METHOD_COLORS["naive_poisson"], lw=1.3, ls=":",
                   label="marginal only (floor)"),
            Line2D([], [], color="#AAAAAA", lw=0.9, ls=(0, (1, 2)),
                   label="chance"),
        ],
        loc="upper center", ncol=5, frameon=False, fontsize=7.5,
        bbox_to_anchor=(0.5, 0.955),
    )

    for j, (num, den, title) in enumerate((
        ("synthetic", "real_matched", "synthetic \u2212 matched real"),
        ("synthetic", "naive_poisson", "synthetic \u2212 marginal baseline"),
    )):
        ax = axes[1][j]
        width = 0.36
        x = np.arange(len(ns))
        vals_all = [0.0]
        for k, modality in enumerate(("rna", "atac")):
            vals = [value_at(util, n, modality, num, "auc_mean", 10000)
                    - value_at(util, n, modality, den, "auc_mean") for n in ns]
            vals_all += [v for v in vals if np.isfinite(v)]
            bars = ax.bar(x + (k - 0.5) * width, vals, width,
                          color=MODALITY_COLORS[modality], label=modality.upper())
            for rect, v in zip(bars, vals):
                if np.isfinite(v):
                    ax.annotate(f"{v:+.2f}",
                                (rect.get_x() + rect.get_width() / 2, v),
                                textcoords="offset points",
                                xytext=(0, 3 if v >= 0 else -9), ha="center",
                                fontsize=6.3)
        lo, hi = min(vals_all), max(vals_all)
        pad = max(0.008, 0.28 * (hi - lo))
        ax.set_ylim(lo - pad, hi + pad)
        ax.axhline(0, color="#666666", lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{n}\n({hold[n]} held out)" for n in ns], fontsize=7)
        ax.set_ylabel("\u0394 AUC", fontsize=8)
        ax.set_title(f"{title} @ 10k", fontsize=9, fontweight="bold")
        ax.tick_params(labelsize=7)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)

    fig.legend(
        handles=[Patch(facecolor=MODALITY_COLORS[m], label=m.upper())
                 for m in ("rna", "atac")],
        loc="lower center", ncol=2, frameon=False, fontsize=7.5,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle("Utility across cohorts", fontsize=12, y=0.99)
    fig.tight_layout(rect=(0, 0.04, 1, 0.885))
    return fig


def heatmap(ax, matrix, rows, cols, title, cmap, vmin, vmax, fmt=".2f",
            centre=None):
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(cols, fontsize=7)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows, fontsize=7)
    ax.set_title(title, fontsize=9, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix[i, j]
            if not np.isfinite(value):
                continue
            colour = ("white" if abs(value - centre) > (vmax - vmin) * 0.32
                      else "#222222")
            ax.text(j, i, format(value, fmt), ha="center", va="center",
                    fontsize=6.8, color=colour)
    return im


def build_heatmaps(steps, util):
    ns = sorted(int(s.lstrip("n")) for s in steps)
    rows = [f"n{n:02d}" for n in ns]
    cols = ["real\n(all)"] + [SIZE_LABELS[s] for s in SIZE_TICKS]

    fig_auc, axes = plt.subplots(1, 2, figsize=(11.2, 3.9))
    for j, modality in enumerate(("rna", "atac")):
        matrix = np.full((len(ns), len(cols)), np.nan)
        for i, n in enumerate(ns):
            matrix[i, 0] = value_at(util, n, modality, "real_trained", "auc_mean")
            for k, size in enumerate(SIZE_TICKS, start=1):
                matrix[i, k] = value_at(util, n, modality, "synthetic",
                                        "auc_mean", size)
        im = heatmap(axes[j], matrix, rows, cols,
                     f"{modality.upper()} \u2014 AUC", "RdBu_r", 0.35, 0.75,
                     centre=CHANCE)
        cbar = fig_auc.colorbar(im, ax=axes[j], fraction=0.035, pad=0.02)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label("AUC (red > 0.5)", fontsize=7)
        axes[j].axvline(0.5, color="white", lw=3)
    fig_auc.suptitle("Utility by cohort and synthetic size", fontsize=12)
    fig_auc.tight_layout(rect=(0, 0, 1, 0.93))

    fig_delta, axes = plt.subplots(1, 2, figsize=(10.2, 3.9))
    cols = [SIZE_LABELS[s] for s in SIZE_TICKS]
    for j, modality in enumerate(("rna", "atac")):
        matrix = np.full((len(ns), len(cols)), np.nan)
        for i, n in enumerate(ns):
            for k, size in enumerate(SIZE_TICKS):
                matrix[i, k] = (value_at(util, n, modality, "synthetic",
                                         "auc_mean", size)
                                - value_at(util, n, modality, "real_matched",
                                           "auc_mean", size))
        im = heatmap(axes[j], matrix, rows, cols,
                     f"{modality.upper()} \u2014 synthetic \u2212 matched real",
                     "RdBu_r", -0.12, 0.12, fmt="+.2f", centre=0.0)
        cbar = fig_delta.colorbar(im, ax=axes[j], fraction=0.035, pad=0.02)
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label("\u0394 AUC (red = synthetic better)", fontsize=7)
    fig_delta.suptitle("Synthetic minus size-matched real", fontsize=12)
    fig_delta.tight_layout(rect=(0, 0, 1, 0.93))
    return fig_auc, fig_delta


def build_all_cohorts(steps):
    fig, axes = plt.subplots(len(steps), 2, figsize=(9.6, 3.3 * len(steps)),
                             squeeze=False)
    for i, step in enumerate(steps):
        curve, raw, folds = load(step)
        for j, modality in enumerate(("rna", "atac")):
            draw_curve(axes[i][j], curve, raw, modality, "auc", compact=True)
            axes[i][j].set_title(
                f"{step} ({cohort_label(step, folds)}) \u2014 {modality.upper()}",
                fontsize=8.5, fontweight="bold")
    fig.legend(handles=legend_handles(compact=True), loc="lower center",
               ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Synthetic-to-real utility across cohorts (AUC)", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    return fig


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()  # Windows viewers lock open files (Errno 22)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"wrote {path}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steps", nargs="+", default=["n08", "n16", "n24"])
    ap.add_argument("--outdir", default=str(DESKTOP / "Utility"))
    args = ap.parse_args()

    steps = [s if s.startswith("n") else f"n{int(s):02d}" for s in args.steps]
    outdir = Path(args.outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    for step in steps:
        curve, raw, folds = load(step)
        save(build_learning_curve(step, curve, raw, folds),
             outdir / f"utility_learning_curve_{step}.png")

    util = comparison_frame()
    if util.empty:
        print("comparison.csv has no utility rows; skipping cross-cohort panels")
        return 0
    save(build_by_cohort(steps, util), outdir / "utility_by_cohort.png")
    fig_auc, fig_delta = build_heatmaps(steps, util)
    save(fig_auc, outdir / "utility_heatmap.png")
    save(fig_delta, outdir / "utility_delta_heatmap.png")
    save(build_all_cohorts(steps),
         outdir / "utility_learning_curve_all_cohorts.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
