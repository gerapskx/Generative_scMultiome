#!/usr/bin/env python3
"""Draw UMAP panels coloured by cell type, sample, donor or condition.

Three modes, all reading the same n24 artifacts:

    python src/plot_umap.py --step n24 --color celltype
    python src/plot_umap.py --step n24 --color sample
    python src/plot_umap.py --step n24 --color celltype --split group

`--style ov` reproduces the omicverse look (`sc.pl.embedding` + adjustText
labels: bold white text, black outline, thin leader lines, no legend).
`--style plain` uses white boxed labels. Categorical columns with more than
`--label-max` levels get a legend instead of on-plot text.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import Bbox

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"
DESKTOP = Path("/mnt/c/Users/gerar/OneDrive/Desktop")

# One colour per cell type, shared with the dashboard so the same population
# keeps the same colour across every figure.
CELLTYPE_COLORS = {
    "Oligodendrocyte": "#4C72B0",
    "Astrocyte": "#DD8452",
    "Microglia": "#55A868",
    "OPC": "#C44E52",
    "Ependymal": "#8172B3",
    "Neuron": "#CCB974",
    "Pericyte": "#64B5CD",
    "Tcell": "#DA8BC3",
}
GROUP_COLORS = {"Aged": "#4C72B0", "PD": "#DD8452"}
FALLBACK = ["#937860", "#8C8C8C", "#B47CC7", "#C7B411"]


def palette(values: list[str], colour_by: str) -> dict[str, str]:
    """Stable colour map: cell type and condition have fixed colours, sample
    uses scanpy's 28-colour default so all 24 donors stay distinct."""
    if colour_by == "celltype":
        out = {}
        spare = iter(FALLBACK)
        for value in values:
            out[value] = CELLTYPE_COLORS.get(value, next(spare, "#7F7F7F"))
        return out
    if colour_by == "group":
        return {v: GROUP_COLORS.get(v, "#7F7F7F") for v in values}
    import seaborn as sns
    base = sns.color_palette("husl", len(values)).as_hex()
    # Permute so neighbouring categories (and the two conditions) are not given
    # adjacent hues: sorted sample IDs would otherwise split Aged/PD by colour.
    order = np.random.default_rng(2026).permutation(len(values))
    return {v: base[order[i]] for i, v in enumerate(values)}


def load_frame(obs: Path, colour_by: str) -> pd.DataFrame:
    frame = pd.read_csv(obs)
    if colour_by == "sample" and "sample" not in frame.columns:
        ingest = obs.parent.parent / "ingest" / "obs.csv"
        key = pd.read_csv(ingest, usecols=["cell", "sample"])
        frame = frame.merge(key, on="cell", how="left")
    return frame.dropna(subset=["umap_1", "umap_2"])


def scatter_panel(ax, sub, colour_by, colours, values, point_size, alpha):
    sizes = sub[colour_by].astype(str).value_counts()
    # Rare populations last so they are not buried under oligodendrocytes.
    for value in sorted(values, key=lambda v: -sizes.get(v, 0)):
        part = sub[sub[colour_by].astype(str) == value]
        ax.scatter(part["umap_1"], part["umap_2"], s=point_size, alpha=alpha,
                   color=colours[value], linewidths=0, rasterized=True,
                   label=value)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal", adjustable="box")
    ax.margins(x=0.10, y=0.05)


def add_legend(ax, fig, colour_by, values, ncol=1):
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return
    fig.legend(handles, labels, loc="center left", bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=7, markerscale=3, ncol=ncol,
               title=colour_by, title_fontsize=8)


def keep_inside(ax, fig, texts, margin=3):
    """Nudge any label that still pokes past the axes back inside."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(renderer)
    for text in texts:
        box = text.get_window_extent(renderer)
        dx = dy = 0.0
        if box.x1 > ax_box.x1 - margin:
            dx = ax_box.x1 - margin - box.x1
        elif box.x0 < ax_box.x0 + margin:
            dx = ax_box.x0 + margin - box.x0
        if box.y1 > ax_box.y1 - margin:
            dy = ax_box.y1 - margin - box.y1
        elif box.y0 < ax_box.y0 + margin:
            dy = ax_box.y0 + margin - box.y0
        if dx or dy:
            px, py = ax.transData.transform(text.get_position())
            text.set_position(ax.transData.inverted().transform((px + dx, py + dy)))


def label_points(ax, fig, labels, fontsize=9, pad=3, step=5, max_radius=140):
    """Place text at each point, nudging outward until no two boxes overlap."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_box = ax.get_window_extent(renderer)
    placed: list[Bbox] = []
    texts = []

    for text, x, y, colour in labels:
        probe = ax.text(x, y, text, fontsize=fontsize, ha="center", va="center")
        fig.canvas.draw()
        bbox = probe.get_window_extent(renderer)
        probe.remove()
        width, height = bbox.width, bbox.height

        cx, cy = ax.transData.transform((x, y))
        chosen = None
        for radius in range(0, max_radius + 1, step):
            if radius == 0:
                offsets = [(0.0, 0.0)]
            else:
                count = max(8, int(2 * math.pi * radius / step))
                offsets = [(radius * math.cos(2 * math.pi * k / count),
                            radius * math.sin(2 * math.pi * k / count))
                           for k in range(count)]
            for dx, dy in offsets:
                px, py = cx + dx, cy + dy
                box = Bbox.from_bounds(px - width / 2 - pad, py - height / 2 - pad,
                                       width + 2 * pad, height + 2 * pad)
                if not (ax_box.x0 <= box.x0 and box.x1 <= ax_box.x1
                        and ax_box.y0 <= box.y0 and box.y1 <= ax_box.y1):
                    continue
                if any(box.overlaps(other) for other in placed):
                    continue
                chosen = (px, py, box)
                break
            if chosen:
                break

        if chosen is None:
            chosen = (cx, cy, Bbox.from_bounds(cx - width / 2, cy - height / 2,
                                               width, height))

        px, py, box = chosen
        placed.append(box)
        dx, dy = px - cx, py - cy
        tx, ty = ax.transData.inverted().transform((px, py))
        if math.hypot(dx, dy) > 10:
            ax.plot([x, tx], [y, ty], color=colour, lw=0.9, alpha=0.9,
                    solid_capstyle="round", zorder=4)
        texts.append(
            ax.text(
                tx, ty, text, fontsize=fontsize, ha="center", va="center",
                color="#1a1a1a", zorder=6, fontweight="medium",
                bbox=dict(boxstyle="round,pad=0.28", facecolor="white",
                          edgecolor=colour, linewidth=1.1, alpha=0.92),
                path_effects=[pe.withStroke(linewidth=2.2, foreground="white")],
            )
        )
    return texts


def plot_ov_labels(ax, fig, frame, colour_by, values, colours, fontsize=8):
    """omicverse style: bold white labels with a black outline, repelled."""
    from adjustText import adjust_text

    texts = []
    for value in values:
        mask = frame[colour_by].astype(str) == value
        x, y = np.median(frame.loc[mask, ["umap_1", "umap_2"]].to_numpy(), axis=0)
        texts.append(
            ax.text(x, y, value, fontsize=fontsize, fontweight="bold", color="white",
                    path_effects=[pe.withStroke(linewidth=2, foreground="black")])
        )
    adjust_text(
        texts, ax=ax,
        arrowprops=dict(arrowstyle="-", color="black", shrinkA=10, shrinkB=10),
        ensure_inside_axes=True,
    )
    keep_inside(ax, fig, texts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", default="n24")
    ap.add_argument("--color", default="celltype",
                    help="celltype, sample, donor, group, cluster")
    ap.add_argument("--split", default="",
                    help="facet by this obs column, e.g. group")
    ap.add_argument("--style", choices=["ov", "plain"], default="ov")
    ap.add_argument("--legend", choices=["auto", "outside", "onplot", "none"],
                    default="auto")
    ap.add_argument("--label", choices=["auto", "none"], default="auto")
    ap.add_argument("--obs", default="", help="override the observation CSV")
    ap.add_argument("--out", default="")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--point-size", type=float, default=1.5)
    ap.add_argument("--alpha", type=float, default=0.55)
    args = ap.parse_args()

    obs = Path(args.obs) if args.obs else (
        SCALE / args.step / "cluster" / "obs_annotated.csv")
    if not obs.exists():
        raise SystemExit(f"missing {obs}")

    frame = load_frame(obs, args.color)
    frame[args.color] = frame[args.color].astype(str)
    values = sorted(frame[args.color].unique())
    colours = palette(values, args.color)

    splits = [("", frame)]
    if args.split:
        frame[args.split] = frame[args.split].astype(str)
        splits = [(v, frame[frame[args.split] == v])
                  for v in sorted(frame[args.split].unique())]

    # Few categories -> labels on the plot; many -> legend instead.
    label = args.label == "auto" and len(values) <= 8
    legend = args.legend
    if legend == "auto":
        legend = "none" if label else "outside"

    if label and len(values) <= 8:
        figsize = (5, 5)
    elif args.split:
        figsize = (5 * len(splits), 5)
    else:
        figsize = (6.2, 5)

    fig, axes = plt.subplots(1, len(splits), figsize=figsize,
                             sharex=True, sharey=True, squeeze=False)
    axes = axes[0]

    for ax, (name, part) in zip(axes, splits):
        scatter_panel(ax, part, args.color, colours, values,
                      args.point_size * (4 if args.style == "ov" else 1),
                      args.alpha)
        if args.split:
            ax.set_title(name, fontsize=10)

    # Layout must be final before labels are placed: tight_layout moves the
    # axes, which would otherwise shift every label off its cluster.
    fig.tight_layout()

    if label:
        for ax, (name, part) in zip(axes, splits):
            if args.style == "ov":
                plot_ov_labels(ax, fig, part, args.color, values, colours,
                                fontsize=7 if args.split else 8)
            else:
                labels = [(v,
                           float(np.median(part.loc[part[args.color] == v, "umap_1"])),
                           float(np.median(part.loc[part[args.color] == v, "umap_2"])),
                           colours[v]) for v in values]
                label_points(ax, fig, labels, fontsize=9)

    if legend == "outside":
        add_legend(axes[-1], fig, args.color, values,
                    ncol=1 if len(values) <= 12 else 2)
    elif legend == "onplot":
        axes[-1].legend(frameon=False, fontsize=7, markerscale=4)

    suffix = "_ovstyle" if args.style == "ov" else ""
    by = f"_by_{args.split}" if args.split else ""
    out = Path(args.out).expanduser() if args.out else (
        DESKTOP / f"GSE193688_UMAP_{args.color}{by}_{args.step}{suffix}.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    print(f"wrote {out}")
    print(f"{len(frame):,} cells, {len(values)} {args.color} categories, "
          f"{len(splits)} panel(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
