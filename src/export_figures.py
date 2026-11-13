#!/usr/bin/env python3
"""Export every dashboard figure as a PNG, one folder per pipeline section.

The Streamlit app only draws figures in memory. This script imports the app module,
replaces its `fig_show` helper with a `savefig` call, and runs each page once per
cohort. Widget selections are driven by a small shim so that pages with a dropdown
(UMAP colour, scaling analysis/metric, subtype cell type/modality) emit one figure
per option instead of only the first.

    python src/export_figures.py --out ~/Desktop/GSE193688_figures
    python src/export_figures.py --cohorts n24 --sections 09_utility
    python src/export_figures.py --list
"""

from __future__ import annotations

import argparse
import contextlib
import io
import logging
import re
import sys
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SCALE = ROOT / "results" / "scale"

# Streamlit logs a "missing ScriptRunContext" warning for every call in bare
# mode. The pages work, so keep the export output readable.
logging.getLogger("streamlit").setLevel(logging.CRITICAL)
logging.getLogger(
    "streamlit.runtime.scriptrunner_utils.script_run_context"
).setLevel(logging.CRITICAL)
warnings.filterwarnings("ignore")

sys.path.insert(0, str(ROOT / "dashboard"))
with contextlib.redirect_stderr(io.StringIO()):
    import app  # noqa: E402

COHORTS = [
    ("n24", SCALE / "n24", "n24 — all 24 donors"),
    ("n16", SCALE / "n16", "n16 — 16 donors"),
    ("n08", SCALE / "n08", "n08 — 8 donors"),
    ("n02", SCALE / "n02", "n02 — 2 donors"),
    ("pilot", ROOT / "results", "pilot — 2 donors (manual annotation)"),
]

SECTIONS = [
    ("01_overview", "Overview", "page_overview"),
    ("02_qc", "QC & composition", "page_qc"),
    ("03_annotation", "Annotation", "page_annotation"),
    ("04_fidelity", "Fidelity", "page_fidelity"),
    ("05_hypothesis", "Hypotheses", "page_hypothesis"),
    ("06_reversal", "Reversal", "page_reversal"),
    ("08_impact", "Cell-type impact", "page_impact"),
    ("09_utility", "Synthetic utility", "page_utility"),
    ("10_subtype", "Subtype transfer", "page_subtype"),
]


def slug(text: str) -> str:
    text = str(text).strip().lower()
    text = text.replace("→", " to ").replace("Δ", "delta").replace("×", "x")
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")[:60]


class _Column:
    """Forward everything to the real column, but route dropdowns to the shim."""

    def __init__(self, parent: "_Shim", real):
        object.__setattr__(self, "_parent", parent)
        object.__setattr__(self, "_real", real)

    def __getattr__(self, name):
        return getattr(self._real, name)

    def __enter__(self):
        self._real.__enter__()
        return self

    def __exit__(self, *exc):
        return self._real.__exit__(*exc)

    def selectbox(self, label, options, *args, **kwargs):
        return self._parent.selectbox(label, options, *args, **kwargs)


class _Shim:
    """Stand-in for the `streamlit` module inside the app.

    The app calls `st.subheader` right before drawing a figure, so the heading is
    used to name the file. `selectbox` consumes a per-run queue, which lets one
    page be called repeatedly with different dropdown choices.
    """

    def __init__(self, real):
        object.__setattr__(self, "_real", real)
        object.__setattr__(self, "queue", [])
        object.__setattr__(self, "selections", [])
        object.__setattr__(self, "heading", "figure")

    def __getattr__(self, name):
        return getattr(self._real, name)

    def title(self, text, *args, **kwargs):
        self.heading = text

    def subheader(self, text, *args, **kwargs):
        self.heading = text

    def columns(self, spec, **kwargs):
        return [_Column(self, r) for r in self._real.columns(spec, **kwargs)]

    def selectbox(self, label, options, *args, **kwargs):
        if self.queue:
            value = self.queue.pop(0)
        elif isinstance(options, dict):
            value = next(iter(options))
        else:
            value = list(options)[0]
        self.selections.append(value)
        return value


def _scaling_rounds() -> list[list[str]]:
    comp = ROOT / "results" / "scale" / "comparison.csv"
    if not comp.exists():
        return [[]]
    frame = pd.read_csv(comp, usecols=["analysis", "metric"])
    pairs = frame.drop_duplicates().sort_values(["analysis", "metric"])
    return [[a, m] for a, m in pairs.itertuples(index=False)]


def _subtype_rounds(root: Path) -> list[list[str]]:
    base = root / "subtype"
    if not base.is_dir():
        return [[]]
    celltypes = sorted(p.name for p in base.iterdir()
                        if p.is_dir() and (p / "rna").is_dir())
    if not celltypes:
        return [[]]
    rounds = []
    for ct in celltypes:
        mods = sorted(p.name for p in (base / ct).iterdir() if p.is_dir())
        for mod in mods:
            rounds.append([ct, mod])
    return rounds


def _rounds(section: str, cohort_dir: Path) -> list[list[str]]:
    if section == "02_qc":
        return [["celltype"], ["donor"], ["group"], ["cluster"]]
    if section == "07_scaling":
        return _scaling_rounds()
    if section == "10_subtype":
        return _subtype_rounds(cohort_dir)
    return [[]]


def export_section(page: str, outdir: Path, rounds: list[list[str]],
                   cohort: str, section: str, dpi: int,
                   manifest: list[dict]) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    func = getattr(app, page)
    written = 0
    state = {"n": 0, "seen": {}}

    for values in rounds:
        shim = _Shim(app.st)
        shim.queue = list(values)
        app.st = shim
        round_keys: dict[tuple, int] = {}

        def save(fig, _shim=shim, _state=state, _round=round_keys):
            heading = slug(_shim.heading) or "figure"
            key = (heading, tuple(str(v) for v in _shim.selections))
            if key in _state["seen"]:
                # The page redraws the same artifact once per dropdown round;
                # keep the first copy and drop the repeats.
                plt.close(fig)
                return None
            _round[key] = _round.get(key, 0) + 1
            _state["n"] += 1
            selections = [slug(v) for v in _shim.selections]
            if section == "07_scaling" and len(selections) == 2:
                # The scaling page sets a "Planned steps" heading before it
                # plots; use the analysis and metric as the name instead.
                heading = "_".join(selections)
                selections = []
            parts = [heading]
            if selections:
                parts.append("-".join(selections))
            if _round[key] > 1:
                parts.append(str(_round[key]))
            name = f"{_state['n']:02d}_" + "_".join(parts) + ".png"
            path = outdir / name
            fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
            width, height = fig.get_size_inches()
            plt.close(fig)
            manifest.append({
                "cohort": cohort,
                "section": section,
                "file": str(path),
                "heading": _shim.heading,
                "selection": " | ".join(str(v) for v in _shim.selections),
                "width_in": round(float(width), 2),
                "height_in": round(float(height), 2),
            })
            return path

        app.fig_show = save
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                func()
        except Exception:
            (outdir / "_ERROR.txt").write_text(traceback.format_exc())
        finally:
            plt.close("all")

        state["seen"].update(round_keys)
        written = state["n"]

    if written == 0 and not (outdir / "_MISSING.txt").exists():
        (outdir / "_MISSING.txt").write_text(
            "No figures were produced. The artifacts for this stage are "
            "missing for this cohort, or the page returned before plotting.\n"
        )
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(Path.home() / "Desktop" / "GSE193688_figures"))
    ap.add_argument("--cohorts", default="n24,n16,n08,n02,pilot",
                    help="comma-separated: n24,n16,n08,n02,pilot")
    ap.add_argument("--sections", default="",
                    help="comma-separated folder names or prefixes; default all")
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--list", action="store_true", help="list sections and exit")
    args = ap.parse_args()

    if args.list:
        for name, label, page in SECTIONS:
            print(f"{name:15s} {label:24s} {page}")
        print(f"{'07_scaling':15s} {'Sample-size scaling':24s} page_scaling")
        return 0

    out = Path(args.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)

    wanted_cohorts = {c.strip() for c in args.cohorts.split(",") if c.strip()}
    wanted_sections = {s.strip() for s in args.sections.split(",") if s.strip()}

    manifest: list[dict] = []
    total = 0

    for name, label, page in SECTIONS:
        if wanted_sections and name not in wanted_sections and label not in wanted_sections:
            continue
        for cohort, cohort_dir, cohort_label in COHORTS:
            if cohort not in wanted_cohorts:
                continue
            app.DATA = cohort_dir
            app.COHORT_LABEL = cohort_label
            section_dir = out / name / cohort
            count = export_section(page, section_dir, _rounds(name, cohort_dir),
                                  cohort, name, args.dpi, manifest)
            total += count
            print(f"{name}/{cohort}: {count} figure(s)")

    # Scaling is cohort-independent: the page itself plots every step.
    if not wanted_sections or "07_scaling" in wanted_sections \
            or "Sample-size scaling" in wanted_sections:
        app.DATA = SCALE / "n24"
        app.COHORT_LABEL = "all steps"
        count = export_section("page_scaling", out / "07_scaling",
                               _scaling_rounds(), "all", "07_scaling",
                               args.dpi, manifest)
        total += count
        print(f"07_scaling: {count} figure(s)")

    if manifest:
        frame = pd.DataFrame(manifest)
        frame.to_csv(out / "index.csv", index=False)

    readme = [
        "# GSE193688 figures",
        "",
        "Generated by `src/export_figures.py` from the dashboard artifacts.",
        f"Format: PNG, {args.dpi} dpi, white background.",
        "",
        "| Folder | What it shows |",
        "|---|---|",
        "| `01_overview` | Cohort summary and cell-type composition |",
        "| `02_qc` | Per-sample QC and UMAP coloured four ways |",
        "| `03_annotation` | Cluster sizes and donor mix |",
        "| `04_fidelity` | Real-vs-synthetic discriminator and label transfer |",
        "| `05_hypothesis` | Power on alternative data |",
        "| `06_reversal` | Model effect vs classifier contribution |",
        "| `07_scaling` | Every scaling metric against donor count |",
        "| `08_impact` | Cell-type impact and separability vs N |",
        "| `09_utility` | Synthetic-to-real learning curve |",
        "| `10_subtype` | Subtype transfer matrices, controls, Jaccard |",
        "",
        "Each section has one subfolder per cohort (`n24`, `n16`, `n08`, `n02`, `pilot`).",
        "`07_scaling` is cohort-independent and sits at the top level.",
        "`index.csv` lists every file with its heading and dropdown selection.",
        "A `_MISSING.txt` means the stage's artifacts do not exist for that cohort;",
        "`_ERROR.txt` records a page that failed to render.",
        "",
    ]
    (out / "README.md").write_text("\n".join(readme))

    print(f"\n{total} figure(s) written to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
