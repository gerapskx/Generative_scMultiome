#!/usr/bin/env python3
"""Streamlit dashboard for the GSE193688 multiome pipeline.

Eleven pages, each reading one stage's artifacts from `results/`. Missing
artifacts are reported instead of crashing, so the app is usable while the
pipeline is still running.

    streamlit run dashboard/app.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st
import yaml

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCALE = RESULTS / "scale"

# The utility figures are shared with `src/plot_utility.py` so the dashboard and
# the exported PNGs can never drift apart.
sys.path.insert(0, str(ROOT / "src"))
import plot_utility as pu  # noqa: E402

# Every page below reads from one cohort directory. The largest step (n24) is the
# default so the dashboard shows all samples, not just the 2-donor pilot; the pilot
# and the intermediate steps stay selectable for comparison.
STEP_DIRS = {
    "n24 — all 24 donors": SCALE / "n24",
    "n16 — 16 donors": SCALE / "n16",
    "n08 — 8 donors": SCALE / "n08",
    "n02 — 2 donors": SCALE / "n02",
    "pilot — 2 donors (manual annotation)": RESULTS,
}

st.set_page_config(page_title="GSE193688 multiome", layout="wide")

CELLTYPE_COLORS = {
    "Oligodendrocyte": "#4C72B0",
    "Astrocyte": "#DD8452",
    "Microglia": "#55A868",
    "OPC": "#C44E52",
}


# --------------------------------------------------------------------- utils --
@st.cache_data(show_spinner=False)
def read_csv(path: str) -> pd.DataFrame | None:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        return pd.read_csv(p, keep_default_na=False)
    except Exception:
        return None


@st.cache_data(show_spinner=False)
def read_json(path: str):
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return None
    try:
        return json.loads(p.read_text())
    except Exception:
        return None


def available(path: str) -> bool:
    p = Path(path)
    return p.exists() and p.stat().st_size > 0


def missing_box(what: str, stage: str) -> None:
    st.warning(f"{what} is not available yet. Run the `{stage}` stage first.")


def fig_show(fig) -> None:
    st.pyplot(fig, width="stretch")
    plt.close(fig)


def num(value, digits: int = 3, signed: bool = False) -> str:
    """Format a number for display, with an em dash for missing values."""
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "\u2014"
    if not np.isfinite(value):
        return "\u2014"
    return f"{value:+.{digits}f}" if signed else f"{value:.{digits}f}"


def obs() -> pd.DataFrame | None:
    return read_csv(str(d("cluster", "obs_annotated.csv")))


def d(*parts: str) -> Path:
    """Resolve a stage artifact inside the selected cohort directory."""
    return DATA.joinpath(*parts)


def n_lines(path: Path) -> int:
    if not path.exists():
        return 0
    with path.open() as fh:
        return sum(1 for _ in fh)


def config() -> dict:
    p = ROOT / "config" / "config.yaml"
    return yaml.safe_load(p.read_text()) if p.exists() else {}


# -------------------------------------------------------------------- pages --
def page_overview() -> None:
    st.title("GSE193688 — substantia nigra multiome")
    o = obs()
    is_pilot = DATA == RESULTS
    if is_pilot:
        st.caption(
            "snRNA + snATAC, CellRanger ARC 1.0.1, GRCh38. Pilot: C16007 (Aged) "
            "and P02284 (PD), one donor per condition."
        )
    else:
        n_don = o["donor"].nunique() if o is not None else "—"
        n_grp = o["group"].nunique() if o is not None else "—"
        st.caption(
            f"snRNA + snATAC, CellRanger ARC 1.0.1, GRCh38. {COHORT_LABEL}: "
            f"{n_don} donors across {n_grp} conditions. Every stage below is "
            "computed on this cohort."
        )

    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Nuclei", f"{len(o):,}" if o is not None else "—")
    c2.metric("Donors", o["donor"].nunique() if o is not None else "—")
    c3.metric("Cell types", o["celltype"].nunique() if o is not None else "—")
    genes = n_lines(d("ingest", "rna_features.tsv"))
    peaks = n_lines(d("ingest", "atac_features.tsv"))
    c4.metric("Genes", f"{genes:,}" if genes else "—")
    c5.metric("Peaks", f"{peaks:,}" if peaks else "—")

    st.subheader("Pipeline status")
    artifacts = {
        "0 ingest": d("ingest", "obs.csv"),
        "1 assemble": (ROOT / "data" / "processed" / "rna.h5ad") if is_pilot
        else d("assemble", "rna.h5ad"),
        "2 cluster": d("cluster", "clusters.csv"),
        "3 annotation": d("cluster", "obs_annotated.csv"),
        "4 split": d("split", "split.csv"),
        "5 features": d("features", "rna_features.tsv"),
        "6 simulation": d("simulation", "manifest.csv"),
        "7 fidelity": d("fidelity", "summary.json"),
        "8 hypothesis": d("hypothesis", "summary.json"),
        "9 reversal": d("reversal", "summary.json"),
        "10 impact": d("impact", "summary.json"),
        "11 utility": d("utility", "summary.json"),
        "12 scaling": SCALE / "comparison.csv",
    }
    status = pd.DataFrame(
        [{"stage": k, "state": "done" if available(str(v)) else "pending",
          "artifact": str(v.relative_to(ROOT))} for k, v in artifacts.items()]
    )
    st.dataframe(status, width="stretch", hide_index=True)

    if o is not None:
        st.subheader("Cell-type composition")
        comp = (o.groupby(["celltype", "group"]).size().rename("nuclei")
                .reset_index())
        fig, ax = plt.subplots(figsize=(7, 3.2))
        sns.barplot(data=comp, x="celltype", y="nuclei", hue="group", ax=ax)
        ax.set_xlabel("")
        ax.set_ylabel("nuclei")
        ax.tick_params(axis="x", rotation=20)
        fig_show(fig)

    if is_pilot:
        st.error(
            "**Interpretation limits.** The pilot is oligodendrocyte-dominated "
            "(~88% of nuclei) and dopaminergic neurons are essentially absent "
            "(`TH` has 0 UMIs in C16007). Aged vs PD is also a donor contrast "
            "because there is one donor per condition. The reversal page is a "
            "donor contrast, not a replicated condition effect."
        )
    else:
        st.info(
            "**Full cohort.** This view is computed on all selected donors. The "
            "n02 step is still a donor contrast (one donor per condition); from "
            "n08 the donor-held-out metrics become estimable. Use the cohort "
            "selector in the sidebar to compare steps."
        )


def page_report() -> None:
    st.title("Report — what each section shows and what it means")
    st.caption(
        f"Cohort: **{COHORT_LABEL}**. Every number below is read from this "
        "cohort's artifacts, so changing the sidebar selector regenerates the report."
    )

    o = obs()
    qc = read_csv(str(d("ingest", "qc_per_sample.csv")))
    disc = read_csv(str(d("fidelity", "discriminator.csv")))
    calib = read_csv(str(d("hypothesis", "calibration.csv")))
    impact = read_csv(str(d("impact", "celltype_impact.csv")))
    curve = read_csv(str(d("utility", "learning_curve.csv")))
    poc = read_json(str(d("reversal", "poc", "summary.json")))
    comp = read_csv(str(SCALE / "comparison.csv"))
    sub_transfer = read_csv(str(d("subtype", "transfer", "transfer.csv")))
    sub_jac = read_csv(str(d("subtype", "transfer", "jaccard.csv")))
    sub_comp = read_csv(str(d("subtype", "comparison.csv")))
    if sub_transfer is not None and not sub_transfer.empty:
        sub_transfer = sub_transfer.copy()
        sub_transfer["value"] = pd.to_numeric(
            sub_transfer["value"], errors="coerce")

    def num(x, nd: int = 3) -> str:
        try:
            v = float(x)
        except (TypeError, ValueError):
            return "—"
        return f"{v:.{nd}f}" if np.isfinite(v) else "—"

    def series(analysis: str, metric: str, **filters) -> dict:
        if comp is None:
            return {}
        s = comp[(comp["analysis"] == analysis) & (comp["metric"] == metric)].copy()
        for k, v in filters.items():
            if k == "synthetic_size":
                s = s[pd.to_numeric(s[k], errors="coerce") == float(v)]
            elif k in s.columns:
                s = s[s[k].astype(str) == str(v)]
        s["value"] = pd.to_numeric(s["value"], errors="coerce")
        return s.groupby("N")["value"].mean().to_dict()

    def sub_auc(ct: str, mod: str, train: str, test: str,
                model: str = "logistic_regression"):
        if sub_transfer is None:
            return None
        g = sub_transfer[(sub_transfer["celltype"] == ct) &
                         (sub_transfer["modality"] == mod) &
                         (sub_transfer["train_domain"] == train) &
                         (sub_transfer["test_domain"] == test) &
                         (sub_transfer["model"] == model) &
                         (sub_transfer["metric"] == "auc") &
                         (sub_transfer["stratum"] == "overall")]
        return g["value"].mean() if not g.empty else None

    blocks: list[str] = []
    add = blocks.append

    # ---------------------------------------------------------------- study --
    add(
        "## 1. What the study asks\n\n"
        "The cohort is **Aged vs PD** substantia nigra (snRNA + snATAC). The pilot "
        "has one donor per condition, so every pilot Aged-vs-PD number is also a "
        "donor contrast. The scaling series repeats **every stage** on nested donor "
        "sets — n02, n08, n16, n24 — with a fixed n24 feature set and five "
        "simulation seeds. The question is not only whether a p-value survives, but "
        "**how the interpretation changes with N**."
    )

    # ------------------------------------------------------------------- QC --
    if o is not None:
        dbl = "—"
        if qc is not None and {"doublets", "cells_in"}.issubset(qc.columns):
            denom = pd.to_numeric(qc["cells_in"], errors="coerce").sum()
            if denom:
                dbl = f"{100 * pd.to_numeric(qc['doublets'], errors='coerce').sum() / denom:.1f}%"
        top = o["celltype"].value_counts()
        top_txt = ", ".join(
            f"{ct} {100 * n / len(o):.0f}%" for ct, n in top.head(3).items()
        )
        add(
            f"## 2. QC and composition\n\n"
            f"**What it shows.** Cells and donors that pass QC, per-sample doublet "
            f"rate, and the cell-type composition.\n\n"
            f"**This cohort.** {len(o):,} nuclei, {o['donor'].nunique()} donors, "
            f"{o['celltype'].nunique()} cell types, doublet rate {dbl}. Composition: "
            f"{top_txt}.\n\n"
            f"**How to read it.** A healthy step has a low doublet rate and a "
            f"composition consistent with the previous step. A large shift in "
            f"composition between steps means the donor set, not biology, is driving "
            f"the downstream contrast.\n\n"
            f"**Caveat.** The cohort is oligodendrocyte-dominated and dopaminergic "
            f"neurons are essentially absent (`TH` ≈ 0), so any PD-neurodegeneration "
            f"narrative is not supported by these samples."
        )

    # ----------------------------------------------------------- annotation --
    if o is not None:
        add(
            "## 3. Clustering and annotation\n\n"
            "**What it shows.** Harmony clusters and their assignment to cell types "
            "from a canonical marker panel (auto annotation; the pilot used a manual gate).\n\n"
            "**How to read it.** Check that cluster sizes and top markers are "
            "consistent with the assigned cell type, and that the map is stable across "
            "steps. A label that appears only at large N is usually a small cluster "
            "gaining enough cells to score.\n\n"
            "**Caveat.** Annotation is marker-score based, not validated by an "
            "independent assay; rare types (Neuron, Pericyte, T cell) are small."
        )

    # --------------------------------------------------------------- split --
    add(
        "## 4. Split and features\n\n"
        "**What it shows.** A cell-level stratified fit/holdout split "
        "(cell type × donor) and the 2,000 RNA genes / 2,000 ATAC peaks used "
        "downstream. Every step reuses the n24 selection so N is not confounded with "
        "feature drift; each step's own selection is kept for the stability metric.\n\n"
        "**How to read it.** Jaccard and log2FC Spearman vs n24 measure whether the "
        "feature set has converged. Low overlap means the smaller step is selecting a "
        "different biology, not just fewer features."
    )

    # ------------------------------------------------------------- fidelity --
    if disc is not None:
        disc["auc"] = pd.to_numeric(disc["auc"], errors="coerce")
        rows = []
        for mod in ("rna", "atac"):
            for model in ("logistic_regression", "random_forest"):
                s = disc[(disc["modality"] == mod)
                          & (disc["comparison"] == "real_vs_synthetic")
                          & (disc["model"] == model)]["auc"]
                rows.append(f"{mod.upper()} {model.split('_')[0]} {num(s.mean())}")
        add(
            "## 5. Simulation fidelity\n\n"
            "**What it measures.** A classifier trained to separate real holdout "
            "cells from scDesign3 synthetic cells. AUC 0.5 means the simulation is "
            "indistinguishable from real data; further from 0.5 means more separable.\n\n"
            f"**This cohort.** Real vs synthetic AUC: {'; '.join(rows)}. "
            "The `identity` control sits at 0.5 by construction and "
            "`positive_control` (Aged vs PD) shows the test has power.\n\n"
            "**How to read it.** Lower is better here. A rising AUC with N is a "
            "negative result: the simulator becomes easier to detect as the real cohort "
            "grows, so it is not a sample-size-scalable substitute for real cells.\n\n"
            "**Caveat.** AUC distance is not a biological fidelity score; it only "
            "says whether a discriminative model can tell the domains apart."
        )

    # ----------------------------------------------------------- hypothesis --
    if calib is not None:
        for col in calib.columns:
            if col.startswith(("power_", "fdp_", "n_sig_")):
                calib[col] = pd.to_numeric(calib[col], errors="coerce")
        cond = calib[calib["contrast_type"] == "condition"]
        bits = []
        for mod in ("rna", "atac"):
            for sett in ("alt", "null"):
                g = cond[(cond["modality"] == mod) & (cond["set"] == sett)]
                if g.empty:
                    continue
                bits.append(
                    f"{mod.upper()} {sett}: power {num(g['power_0.05'].mean())}, "
                    f"FDP {num(g['fdp_0.05'].mean())}, "
                    f"discoveries {num(g['n_sig_0.05'].mean(), 1)}"
                )
        add(
            "## 6. Null / alternative calibration\n\n"
            "**What it measures.** Differential expression is run on the `alt` "
            "simulation (which carries the fitted condition effect) and the `null` "
            "simulation (condition coefficients zeroed). Ground truth is the same "
            "test on the real fit cells (FDR 0.05, |log2FC| ≥ 0.5).\n\n"
            f"**This cohort (FDR 0.05).** {'; '.join(bits)}.\n\n"
            "**How to read it.** The `null` side must stay silent — near-zero "
            "discoveries and near-zero FC correlation. The `alt` side is power and "
            "FDP against the real calls: power should rise and FDP should fall if "
            "the simulation reproduces the real effect.\n\n"
            "**Caveat.** FDP is undefined when a simulation makes no "
            "discoveries, so a null FDP of 1.0 with zero discoveries is a "
            "degenerate value, not a failure."
        )

    # ------------------------------------------------------------- reversal --
    if poc:
        rec = poc.get("recovery", {})
        con = poc.get("donor_confound", {})
        add(
            "## 7. Reversal proof of concept\n\n"
            "**What it measures.** Three truth-known checks: (i) does the `alt` "
            "simulation reproduce the real PD − Aged effect ranking "
            "(Spearman ρ, top-100 Jaccard)? (ii) is the `null` silent "
            "(ρ ≈ 0, few effects above the real q95)? (iii) how confounded are "
            "condition and donor (donor separability, donor-grouped condition "
            "accuracy, donor-level permutation)?\n\n"
            f"**This cohort.** Recovery ρ {num(rec.get('spearman_real_alt'))}, "
            f"top-100 Jaccard {num(rec.get('topn_jaccard_alt'))}; null ρ "
            f"{num(rec.get('spearman_real_null'))}, null effects above real q95 "
            f"{num(rec.get('frac_null_above_real_q95'))}; donor balanced "
            f"accuracy {num(con.get('donor_balanced_accuracy'))}, "
            f"donor-grouped condition accuracy "
            f"{num(con.get('condition_balanced_accuracy'))} "
            f"({con.get('identifiable_strata', '—')}/{con.get('strata', '—')} "
            f"strata identifiable).\n\n"
            "**How to read it.** High recovery and near-zero null correlation "
            "mean the analysis is wired correctly. Condition accuracy near 0.5 on "
            "unseen donors means the condition signal does **not** transfer across "
            "donors; donor accuracy near 1.0 means donor identity dominates.\n\n"
            "**Caveat.** At n02 condition and donor are the same axis, so the "
            "apparently high recovery is the donor contrast being replayed and the "
            "condition metric is not identifiable."
        )

    # --------------------------------------------------------------- impact --
    if impact is not None:
        impact["de_burden"] = pd.to_numeric(impact["de_burden"], errors="coerce")
        impact["donor_auc"] = pd.to_numeric(impact["donor_auc"], errors="coerce")
        lines = []
        for mod in ("rna", "atac"):
            g = impact[impact["modality"] == mod].dropna(subset=["de_burden"])
            if g.empty:
                continue
            top = g.nlargest(2, "de_burden")
            lines.append(
                f"{mod.upper()}: " + ", ".join(
                    f"{r['celltype']} (burden {num(r['de_burden'])}, "
                    f"donor AUC {num(r['donor_auc'])})"
                    for _, r in top.iterrows()
                )
            )
        add(
            "## 8. Cell-type impact\n\n"
            "**What it measures.** Which cell types the condition affects most: "
            "the share of testable features that are DE, mean |log2FC|, a "
            "donor-held-out condition classifier, and the variance explained by "
            "condition vs donor.\n\n"
            f"**This cohort.** {' | '.join(lines)}.\n\n"
            "**How to read it.** High DE burden = many features move. Donor-held-out "
            "AUC near 0.5 = the condition signal does not generalise to unseen "
            "donors, so the burden is donor-driven. `eta2_donor` >> `eta2_condition` "
            "is the same warning.\n\n"
            "**Caveat.** At n02 `donor_auc` is undefined because condition == "
            "donor; the metric only becomes a condition test from n08."
        )

    # -------------------------------------------------------------- utility --
    if curve is not None:
        curve["auc_mean"] = pd.to_numeric(curve["auc_mean"], errors="coerce")
        curve["synthetic_size"] = pd.to_numeric(curve["synthetic_size"],
                                                 errors="coerce")
        lines = []
        for mod in ("rna", "atac"):
            g = curve[curve["modality"] == mod]
            real = g[g["method"] == "real_trained"]["auc_mean"].mean()
            syn = g[g["method"] == "synthetic"]
            if syn.empty:
                continue
            syn = syn.dropna(subset=["auc_mean"])
            if syn.empty:
                continue
            big = syn.loc[syn["synthetic_size"].idxmax()]
            lines.append(
                f"{mod.upper()}: real-trained {num(real)}, "
                f"synthetic {int(big['synthetic_size']):,} cells {num(big['auc_mean'])}"
            )

        # Full learning curve, not just the endpoint.
        curve["bal_mean"] = pd.to_numeric(curve.get("bal_mean"), errors="coerce")
        order = {"real_trained": 0, "synthetic": 1, "naive_poisson": 2,
                 "permuted_labels": 3, "majority": 4}
        table = curve.copy()
        table["_o"] = table["method"].map(order).fillna(9)
        table = table.sort_values(["modality", "_o", "synthetic_size"])
        util_table = ("| modality | method | synthetic cells | AUC | balanced acc. |\n"
                      "|---|---|---|---|---|\n")
        for _, r in table.iterrows():
            util_table += (
                f"| {r['modality'].upper()} | {r['method']} | "
                f"{int(r['synthetic_size']):,} | {num(r['auc_mean'])} | "
                f"{num(r['bal_mean'])} |\n"
            )

        add(
            "## 9. Synthetic-to-real utility (ML PoC)\n\n"
            "**What it measures.** A condition classifier is trained on synthetic "
            "cells generated from training donors and tested on **held-out real "
            "donors**. The learning curve varies synthetic training size from 500 to "
            "10,000.\n\n"
            f"**This cohort.** {'; '.join(lines)}.\n\n"
            + util_table +
            "\n**How to read it.** `real_trained` is the ceiling, `naive_poisson` the "
            "marginal-only lower bound, and `permuted_labels`/`majority` are chance. "
            "A synthetic curve that rises with size but stays below the real ceiling "
            "means the simulation carries some condition signal but not all of it.\n\n"
            "**Caveat.** The ceiling is a single donor-held-out fold, so it moves "
            "with the composition of the held-out donors; compare it across steps with "
            "that in mind. n02 is `not_estimable` (one donor per condition)."
        )

    # -------------------------------------------------------------- scaling --
    if comp is not None:
        rec = series("reversal_poc", "spearman_real_alt")
        donor = series("reversal_poc", "donor_balanced_accuracy")
        cond = series("reversal_poc", "condition_balanced_accuracy")
        fid = series("fidelity", "auc_logistic_regression", modality="rna")
        util = series("utility", "auc_mean", label="synthetic",
                      synthetic_size=10000, modality="rna")
        stab = series("stability", "jaccard_top100_vs_fixed")
        ns = sorted(rec)
        if ns:
            head = "| N | recovery ρ | donor bal. acc. | condition bal. acc. | RNA fid. AUC | RNA syn-10k AUC | feature Jaccard |\n"
            head += "|---|---|---|---|---|---|---|\n"
            for n in ns:
                # the fixed set is n24, so its self-overlap is 1 by definition
                jac = stab.get(n, 1.0 if n == max(ns) else None)
                head += (
                    f"| {n} | {num(rec.get(n))} | {num(donor.get(n))} | "
                    f"{num(cond.get(n))} | {num(fid.get(n))} | "
                    f"{num(util.get(n))} | {num(jac)} |\n"
                )
            add(
                "## 10. Sample-size scaling\n\n"
                "**What it shows.** Every metric above recomputed at n02, n08, n16 "
                "and n24, with the feature set fixed to n24.\n\n"
                + head +
                "\n**How to read it.** Metrics that improve with N are gaining "
                "donor-level information. A recovery ρ that falls with N, condition "
                "accuracy stuck near 0.5 and donor accuracy falling toward 0.5 all say "
                "the same thing: the apparent pilot signal was donor identity, and "
                "adding donors makes the condition effect harder — but more honest — to "
                "estimate. Rising feature Jaccard means the selection is converging."
            )

    # ------------------------------------------------------ subtype transfer --
    if sub_transfer is not None and not sub_transfer.empty:
        tr = sub_transfer.copy()
        tr["value"] = pd.to_numeric(tr["value"], errors="coerce")
        tr = tr[(tr["metric"] == "auc") & (tr["stratum"] == "overall")]

        def sub_auc(ct, mod, train, test, model="logistic_regression"):
            g = tr[(tr["celltype"] == ct) & (tr["modality"] == mod) &
                   (tr["train_domain"] == train) & (tr["test_domain"] == test) &
                   (tr["model"] == model)]
            return g["value"].mean() if not g.empty else None

        ari = {}
        if sub_comp is not None and not sub_comp.empty:
            for _, r in sub_comp[sub_comp["analysis"] == "concordance"].iterrows():
                if r["metric"] == "rna_atac_ari":
                    ari[str(r["celltype"])] = r["value"]
        jac = {}
        if sub_jac is not None and not sub_jac.empty:
            for _, r in sub_jac.iterrows():
                jac[(str(r["celltype"]), str(r["modality"]))] = r["jaccard"]

        head = ("| cell type | modality | real→real | syn alt→real | "
                "real→syn alt | syn null→real | null→null | "
                "var-feature Jaccard | RNA/ATAC ARI |\n")
        head += "|---|---|---|---|---|---|---|---|\n"
        for ct in sorted(tr["celltype"].unique()):
            for mod in sorted(tr[tr["celltype"] == ct]["modality"].unique()):
                head += (
                    f"| {ct} | {mod} | "
                    f"{num(sub_auc(ct, mod, 'real', 'real'))} | "
                    f"{num(sub_auc(ct, mod, 'syn_alt', 'real'))} | "
                    f"{num(sub_auc(ct, mod, 'real', 'syn_alt'))} | "
                    f"{num(sub_auc(ct, mod, 'syn_null', 'real'))} | "
                    f"{num(sub_auc(ct, mod, 'syn_null', 'syn_null'))} | "
                    f"{num(jac.get((ct, mod)))} | "
                    f"{num(ari.get(ct))} |\n"
                )
        add(
            "## 11. Subtype transfer\n\n"
            "**What it shows.** The two disease-impacted cell types "
            "(oligodendrocyte, microglia) are subclustered into k=2 on training "
            "donors only, separately in RNA and ATAC. scDesign3 then simulates "
            "balanced 2×2 subtype × condition data in three modes: **alt** "
            "(`subtype + condition`), **null** (subtype collapsed, `condition` "
            "refit) and **null_zero** (alt fit with subtype coefficients zeroed). "
            "A logistic classifier is trained in one domain and tested in another, "
            "with donors held out on the real side.\n\n"
            + head +
            "\n**How to read it.** `real→real` is the donor-disjoint ceiling. "
            "`syn alt→real` tests whether the simulator preserves the subtype "
            "axis well enough for a real-trained problem; `real→syn alt` is the "
            "reverse direction. A large positive asymmetry means the synthetic "
            "cells are easier to classify than real cells. The null simulations "
            "should collapse to chance when trained on synthetic and tested on "
            "real. `var-feature Jaccard` is the overlap between the subtype "
            "variable features and the fixed n24 feature set. `RNA/ATAC ARI` is "
            "the agreement of the two modality-specific subtype labels on shared "
            "cells; it is a ceiling on cross-modality transfer, not a result.\n\n"
            "**Caveat.** The Microglia split is weak in both modalities "
            "(silhouette ≈0.16–0.18) and RNA/ATAC agreement is near zero, so "
            "it is not a validated shared subtype. The Oligodendrocyte split is "
            "stronger in RNA than ATAC (silhouette 0.26 vs 0.14) and ARI ≈0.19. "
            "Treat these as unsupervised axes, not known developmental subtypes."
        )

    # ---------------------------------------------------------- limitations --
    add(
        "## 12. Limitations to keep in mind\n\n"
        "- **No dopaminergic neurons.** `TH` is essentially absent, so the "
        "cell-type labels are glial and the PD narrative is not supported.\n"
        "- **Oligodendrocyte dominance.** ~79–88% of nuclei are oligodendrocyte "
        "lineage, so power is concentrated there.\n"
        "- **Donor confounding.** At n02 condition and donor are the same axis; "
        "from n08 the donor-held-out metrics are the honest ones.\n"
        "- **Simulation is a fixed-effect approximation.** scDesign3 has no "
        "cell-level covariates, so its condition effect is group-level and its "
        "fidelity does not improve with N.\n"
        "- **Feature set fixed at n24.** Smaller steps reuse the n24 genes/peaks, "
        "so the curve isolates N rather than feature drift.\n"
        "- **One utility fold per step.** The real-trained ceiling moves with "
        "the held-out donor composition.\n"
        "- **Subtypes are unsupervised axes.** The k=2 splits are not validated "
        "developmental subtypes; Microglia RNA and ATAC labels disagree (ARI ≈0), "
        "so cross-modality transfer there is a negative result, not a hidden signal."
    )

    report = "\n\n".join(blocks)
    st.markdown(report)
    st.download_button(
        "Download this report (Markdown)",
        report,
        file_name=f"gse193688_report_{COHORT_LABEL.split(' ')[0]}.md",
        mime="text/markdown",
    )


def page_qc() -> None:
    st.title("QC and composition")
    qc = read_csv(str(d("ingest", "qc_per_sample.csv")))
    o = obs()

    if qc is not None:
        st.subheader("Per-sample QC")
        st.dataframe(qc, width="stretch", hide_index=True)
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.2))
        for ax, col in zip(axes, ["cells_in", "cells_out", "doublets"]):
            if col in qc.columns:
                sns.barplot(data=qc, x="sample", y=col, ax=ax, color="#4C72B0")
                ax.set_title(col)
                ax.set_xlabel("")
                ax.tick_params(axis="x", rotation=30)
        fig.tight_layout()
        fig_show(fig)

    if o is None:
        missing_box("Annotated observation table", "annotate")
        return

    st.subheader("Cell type x condition")
    tab = pd.crosstab(o["celltype"], o["group"])
    st.dataframe(tab, width="stretch")

    st.subheader("UMAP")
    if {"umap_1", "umap_2"}.issubset(o.columns):
        colour = st.selectbox("Colour by", ["celltype", "donor", "group", "cluster"])
        fig, ax = plt.subplots(figsize=(6.5, 5))
        if colour == "celltype":
            for ct, sub in o.groupby("celltype"):
                ax.scatter(sub["umap_1"], sub["umap_2"], s=5, alpha=0.7,
                           label=ct, color=CELLTYPE_COLORS.get(ct))
            ax.legend(markerscale=3, fontsize=8, frameon=False)
        else:
            for val, sub in o.groupby(colour):
                ax.scatter(sub["umap_1"], sub["umap_2"], s=5, alpha=0.6, label=str(val))
            ax.legend(markerscale=3, fontsize=8, frameon=False)
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
        fig_show(fig)


def page_annotation() -> None:
    st.title("Clustering and annotation")
    clusters = read_csv(str(d("cluster", "clusters.csv")))
    sizes = read_csv(str(d("cluster", "cluster_sizes.csv")))
    markers = read_csv(str(d("cluster", "cluster_markers.csv")))
    cmap = read_csv(str(d("cluster", "celltype_map.csv")))
    if cmap is None:
        cmap = read_csv(str(ROOT / "config" / "celltype_map.csv"))

    if clusters is None:
        missing_box("Cluster assignment", "cluster")
        return

    if cmap is not None:
        st.subheader("Cluster → cell type")
        st.dataframe(cmap, width="stretch", hide_index=True)
        if DATA == RESULTS:
            st.caption("Manual gate used for the pilot.")
        else:
            st.caption("Marker-panel auto annotation (`src/auto_annotate.py`).")

    if sizes is not None:
        st.subheader("Cluster sizes and donor mix")
        st.dataframe(sizes, width="stretch", hide_index=True)
        donor_cols = [c for c in sizes.columns if c not in ("cluster", "total")]
        if donor_cols:
            plot = sizes.set_index("cluster")[donor_cols]
            fig, ax = plt.subplots(figsize=(11, 3.2))
            plot.plot(kind="bar", stacked=True, ax=ax, width=0.85)
            ax.set_ylabel("nuclei")
            ax.set_xlabel("")
            ax.legend(frameon=False, fontsize=8)
            fig_show(fig)

    if markers is not None:
        st.subheader("Top markers per cluster")
        cl = st.selectbox("Cluster", sorted(markers["cluster"].unique()))
        sub = markers[(markers["cluster"] == cl) & (markers["p_val_adj"] < 0.05)]
        sub = sub.nlargest(15, "avg_log2FC")[["gene", "avg_log2FC", "pct.1", "pct.2"]]
        st.dataframe(sub, width="stretch", hide_index=True)
        st.caption(
            "Top-log2FC markers are often lncRNA/MT/ribosomal; the manual map "
            "above was assigned from canonical marker-panel scores, not this table."
        )


def page_fidelity() -> None:
    st.title("Simulation fidelity")
    disc = read_csv(str(d("fidelity", "discriminator.csv")))
    transfer = read_csv(str(d("fidelity", "label_transfer.csv")))
    per_class = read_csv(str(d("fidelity", "per_class.csv")))

    if disc is None:
        missing_box("Fidelity results", "fidelity")
        return

    st.subheader("Discriminator: real vs synthetic")
    agg = (disc.groupby(["modality", "comparison", "model"])["auc"]
           .agg(["mean", "std", "count"]).reset_index())
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.6), sharey=True)
    for ax, mod in zip(axes, ["rna", "atac"]):
        sub = agg[agg["modality"] == mod]
        if sub.empty:
            continue
        sns.barplot(data=sub, x="comparison", y="mean", hue="model", ax=ax)
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set_ylim(0, 1.02)
        ax.set_title(mod.upper())
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=25)
    fig_show(fig)
    st.caption(
        "A faithful simulator sits at AUC 0.5. `identity` is the null control, "
        "`positive_control` (Aged vs PD) shows the test has power, and "
        "`naive_poisson` shows how much a marginal-only simulator loses."
    )

    if transfer is not None:
        st.subheader("Label transfer (cell type)")
        agg = (transfer.groupby(["modality", "direction", "model"])["balanced_accuracy"]
               .agg(["mean", "std", "count"]).reset_index())
        fig, ax = plt.subplots(figsize=(9, 3.4))
        sns.barplot(data=agg, x="direction", y="mean", hue="model", ax=ax)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("balanced accuracy")
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)
        fig_show(fig)
        st.dataframe(agg.round(3), width="stretch", hide_index=True)

    if per_class is not None:
        st.subheader("Per-class recall (synthetic -> real)")
        sub = per_class[per_class["direction"] == "synthetic_to_real"]
        if not sub.empty:
            pivot = sub.pivot_table(index=["modality", "celltype"],
                                    columns="model", values="recall", aggfunc="mean")
            fig, ax = plt.subplots(figsize=(6, max(2, 0.4 * len(pivot))))
            sns.heatmap(pivot, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1, ax=ax)
            fig_show(fig)


def page_hypothesis() -> None:
    st.title("Null / alternative hypothesis benchmark")
    calib = read_csv(str(d("hypothesis", "calibration.csv")))
    truth = read_csv(str(d("hypothesis", "truth.csv")))

    if calib is None:
        missing_box("Hypothesis results", "hypothesis")
        return

    # Empty cells are kept as "" by read_csv; coerce before plotting/rounding.
    for col in calib.columns:
        if col.startswith(("power_", "fdp_", "n_sig_")):
            calib[col] = pd.to_numeric(calib[col], errors="coerce")

    fdr_cols = [c for c in calib.columns if c.startswith("power_")]
    fdr_vals = [float(c.split("_")[1]) for c in fdr_cols]

    st.subheader("Power on alternative data")
    alt = calib[calib["set"] == "alt"]
    if not alt.empty and fdr_cols:
        long = alt.melt(id_vars=["modality", "contrast_type", "contrast"],
                        value_vars=fdr_cols, var_name="fdr", value_name="power")
        long["fdr"] = long["fdr"].str.replace("power_", "").astype(float)
        agg = long.groupby(["modality", "contrast_type", "fdr"])["power"].mean().reset_index()
        fig, axes = plt.subplots(1, 2, figsize=(11, 3.4), sharey=True)
        for ax, mod in zip(axes, ["rna", "atac"]):
            sub = agg[agg["modality"] == mod]
            if sub.empty:
                continue
            sns.lineplot(data=sub, x="fdr", y="power", hue="contrast_type",
                         marker="o", ax=ax)
            ax.set_ylim(0, 1.02)
            ax.set_title(mod.upper())
            ax.set_xlabel("target FDR")
            ax.set_ylabel("mean power")
        fig_show(fig)

    st.subheader("Null calibration")
    null = calib[calib["set"] == "null"]
    if not null.empty and fdr_cols:
        cols = ["modality", "contrast_type", "contrast"] + fdr_cols
        st.dataframe(null[cols].round(3), width="stretch", hide_index=True)
        st.caption(
            "`condition` contrasts must not be called under `null` (zeroed "
            "condition coefficients). `celltype` contrasts should retain power "
            "under `null`, because the null preserves cell-type structure."
        )

    if truth is not None:
        st.subheader("Real-data reference calls")
        st.dataframe(truth, width="stretch", hide_index=True)


def page_reversal() -> None:
    st.title("Counterfactual PD → Aged reversal")
    combined = read_csv(str(d("reversal", "combined.csv")))
    top = read_csv(str(d("reversal", "top_features.csv")))
    model_summary = read_json(str(d("reversal", "model_summary.json")))
    summary = read_json(str(d("reversal", "summary.json")))

    if combined is None:
        missing_box("Reversal results", "reversal")
        return

    if DATA == RESULTS:
        st.error(
            "This is a **donor contrast**, not a replicated condition effect: one "
            "donor per condition. Treat the ranking as hypothesis-generating."
        )
    else:
        st.info(
            "Fitted on every donor in this cohort. The counterfactual is "
            "model-based, so read it as a condition effect estimated across "
            "donors rather than a per-donor experiment."
        )

    if summary:
        strata = pd.DataFrame(summary.get("strata", []))
        if not strata.empty:
            st.subheader("Agreement between the two ranking methods")
            st.dataframe(strata.round(3), width="stretch", hide_index=True)
            st.caption(
                "`spearman_model_vs_classifier` is the agreement between the "
                "scDesign3 counterfactual and the classifier contribution. "
                "`spearman_reverse_direction` should be negative; "
                "`spearman_permuted` should be near zero."
            )

    st.subheader("Model Δ vs classifier contribution")
    for (modality, celltype), sub in combined.groupby(["modality", "celltype"]):
        fig, ax = plt.subplots(figsize=(5.2, 4))
        ax.scatter(sub["log2fc"], sub["contribution"], s=8, alpha=0.5)
        ax.axhline(0, color="grey", lw=0.8)
        ax.axvline(0, color="grey", lw=0.8)
        ax.set_xlabel("model log2FC (PD → Aged)")
        ax.set_ylabel("classifier contribution (toward Aged)")
        ax.set_title(f"{modality.upper()} / {celltype}")
        fig_show(fig)

    if top is not None:
        st.subheader("Top features by method")
        c1, c2 = st.columns(2)
        with c1:
            method = st.selectbox("Method", sorted(top["method"].unique()))
        with c2:
            mod = st.selectbox("Modality", sorted(top["modality"].unique()))
        sub = top[(top["method"] == method) & (top["modality"] == mod)].head(25)
        cols = [c for c in ["rank", "celltype", "feature", "log2fc", "contribution"]
                if c in sub.columns]
        st.dataframe(sub[cols], width="stretch", hide_index=True)

    if model_summary:
        st.subheader("Model controls")
        strata = pd.DataFrame(model_summary.get("strata", []))
        if not strata.empty:
            keep = [c for c in ["modality", "celltype", "reverse_correlation",
                                "observed_sd", "permuted_sd", "permuted_q95",
                                "frac_observed_above_perm_q95"] if c in strata.columns]
            st.dataframe(strata[keep].round(3), width="stretch",
                         hide_index=True)
            st.caption(
                "`reverse_correlation` is −1 by construction. "
                "`frac_observed_above_perm_q95` is the share of features whose "
                "counterfactual effect exceeds the 95th percentile of the "
                "permuted-condition null."
            )


def page_scaling() -> None:
    st.title("Sample-size scaling")
    plan = read_csv(str(ROOT / "config" / "scaling" / "plan.csv"))
    comp = read_csv(str(SCALE / "comparison.csv"))

    st.caption(
        "Nested Aged/PD cohorts. Every step is a superset of the previous one, "
        "the pilot donors are always included, and the feature set is fixed to "
        "the largest step so N is not confounded with feature drift. "
        "This page always plots every step; the sidebar cohort selector does "
        "not filter it."
    )

    if plan is not None:
        st.subheader("Planned steps")
        st.dataframe(plan, width="stretch", hide_index=True)

    if comp is None:
        missing_box("Scaling comparison", "scaling_compare")
        return

    comp["value"] = pd.to_numeric(comp["value"], errors="coerce")
    analyses = sorted(comp["analysis"].unique())
    analysis = st.selectbox("Analysis", analyses)
    sub = comp[comp["analysis"] == analysis]

    if analysis == "utility":
        # The utility rows carry a synthetic-size axis; averaging it away (the
        # generic N-vs-value panel below) hides the learning curve entirely.
        util = sub.copy()
        util["N"] = pd.to_numeric(util["N"], errors="coerce")
        steps = [f"n{int(n):02d}" for n in sorted(util["N"].dropna().unique())
                 if (SCALE / f"n{int(n):02d}" / "utility" /
                     "learning_curve.csv").exists()]
        if not steps:
            st.info("No utility learning curves for these steps yet.")
            return
        st.subheader("Utility across cohorts")
        fig_show(pu.build_by_cohort(steps, util))
        fig_auc, fig_delta = pu.build_heatmaps(steps, util)
        fig_show(fig_auc)
        fig_show(fig_delta)
        st.caption(
            "Each point is the AUC at the largest synthetic size (10k), with "
            "the size-matched real control beside it; held-out donor counts are "
            "under the ticks. n16 is below chance for both modalities and is "
            "shown as such. See the Synthetic utility page for the full curve."
        )
        return

    metrics = sorted(sub["metric"].unique())
    metric = st.selectbox("Metric", metrics)
    d = sub[sub["metric"] == metric].copy()
    if d.empty:
        st.info("No values for this metric yet.")
        return

    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    d["synthetic_size"] = pd.to_numeric(d.get("synthetic_size"), errors="coerce")

    if analysis in ("utility", "utility_iter") and d["synthetic_size"].notna().any():
        # Keep the iteration axis: averaging over synthetic size would hide the curve.
        agg = (d.dropna(subset=["synthetic_size"])
                .groupby(["N", "modality", "label", "synthetic_size"])["value"]
                .mean().reset_index())
        fig, ax = plt.subplots(figsize=(9, 4.2))
        sns.lineplot(data=agg, x="synthetic_size", y="value", hue="modality",
                     style="label", marker="o", ax=ax)
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set_xlabel("synthetic training cells")
        ax.set_ylabel(metric)
        ax.set_title(f"{analysis}: {metric} (by synthetic size)")
        fig_show(fig)
        st.dataframe(agg.round(4), width="stretch", hide_index=True)
    else:
        agg = d.groupby(["N", "modality"])["value"].mean().reset_index()
        fig, ax = plt.subplots(figsize=(8, 3.8))
        sns.lineplot(data=agg, x="N", y="value", hue="modality", marker="o", ax=ax)
        ax.set_xlabel("donors (N)")
        ax.set_ylabel(metric)
        ax.set_title(f"{analysis}: {metric}")
        fig_show(fig)

        st.subheader("All values")
        st.dataframe(d.round(4), width="stretch", hide_index=True)

    st.caption(
        "A metric that improves with N is gaining donor-level information. A "
        "metric that improves only with synthetic size is resampling the same "
        "donor. `donor_auc` is undefined at n02, where condition == donor."
    )


def page_impact() -> None:
    st.title("Cell-type impact of disease")
    st.caption(f"Cohort: {COHORT_LABEL}.")
    steps = sorted([p.name for p in SCALE.glob("n*") if p.is_dir()])
    if not steps:
        missing_box("Scaling steps", "plan_scaling")
        return

    impact = read_csv(str(d("impact", "celltype_impact.csv")))
    if impact is None:
        missing_box(f"Cell-type impact for {COHORT_LABEL}", "impact")
        return

    for col in ("de_burden", "mean_abs_lfc_sig", "donor_auc",
                 "condition_variance_share", "pd_fraction",
                 "reversal_mean_abs_lfc"):
        if col in impact.columns:
            impact[col] = pd.to_numeric(impact[col], errors="coerce")

    st.subheader("Ranked impact")
    fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
    sns.barplot(data=impact, x="celltype", y="de_burden", hue="modality", ax=axes[0])
    axes[0].set_ylabel("significant / testable")
    axes[0].set_xlabel("")
    axes[0].tick_params(axis="x", rotation=25)
    sns.barplot(data=impact, x="celltype", y="donor_auc", hue="modality", ax=axes[1])
    axes[1].axhline(0.5, color="black", ls="--", lw=1)
    axes[1].set_ylim(0, 1.02)
    axes[1].set_ylabel("donor-held-out AUC")
    axes[1].set_xlabel("")
    axes[1].tick_params(axis="x", rotation=25)
    fig.tight_layout()
    fig_show(fig)

    st.dataframe(impact.round(3), width="stretch", hide_index=True)

    st.subheader("Donor-held-out separability vs N")
    frames = []
    for s in steps:
        p = SCALE / s / "impact" / "celltype_impact.csv"
        if not p.exists():
            continue
        t = pd.read_csv(p)
        t["N"] = int(s.replace("n", ""))
        frames.append(t)
    if frames:
        all_imp = pd.concat(frames, ignore_index=True)
        all_imp["donor_auc"] = pd.to_numeric(all_imp["donor_auc"], errors="coerce")
        fig, ax = plt.subplots(figsize=(8, 3.8))
        sns.lineplot(data=all_imp, x="N", y="donor_auc", hue="celltype",
                     style="modality", marker="o", ax=ax)
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("donor-held-out AUC")
        fig_show(fig)

    st.error(
        "At n02 there is one donor per condition, so `donor_auc` is undefined "
        "and every cell-level contrast is also a donor contrast. The metric only "
        "becomes a condition test from n08 onward."
    )


def page_utility() -> None:
    st.title("Synthetic-to-real utility")
    st.caption(
        f"Cohort: **{COHORT_LABEL}**. A classifier is trained on synthetic cells "
        "and tested on donors that were held out before simulation, so a high "
        "score means the synthetic data carries donor-level condition signal, "
        "not that it memorised the training donors."
    )

    curve = read_csv(str(d("utility", "learning_curve.csv")))
    raw = read_csv(str(d("utility", "utility.csv")))
    folds = read_csv(str(d("utility", "sim", "folds.csv")))
    summary = read_json(str(d("utility", "summary.json")))
    if curve is None:
        missing_box(f"Utility results for {COHORT_LABEL}", "utility")
        return

    for col in ("auc_mean", "auc_sd", "bal_mean", "bal_sd"):
        if col in curve.columns:
            curve[col] = pd.to_numeric(curve[col], errors="coerce")
    if "synthetic_size" in curve.columns:
        curve["synthetic_size"] = pd.to_numeric(curve["synthetic_size"],
                                                errors="coerce")

    st.subheader("At the largest synthetic size")
    for col, modality in zip(st.columns(2), ("rna", "atac")):
        ceiling = pu.single(curve, modality, "real_trained", "auc")
        floor = pu.single(curve, modality, "naive_poisson", "auc")
        x, mean, _ = pu.method_curve(curve, modality, "synthetic", "auc")
        _, matched, _ = pu.method_curve(curve, modality, "real_matched", "auc")
        syn = float(mean[-1]) if len(mean) else float("nan")
        mat = float(matched[-1]) if len(matched) else float("nan")
        top = f"{pu.SIZE_LABELS.get(int(x[-1]), int(x[-1]))}" if len(x) else "—"
        col.markdown(f"**{modality.upper()}** · synthetic at {top}")
        col.metric("synthetic", num(syn),
                   f"{num(syn - pu.CHANCE, signed=True)} vs chance")
        col.metric("matched real (same n)", num(mat),
                   f"{num(syn - mat, signed=True)} for synthetic")
        col.metric("real donors (ceiling)", num(ceiling),
                   f"{num(syn - ceiling, signed=True)} gap")
        col.caption(f"marginal-only floor: {num(floor)} AUC")

    st.subheader("Learning curve")
    fig = pu.build_learning_curve(
        DATA.name if DATA.name.startswith("n") else "cohort",
        curve, raw, folds,
        title=f"Synthetic-to-real utility \u2014 {COHORT_LABEL}")
    fig_show(fig)
    st.caption(
        "x is the number of synthetic training cells, on a log2 axis. "
        "`synthetic` and `matched real` are both trained on exactly that many "
        "cells, so their gap is the value of the synthetic data at equal "
        "sample size. The red dashed and orange dotted lines are reference "
        "levels (`real donors` uses every real training cell, `marginal only` "
        "is the naive baseline); neither is a point on the size axis. Grey is "
        "the same synthetic cells with shuffled labels."
    )

    util = pu.comparison_frame() if (SCALE / "comparison.csv").exists() else None
    if util is not None and not util.empty:
        util["N"] = pd.to_numeric(util["N"], errors="coerce")
        steps = [f"n{int(n):02d}" for n in
                 sorted(util["N"].dropna().unique())
                 if (SCALE / f"n{int(n):02d}" / "utility" /
                     "learning_curve.csv").exists()]
        if len(steps) > 1:
            st.subheader("Across cohorts")
            st.caption(
                "Held-out donor counts are shown under each tick. n16 is below "
                "chance for both modalities: with 5 held-out donors the "
                "condition contrast does not transfer, so the dip is real and "
                "is left visible rather than smoothed away."
            )
            fig_show(pu.build_by_cohort(steps, util))
            fig_auc, fig_delta = pu.build_heatmaps(steps, util)
            fig_show(fig_auc)
            st.caption(
                "Colour is centred on chance (0.5): red is above chance, blue "
                "below. The first column is the real-donor ceiling for that "
                "cohort."
            )
            fig_show(fig_delta)
            st.caption(
                "\u0394 AUC = synthetic \u2212 size-matched real, at identical "
                "training size. Red means the synthetic cells beat the same "
                "number of real cells; blue means they add nothing over real "
                "cells of the same count."
            )
    else:
        st.info(
            "`results/scale/comparison.csv` has no utility rows yet; run "
            "`scaling_compare` for the cross-cohort panels."
        )

    st.subheader("All rows")
    st.dataframe(curve.round(3), width="stretch", hide_index=True)

    if summary and summary.get("status") == "not_estimable":
        st.warning(
            "Not estimable: at least one condition has fewer than two donors, so "
            "a donor-disjoint train/test split cannot be formed."
        )

    st.caption(
        "Synthetic cells are drawn from the scDesign3 fit on training donors "
        "only, then subsampled without replacement. `real_trained` is the upper "
        "bound, `real_matched` is the same-n real control, `naive_poisson` is a "
        "marginal-only lower bound, and `permuted_labels` / `majority` are "
        "chance controls. A utility that rises with N is gaining donor-level "
        "information; one that rises only with synthetic size is resampling the "
        "same donors."
    )


def page_subtype() -> None:
    st.title("Subtype transfer")
    st.caption(
        f"Cohort: **{COHORT_LABEL}**. k=2 subclusters are defined on training "
        "donors only, separately in RNA and ATAC. Synthetic cells come from "
        "scDesign3 with `subtype + condition` (alt), a subtype-collapsed refit "
        "(null), or the alt fit with subtype coefficients zeroed (null_zero)."
    )

    transfer = read_csv(str(d("subtype", "transfer", "transfer.csv")))
    jac = read_csv(str(d("subtype", "transfer", "jaccard.csv")))
    comp = read_csv(str(d("subtype", "comparison.csv")))
    if transfer is None:
        missing_box("Subtype transfer", "subtype_transfer")
        return

    transfer["value"] = pd.to_numeric(transfer["value"], errors="coerce")
    celltypes = sorted(transfer["celltype"].dropna().unique())
    modalities = sorted(transfer["modality"].dropna().unique())
    if not celltypes or not modalities:
        st.info("No subtype transfer rows yet.")
        return
    c1, c2 = st.columns(2)
    ct = c1.selectbox("Cell type", celltypes)
    mod = c2.selectbox("Modality", modalities)
    sub = transfer[(transfer["celltype"] == ct) &
                    (transfer["modality"] == mod) &
                    (transfer["metric"] == "auc") &
                    (transfer["stratum"] == "overall")]

    order = ["real", "syn_alt", "syn_null", "syn_null_zero"]
    st.subheader("Transfer matrix (AUC)")
    lr = sub[sub["model"] == "logistic_regression"]
    if not lr.empty:
        piv = lr.pivot_table(index="train_domain", columns="test_domain",
                             values="value", aggfunc="mean")
        rows = [x for x in order if x in piv.index]
        cols = [x for x in order if x in piv.columns]
        piv = piv.reindex(index=rows, columns=cols)
        fig, ax = plt.subplots(figsize=(5.6, 4.2))
        sns.heatmap(piv, annot=True, fmt=".2f", cmap="Blues", vmin=0.4,
                     vmax=1.0, ax=ax)
        ax.set_xlabel("test domain")
        ax.set_ylabel("train domain")
        fig_show(fig)
        st.dataframe(piv.round(3), width="stretch")

    st.subheader("Headline pairs and controls")
    pairs = [
        ("real", "real", "logistic_regression", "real → real (ceiling)"),
        ("syn_alt", "real", "logistic_regression", "syn alt → real"),
        ("real", "syn_alt", "logistic_regression", "real → syn alt"),
        ("syn_null", "real", "logistic_regression", "syn null → real"),
        ("syn_null_zero", "real", "logistic_regression",
         "syn null_zero → real"),
        ("real", "real", "random_forest", "real → real (RF)"),
        ("syn_alt", "real", "random_forest", "syn alt → real (RF)"),
        ("real", "real", "donor_only", "donor-only ceiling"),
        ("real", "real", "condition_only", "condition-only control"),
        ("real", "real", "permuted_labels", "permuted labels"),
        ("real", "real", "majority", "majority"),
    ]
    rows = []
    for train, test, model, label in pairs:
        g = sub[(sub["train_domain"] == train) &
                (sub["test_domain"] == test) &
                (sub["model"] == model)]
        if g.empty:
            continue
        rows.append({
            "pair": label,
            "AUC": g["value"].mean(),
            "SD": g["value"].std(),
            "n seeds": g["value"].size,
        })
    if rows:
        st.dataframe(pd.DataFrame(rows).round(3), width="stretch",
                     hide_index=True)

    cond = transfer[(transfer["celltype"] == ct) &
                    (transfer["modality"] == mod) &
                    (transfer["metric"] == "auc") &
                    (transfer["model"] == "logistic_regression") &
                    (transfer["stratum"].astype(str).str.startswith("condition:"))]
    if not cond.empty:
        st.subheader("Condition-stratified AUC")
        fig, ax = plt.subplots(figsize=(7.5, 3.6))
        sns.barplot(data=cond, x="stratum", y="value", hue="train_domain",
                    ax=ax)
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("")
        ax.tick_params(axis="x", rotation=20)
        fig_show(fig)

    donor = transfer[(transfer["celltype"] == ct) &
                      (transfer["modality"] == mod) &
                      (transfer["metric"] == "auc") &
                      (transfer["model"] == "logistic_regression") &
                      (transfer["stratum"].astype(str).str.startswith("donor:")) &
                      (transfer["test_domain"] == "real")]
    if not donor.empty:
        st.subheader("Donor-stratified AUC")
        agg = (donor.groupby(["train_domain", "stratum"])["value"]
               .mean().reset_index())
        fig, ax = plt.subplots(figsize=(8.5, 3.6))
        sns.lineplot(data=agg, x="stratum", y="value", hue="train_domain",
                     marker="o", ax=ax)
        ax.axhline(0.5, color="black", ls="--", lw=1)
        ax.set_ylim(0, 1.02)
        ax.set_ylabel("AUC")
        ax.set_xlabel("held-out donor")
        ax.tick_params(axis="x", rotation=70)
        fig_show(fig)

    if jac is not None and not jac.empty:
        st.subheader("Variable-feature overlap")
        st.dataframe(jac.round(3), width="stretch", hide_index=True)

    if comp is not None and not comp.empty:
        st.subheader("QC and concordance")
        sub2 = comp[(comp["celltype"] == ct) &
                     (comp["analysis"].isin(["qc", "concordance"]))]
        st.dataframe(sub2.round(4), width="stretch", hide_index=True)

    st.caption(
        "`real → real` is the donor-disjoint ceiling. `syn alt → real` asks whether "
        "the alternative simulation preserves the subtype axis; `real → syn alt` is the "
        "reverse direction. The two null simulations should sit near chance when trained "
        "on synthetic and tested on real. `donor_only` is a donor-mean ceiling, and "
        "`condition_only` checks whether the split is just condition. Microglia RNA and "
        "ATAC labels can disagree, so cross-modality concordance is reported per cell type."
    )


def page_methods() -> None:
    st.title("Methods and caveats")
    cfg = config()

    st.subheader("Configuration")
    if cfg:
        st.json(cfg, expanded=False)

    st.subheader("Decisions")
    st.markdown(
        """
- **Clustering** uses Harmony over donor because donor and condition are
  perfectly confounded. The embedding is used only for clustering and
  annotation; every statistic downstream uses raw counts.
- **Feature selection** is fit on the fit split only: 2,000 RNA genes by
  SCTransform-v2 residual variance and 2,000 ATAC peaks by TF-IDF variance,
  each with a 5% detection floor (ATAC requires it in both groups).
- **scDesign3** fits RNA and ATAC separately — there is no joint cross-modality
  copula — and pairs modalities only through shared covariates. The interaction
  model is `celltype + condition + celltype:condition`.
- **Null** simulations zero every `condition` coefficient in the fitted marginal
  models, so cell-type structure survives but the condition effect does not.
- **Fidelity** uses a real-vs-synthetic discriminator and cross-domain label
  transfer, with identity, positive (Aged vs PD) and naive-Poisson controls.
- **Reversal** ranks features two ways: the deterministic scDesign3
  counterfactual (`extract_para`) and a classifier contribution. Only strata
  with ≥100 cells per condition are reported.
- **Cell-type impact** ranks cell types by condition DE burden, a
  donor-held-out condition classifier, the share of variance explained by
  condition versus donor, and the counterfactual magnitude.
- **Synthetic utility** trains a condition classifier on synthetic cells from
  training donors and tests it on held-out real donors. Real-trained,
  size-matched real, naive-Poisson, permuted-label and majority baselines
  bracket the result; the size-matched control is trained on the same number of
  real cells as each synthetic size, so the synthetic-minus-real gap is not
  confounded with less training data.
- **Scaling** repeats every stage on nested donor sets (n02, n08, n16, n24) with
  a fixed n24 feature set and 5 simulation seeds.
- **Subtype transfer** subclusters Oligodendrocyte and Microglia into k=2 on
  training donors only, filters doublets and the opposite-lineage marker, regresses
  depth out of the embedding, and winsorises Harmony dimensions before k-means.
  scDesign3 is fit with `subtype + condition` (alt), a collapsed `condition`
  refit (null) and a zero-subtype-coefficient sensitivity null. Transfer is
  evaluated in both directions with donor-disjoint real test cells.
"""
    )

    st.subheader("Limitations")
    st.markdown(
        """
- **The pilot has two donors, one per condition.** Every pilot Aged-vs-PD
  number is also a donor contrast; there is no biological replication. The
  scaling steps add donors and report the donor-held-out AUC precisely so this
  changes.
- **No dopaminergic neurons.** `TH` has zero UMIs in C16007; across all 34
  samples the maximum is 35 `TH`+ nuclei. The available classes are glial.
- **Oligodendrocyte-dominated.** ~88% of pilot nuclei are oligodendrocyte
  lineage, so power is concentrated there and other classes are reported
  descriptively.
- **Marginal condition effect.** The scDesign3 mean model has no cell-level
  covariates, so its counterfactual is group-level; per-cell variation comes
  from the classifier branch.
- **Donor-disjoint utility starts at n08.** At n02 there is one donor per
  condition, so a donor-disjoint train/test split cannot be formed and the
  utility table is marked `not_estimable`.
- **Feature set fixed at n24.** Smaller steps reuse the n24 genes/peaks so the
  scaling curve isolates sample size; the per-step Jaccard overlap is reported.
- **Reversal permutation is pilot-only.** Refitting the marginal model with
  permuted labels once per cell type is affordable at n02 but not at n24, so
  scaling steps report the deterministic counterfactual without that control.
- **Subtypes are unsupervised axes, not known biology.** The k=2 splits are
  defined on training donors and are weak (silhouette 0.14–0.26). Microglia RNA
  and ATAC labels disagree (ARI ≈0), so cross-modality transfer there is a
  negative result; Oligodendrocyte ARI is ≈0.19. Do not read the split as a
  validated developmental subtype.
"""
    )

    st.subheader("Reproduce")
    st.code(
        "snakemake -n\n"
        "snakemake -j 8 --use-conda\n"
        "streamlit run dashboard/app.py",
        language="bash",
    )


PAGES = {
    "Overview": page_overview,
    "Report": page_report,
    "QC & composition": page_qc,
    "Annotation": page_annotation,
    "Fidelity": page_fidelity,
    "Hypotheses": page_hypothesis,
    "Reversal": page_reversal,
    "Cell-type impact": page_impact,
    "Synthetic utility": page_utility,
    "Sample-size scaling": page_scaling,
    "Subtype transfer": page_subtype,
    "Methods & caveats": page_methods,
}

st.sidebar.title("GSE193688")
st.sidebar.caption("snRNA + snATAC · substantia nigra")

cohorts = {k: v for k, v in STEP_DIRS.items()
           if (v / "cluster" / "obs_annotated.csv").exists()}
if not cohorts:
    cohorts = {"pilot — 2 donors (manual annotation)": RESULTS}
COHORT_LABEL = st.sidebar.selectbox("Cohort", list(cohorts))
DATA = cohorts[COHORT_LABEL]

choice = st.sidebar.radio("Page", list(PAGES), label_visibility="collapsed")
st.sidebar.divider()
st.sidebar.caption(
    f"Showing **{COHORT_LABEL}**.\n\n"
    "Scaling steps: n02 → n08 → n16 → n24.\n\n"
    "Condition is confounded with donor at n02."
)
PAGES[choice]()
