"""Snakemake entry point. Run from the project root:

    snakemake -n                      # dry run
    snakemake -j 8 --use-conda        # real run

The pilot (2 donors) and every scaling step (n02, n08, n16, n24) use the same
rules. A step is selected with a generated config file:

    snakemake --configfile config/scaling/n08.yaml -j 8 --use-conda

`src/plan_scaling.py` writes `config/subsets/nXX.csv` and
`config/scaling/nXX.yaml`, so no rule needs to know the sample size.

The annotation gate is intentionally a hard stop in the pilot: `snakemake` runs
clustering, writes `config/celltype_map.csv` as a template, and fails on `annotate`
until a human fills in the `celltype` column. Scaling steps set
`cluster.auto_annotate: true` and use `src/auto_annotate.py` instead.
"""

configfile: "config/config.yaml"

OUT = config.get("results_dir", "results")
PROCESSED = config.get("processed_dir", "data/processed")
SHEET = config["samples"]["sheet"]
INGEST = f"{OUT}/ingest"
CLUSTER = f"{OUT}/cluster"
SPLIT = f"{OUT}/split"
FEATURES = f"{OUT}/features"
SIM = f"{OUT}/simulation"
FIDELITY = f"{OUT}/fidelity"
HYPOTHESIS = f"{OUT}/hypothesis"
REVERSAL = f"{OUT}/reversal"
IMPACT = f"{OUT}/impact"
UTILITY = f"{OUT}/utility"
DASHBOARD = f"{OUT}/dashboard"

INGEST_FILES = [
    f"{INGEST}/{name}"
    for name in (
        "obs.csv",
        "rna_counts.mtx.gz",
        "atac_counts.mtx.gz",
        "rna_features.tsv",
        "atac_features.tsv",
        "consensus_peaks.bed",
        "qc_per_sample.csv",
        "sessionInfo.txt",
    )
]

AUTO = bool(config["cluster"].get("auto_annotate", False))
MAP = config["cluster"].get("map", "config/celltype_map.csv")
if AUTO:
    # never clobber the human-curated pilot map
    MAP = f"{CLUSTER}/celltype_map.csv"
FIXED_FEATURES = config["features"].get("fixed_from", "")
SCALING = config.get("scaling", {})
SUBSET = bool(SCALING.get("subset", False))
# The h5ad exports are an archive/QC artifact; no downstream rule reads them.
# Scaling steps turn this off to avoid materialising a 189k x 37k sparse matrix.
ASSEMBLE = bool(SCALING.get("assemble", True))
FULL_INGEST = SCALING.get("full_ingest", "")
INGEST_ONLY = bool(SCALING.get("ingest_only", False))
REGENERATE_SHEET = bool(config["samples"].get("regenerate", True))

ALL_TARGETS = [
    f"{FEATURES}/rna_features.tsv",
    f"{FEATURES}/atac_features.tsv",
    f"{SIM}/manifest.csv",
    f"{FIDELITY}/summary.json",
    f"{HYPOTHESIS}/summary.json",
    f"{DASHBOARD}/manifest.json",
]
if ASSEMBLE:
    ALL_TARGETS[:0] = [f"{PROCESSED}/rna.h5ad", f"{PROCESSED}/atac.h5ad"]
REVERSAL_ENABLED = bool(config.get("reversal", {}).get("enabled", False))
if REVERSAL_ENABLED:
    ALL_TARGETS.append(f"{REVERSAL}/summary.json")
if config.get("reversal", {}).get("poc", False):
    ALL_TARGETS.append(f"{REVERSAL}/poc/summary.json")
if config.get("impact", {}).get("enabled", False):
    ALL_TARGETS.append(f"{IMPACT}/summary.json")
if config.get("utility", {}).get("enabled", False):
    ALL_TARGETS.append(f"{UTILITY}/summary.json")

# ------------------------------------------------- subtype transfer ----
SUBTYPE = config.get("subtype", {})
SUBTYPE_ENABLED = bool(SUBTYPE.get("enabled", False))
SUBTYPE_ROOT = f"{OUT}/subtype"
SUBTYPE_TRANSFER = f"{SUBTYPE_ROOT}/transfer"
SUBTYPE_CELLTYPES = list(SUBTYPE.get("celltypes", []))
SUBTYPE_MODALITIES = list(SUBTYPE.get("modalities", []))
SUBTYPE_MODES = list(SUBTYPE.get("modes", ["alt", "null", "null_zero"]))
SUBTYPE_SEEDS = [int(s) for s in SUBTYPE.get("seeds", [])]
SUBTYPE_CLUSTER_FILES: list[str] = []
SUBTYPE_SIM_FILES: list[str] = []
for _ct in SUBTYPE_CELLTYPES:
    for _mod in SUBTYPE_MODALITIES:
        _base = f"{SUBTYPE_ROOT}/{_ct}/{_mod}"
        SUBTYPE_CLUSTER_FILES += [
            f"{_base}/labels.csv", f"{_base}/embedding.csv",
            f"{_base}/markers.csv", f"{_base}/features.csv",
            f"{_base}/filtered_cells.csv", f"{_base}/qc.json",
        ]
        SUBTYPE_SIM_FILES.append(f"{_base}/simulation/manifest.csv")
        for _mode in SUBTYPE_MODES:
            SUBTYPE_SIM_FILES.append(f"{_base}/simulation/{_mode}/obs.csv")
            SUBTYPE_SIM_FILES += [
                f"{_base}/simulation/{_mode}/seed{_s}_counts.mtx.gz"
                for _s in SUBTYPE_SEEDS
            ]
SUBTYPE_FILES = [
    f"{SUBTYPE_ROOT}/donors.csv",
    f"{SUBTYPE_ROOT}/dropped.csv",
    f"{SUBTYPE_ROOT}/split_summary.csv",
] + SUBTYPE_CLUSTER_FILES + SUBTYPE_SIM_FILES + [
    f"{SUBTYPE_TRANSFER}/transfer.csv",
    f"{SUBTYPE_TRANSFER}/transfer_summary.csv",
    f"{SUBTYPE_TRANSFER}/jaccard.csv",
    f"{SUBTYPE_TRANSFER}/summary.json",
    f"{SUBTYPE_ROOT}/comparison.csv",
    f"{SUBTYPE_ROOT}/summary.json",
]
if SUBTYPE_ENABLED:
    ALL_TARGETS.append(f"{SUBTYPE_ROOT}/comparison.csv")
    ALL_TARGETS.append(f"{SUBTYPE_ROOT}/summary.json")
if INGEST_ONLY:
    ALL_TARGETS = list(INGEST_FILES)


rule all:
    input:
        ALL_TARGETS


if REGENERATE_SHEET:
    rule sample_sheet:
        """Regenerate the sample sheet from the files in the data directory."""
        output:
            SHEET,
        conda:
            "envs/py_env.yaml"
        script:
            "src/make_sample_sheet.py"


if SUBSET:
    rule subset_ingest:
        """Subset the one-time full ingest to this step's donors."""
        input:
            f"{FULL_INGEST}/obs.csv",
            f"{FULL_INGEST}/rna_counts.mtx.gz",
            f"{FULL_INGEST}/atac_counts.mtx.gz",
            f"{FULL_INGEST}/rna_features.tsv",
            f"{FULL_INGEST}/atac_features.tsv",
            f"{FULL_INGEST}/consensus_peaks.bed",
            f"{FULL_INGEST}/qc_per_sample.csv",
            f"{FULL_INGEST}/sessionInfo.txt",
        output:
            INGEST_FILES,
        conda:
            "envs/py_env.yaml"
        params:
            full=FULL_INGEST,
            samples=SHEET,
            out=INGEST,
        shell:
            """
            python src/subset_ingest.py \
                --full {params.full} \
                --samples {params.samples} \
                --out {params.out}
            """
else:
    rule ingest:
        """QC + doublets + consensus peaks + exact ATAC recount. Needs R."""
        input:
            SHEET,
        output:
            INGEST_FILES,
        conda:
            "envs/r_env.yaml"
        params:
            samples=SHEET,
            min_features_rna=config["qc"]["min_features_rna"],
            max_features_rna=config["qc"]["max_features_rna"],
            min_counts_rna=config["qc"]["min_counts_rna"],
            max_percent_mt=config["qc"]["max_percent_mt"],
            min_counts_atac=config["qc"]["min_counts_atac"],
            peak_min_width=config["peaks"]["min_width"],
            peak_max_width=config["peaks"]["max_width"],
            consensus_from=config["samples"]["consensus_from"],
            doublets=str(config["qc"]["doublets"]["enabled"]).lower(),
        shell:
            """
            Rscript R/00_ingest.R \
                --samples {params.samples} \
                --out {INGEST} \
                --min-features-rna {params.min_features_rna} \
                --max-features-rna {params.max_features_rna} \
                --min-counts-rna {params.min_counts_rna} \
                --max-percent-mt {params.max_percent_mt} \
                --min-counts-atac {params.min_counts_atac} \
                --peak-min-width {params.peak_min_width} \
                --peak-max-width {params.peak_max_width} \
                --consensus-from {params.consensus_from} \
                --doublets {params.doublets}
            """


rule subset_features:
    """Stream only the selected features out of the ingest matrices.

    Simulation only needs the ~2,000 selected genes/peaks, but reading the full
    ATAC matrix (280M nonzeros at n24) peaks at several GB and repeatedly hit
    the OOM killer on this VM. Pre-extracting the selected rows keeps every
    downstream R job small.
    """
    input:
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
    output:
        rna=f"{INGEST}/rna_counts_sel.mtx.gz",
        atac=f"{INGEST}/atac_counts_sel.mtx.gz",
        rna_feats=f"{INGEST}/rna_features_sel.tsv",
        atac_feats=f"{INGEST}/atac_features_sel.tsv",
    conda:
        "envs/py_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
    shell:
        """
        python src/subset_features.py \
            --ingest {params.ingest} \
            --features {params.features}
        """


if ASSEMBLE:
    rule assemble:
        """MTX -> AnnData. Asserts RNA/ATAC cell order and integer counts."""
        input:
            INGEST_FILES,
        output:
            rna=f"{PROCESSED}/rna.h5ad",
            atac=f"{PROCESSED}/atac.h5ad",
        conda:
            "envs/py_env.yaml"
        shell:
            "python src/assemble.py --ingest {INGEST} --out {PROCESSED}"


rule cluster:
    """RNA PCA + ATAC LSI + joint WNN clustering, plus marker table."""
    input:
        INGEST_FILES,
    output:
        clusters=f"{CLUSTER}/clusters.csv",
        markers=f"{CLUSTER}/cluster_markers.csv",
        sizes=f"{CLUSTER}/cluster_sizes.csv",
        session=f"{CLUSTER}/sessionInfo.txt",
    conda:
        "envs/r_env.yaml"
    params:
        ingest=INGEST,
        map=MAP,
        resolution=config["cluster"]["resolution"],
        dims_rna=config["cluster"]["dims_rna"],
        dims_atac=config["cluster"]["dims_atac"],
        n_hvg=config["cluster"]["n_hvg"],
        integrate=config["cluster"]["integrate"],
        seed=config["split"]["seed"],
    shell:
        """
        Rscript R/01_cluster.R \
            --ingest {params.ingest} \
            --out {CLUSTER} \
            --map {params.map} \
            --resolution {params.resolution} \
            --dims-rna {params.dims_rna} \
            --dims-atac {params.dims_atac} \
            --n-hvg {params.n_hvg} \
            --integrate {params.integrate} \
            --seed {params.seed}
        """


if AUTO:
    rule annotate:
        """Assign cell types from canonical marker scores (scaling runs)."""
        input:
            INGEST_FILES,
            clusters=f"{CLUSTER}/clusters.csv",
            markers=f"{CLUSTER}/cluster_markers.csv",
        output:
            obs=f"{CLUSTER}/obs_annotated.csv",
            cellmap=f"{CLUSTER}/celltype_map.csv",
            scores=f"{CLUSTER}/cluster_annotation_scores.csv",
        conda:
            "envs/py_env.yaml"
        params:
            ingest=INGEST,
            clusters=f"{CLUSTER}/clusters.csv",
            markers=f"{CLUSTER}/cluster_markers.csv",
            out=CLUSTER,
            min_score=config["cluster"].get("auto_min_score", 0.5),
        shell:
            """
            python src/auto_annotate.py \
                --ingest {params.ingest} \
                --clusters {params.clusters} \
                --markers {params.markers} \
                --out {params.out} \
                --min-score {params.min_score}
            """
else:
    rule annotate:
        """Apply the manual cell-type map. Fails until every cluster is labelled."""
        input:
            f"{CLUSTER}/clusters.csv",
        output:
            f"{CLUSTER}/obs_annotated.csv",
        conda:
            "envs/py_env.yaml"
        params:
            clusters=f"{CLUSTER}/clusters.csv",
            map=MAP,
            out=f"{CLUSTER}/obs_annotated.csv",
        shell:
            "python src/annotate.py --clusters {params.clusters} --map {params.map} --out {params.out}"


rule split:
    """Stratified fit/holdout split by celltype x donor."""
    input:
        f"{CLUSTER}/obs_annotated.csv",
    output:
        f"{SPLIT}/split.csv",
        f"{SPLIT}/fit_cells.txt",
        f"{SPLIT}/holdout_cells.txt",
    conda:
        "envs/py_env.yaml"
    params:
        obs=f"{CLUSTER}/obs_annotated.csv",
        out=SPLIT,
        fraction=config["split"]["holdout_fraction"],
        seed=config["split"]["seed"],
    shell:
        """
        python src/split.py \
            --obs {params.obs} \
            --out {params.out} \
            --fraction {params.fraction} \
            --seed {params.seed}
        """


if FIXED_FEATURES:
    rule features:
        """Reuse the n24 feature set; also record this step's own selection."""
        input:
            INGEST_FILES,
            f"{SPLIT}/split.csv",
            f"{CLUSTER}/obs_annotated.csv",
            f"{FIXED_FEATURES}/rna_features.tsv",
            f"{FIXED_FEATURES}/atac_features.tsv",
            f"{FIXED_FEATURES}/feature_params.json",
        output:
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{FEATURES}/feature_params.json",
            f"{FEATURES}/selected/rna_features.tsv",
            f"{FEATURES}/selected/atac_features.tsv",
            f"{FEATURES}/selected/feature_params.json",
        conda:
            "envs/r_env.yaml"
        params:
            ingest=INGEST,
            split=f"{SPLIT}/split.csv",
            obs=f"{CLUSTER}/obs_annotated.csv",
            out=FEATURES,
            fixed=FIXED_FEATURES,
            n_rna=config["features"]["rna"]["n"],
            n_atac=config["features"]["atac"]["n"],
            min_fraction_rna=config["features"]["rna"]["min_fraction"],
            min_fraction_atac=config["features"]["atac"]["min_fraction"],
            seed=config["split"]["seed"],
        shell:
            """
            mkdir -p {params.out}/selected
            Rscript R/02_features.R \
                --ingest {params.ingest} \
                --split {params.split} \
                --obs {params.obs} \
                --out {params.out}/selected \
                --n-rna {params.n_rna} \
                --n-atac {params.n_atac} \
                --min-fraction-rna {params.min_fraction_rna} \
                --min-fraction-atac {params.min_fraction_atac} \
                --seed {params.seed}
            cp {params.fixed}/rna_features.tsv {params.out}/rna_features.tsv
            cp {params.fixed}/atac_features.tsv {params.out}/atac_features.tsv
            cp {params.fixed}/feature_params.json {params.out}/feature_params.json
            """
else:
    rule features:
        """SCT residuals (RNA) + TF-IDF variance (ATAC), fit cells only."""
        input:
            INGEST_FILES,
            f"{SPLIT}/split.csv",
            f"{CLUSTER}/obs_annotated.csv",
        output:
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{FEATURES}/feature_params.json",
        conda:
            "envs/r_env.yaml"
        params:
            ingest=INGEST,
            split=f"{SPLIT}/split.csv",
            obs=f"{CLUSTER}/obs_annotated.csv",
            out=FEATURES,
            n_rna=config["features"]["rna"]["n"],
            n_atac=config["features"]["atac"]["n"],
            min_fraction_rna=config["features"]["rna"]["min_fraction"],
            min_fraction_atac=config["features"]["atac"]["min_fraction"],
            seed=config["split"]["seed"],
        shell:
            """
            Rscript R/02_features.R \
                --ingest {params.ingest} \
                --split {params.split} \
                --obs {params.obs} \
                --out {params.out} \
                --n-rna {params.n_rna} \
                --n-atac {params.n_atac} \
                --min-fraction-rna {params.min_fraction_rna} \
                --min-fraction-atac {params.min_fraction_atac} \
                --seed {params.seed}
            """


# --------------------------------------------------------------- simulation --
SIM_MODES = ["alt", "null"]
SIM_SEEDS = [int(s) for s in config["simulation"]["seeds"]]
SIM_FILES = [f"{SIM}/manifest.csv"]
SIM_FILES += [
    f"{SIM}/{m}/{mode}/seed{s}_counts.mtx.gz"
    for m in ("rna", "atac") for mode in SIM_MODES for s in SIM_SEEDS
]
SIM_FILES += [f"{SIM}/{m}/{mode}/obs.csv" for m in ("rna", "atac") for mode in SIM_MODES]
SIM_FILES += [f"{SIM}/fit_{m}.rds" for m in ("rna", "atac")]


rule simulate:
    """Fit scDesign3 once per modality, then simulate alt (observed) and null."""
    input:
        INGEST_FILES,
        f"{INGEST}/rna_counts_sel.mtx.gz",
        f"{INGEST}/atac_counts_sel.mtx.gz",
        f"{INGEST}/rna_features_sel.tsv",
        f"{INGEST}/atac_features_sel.tsv",
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    output:
        SIM_FILES,
    conda:
        "envs/r_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
        split=f"{SPLIT}/split.csv",
        obs=f"{CLUSTER}/obs_annotated.csv",
        out=SIM,
        family_use=config["simulation"]["family_use"],
        mu_formula=config["simulation"]["mu_formula"],
        sigma_formula=config["simulation"]["sigma_formula"],
        corr_formula=config["simulation"]["corr_formula"],
        n_cores=config["simulation"]["n_cores"],
        seeds=",".join(str(s) for s in SIM_SEEDS),
        modes=",".join(SIM_MODES),
        match_composition=str(config["simulation"]["match_composition"]).lower(),
    shell:
        """
        Rscript R/03_simulate.R \
            --ingest {params.ingest} \
            --features {params.features} \
            --split {params.split} \
            --obs {params.obs} \
            --out {params.out} \
            --family-use {params.family_use} \
            --mu-formula "{params.mu_formula}" \
            --sigma-formula "{params.sigma_formula}" \
            --corr-formula "{params.corr_formula}" \
            --n-cores {params.n_cores} \
            --seeds {params.seeds} \
            --modes {params.modes} \
            --match-composition {params.match_composition}
        """


# --------------------------------------------------------------- evaluation --
rule fidelity:
    """Discriminator AUC + cross-domain label transfer on the alt simulation."""
    input:
        SIM_FILES,
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    output:
        disc=f"{FIDELITY}/discriminator.csv",
        transfer=f"{FIDELITY}/label_transfer.csv",
        per_class=f"{FIDELITY}/per_class.csv",
        summary=f"{FIDELITY}/summary.json",
    conda:
        "envs/py_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
        split=f"{SPLIT}/split.csv",
        obs=f"{CLUSTER}/obs_annotated.csv",
        sim=SIM,
        out=FIDELITY,
        models=",".join(config["evaluation"]["models"]),
        folds=config["evaluation"]["cv_folds"],
    shell:
        """
        python src/fidelity.py \
            --ingest {params.ingest} \
            --features {params.features} \
            --split {params.split} \
            --obs {params.obs} \
            --sim {params.sim} \
            --out {params.out} \
            --mode alt \
            --models {params.models} \
            --folds {params.folds}
        """


rule hypothesis:
    """Null calibration and alt power for the condition and celltype contrasts."""
    input:
        SIM_FILES,
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    output:
        tests=f"{HYPOTHESIS}/tests.csv",
        calibration=f"{HYPOTHESIS}/calibration.csv",
        truth=f"{HYPOTHESIS}/truth.csv",
        summary=f"{HYPOTHESIS}/summary.json",
    conda:
        "envs/py_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
        split=f"{SPLIT}/split.csv",
        obs=f"{CLUSTER}/obs_annotated.csv",
        sim=SIM,
        out=HYPOTHESIS,
        fdr_grid=",".join(str(v) for v in config["hypothesis"]["target_fdr"]),
        test=config["hypothesis"]["de_tests"][0],
    shell:
        """
        python src/hypothesis.py \
            --ingest {params.ingest} \
            --features {params.features} \
            --split {params.split} \
            --obs {params.obs} \
            --sim {params.sim} \
            --out {params.out} \
            --fdr-grid {params.fdr_grid} \
            --test {params.test}
        """


# ----------------------------------------------------------------- reversal --
rule reversal_model:
    """Counterfactual PD -> Aged from the fitted scDesign3 marginal models."""
    input:
        f"{SIM}/fit_rna.rds",
        f"{SIM}/fit_atac.rds",
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    output:
        delta=f"{REVERSAL}/model_delta.csv",
        summary=f"{REVERSAL}/model_summary.json",
    conda:
        "envs/r_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
        split=f"{SPLIT}/split.csv",
        obs=f"{CLUSTER}/obs_annotated.csv",
        sim=SIM,
        out=REVERSAL,
        target=config["reversal"]["source"],
        reference=config["reversal"]["target"],
        min_cells=config["reversal"]["min_cells_per_stratum"],
        mu_formula=config["reversal"]["model"],
        n_cores=config["simulation"]["n_cores"],
        seed=config["split"]["seed"],
        permute=str(config["reversal"].get(
            "permute", "permuted_condition" in config["reversal"]["controls"]
        )).lower(),
    shell:
        """
        Rscript R/04_reversal.R \
            --ingest {params.ingest} \
            --features {params.features} \
            --split {params.split} \
            --obs {params.obs} \
            --sim {params.sim} \
            --out {params.out} \
            --target {params.target} \
            --reference {params.reference} \
            --min-cells-per-stratum {params.min_cells} \
            --mu-formula "{params.mu_formula}" \
            --n-cores {params.n_cores} \
            --seed {params.seed} \
            --permute {params.permute}
        """


rule reversal_classifier:
    """Classifier contribution to PD -> Aged, plus agreement with the model."""
    input:
        f"{REVERSAL}/model_delta.csv",
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    output:
        contribution=f"{REVERSAL}/classifier_contribution.csv",
        combined=f"{REVERSAL}/combined.csv",
        top=f"{REVERSAL}/top_features.csv",
        summary=f"{REVERSAL}/summary.json",
    conda:
        "envs/py_env.yaml"
    params:
        ingest=INGEST,
        features=FEATURES,
        split=f"{SPLIT}/split.csv",
        obs=f"{CLUSTER}/obs_annotated.csv",
        reversal=REVERSAL,
        out=REVERSAL,
        model="logistic_regression",
        folds=config["evaluation"]["cv_folds"],
        target=config["reversal"]["source"],
        reference=config["reversal"]["target"],
        seed=config["split"]["seed"],
    shell:
        """
        python src/reversal.py \
            --ingest {params.ingest} \
            --features {params.features} \
            --split {params.split} \
            --obs {params.obs} \
            --reversal {params.reversal} \
            --out {params.out} \
            --model {params.model} \
            --folds {params.folds} \
            --target {params.target} \
            --reference {params.reference} \
            --seed {params.seed}
        """


# --------------------------------------------------------- cell-type impact --
if config.get("impact", {}).get("enabled", False):
    _impact_input = [
        INGEST_FILES,
        f"{FEATURES}/rna_features.tsv",
        f"{FEATURES}/atac_features.tsv",
        f"{SPLIT}/split.csv",
        f"{CLUSTER}/obs_annotated.csv",
    ]
    if config.get("reversal", {}).get("enabled", False):
        _impact_input.append(f"{REVERSAL}/model_delta.csv")

    rule impact:
        """Rank cell types by DE burden, donor-held-out separability and variance."""
        input:
            _impact_input,
        output:
            table=f"{IMPACT}/celltype_impact.csv",
            de=f"{IMPACT}/de_features.csv",
            summary=f"{IMPACT}/summary.json",
        conda:
            "envs/py_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            split=f"{SPLIT}/split.csv",
            obs=f"{CLUSTER}/obs_annotated.csv",
            reversal=REVERSAL,
            out=IMPACT,
            fdr=config["impact"]["fdr"],
            min_cells=config["impact"]["min_cells_per_stratum"],
            folds=config["impact"]["cv_folds"],
            seed=config["split"]["seed"],
            has_reversal=str(config.get("reversal", {}).get("enabled", False)).lower(),
        shell:
            """
            python src/celltype_impact.py \
                --ingest {params.ingest} \
                --features {params.features} \
                --split {params.split} \
                --obs {params.obs} \
                --reversal {params.reversal} \
                --out {params.out} \
                --fdr {params.fdr} \
                --min-cells {params.min_cells} \
                --folds {params.folds} \
                --seed {params.seed} \
                --has-reversal {params.has_reversal}
            """


# ------------------------------------------------------ synthetic utility ----
if config.get("utility", {}).get("enabled", False):
    rule utility_sim:
        """Fit scDesign3 on training donors and simulate donor-disjoint cohorts."""
        input:
            INGEST_FILES,
            f"{INGEST}/rna_counts_sel.mtx.gz",
            f"{INGEST}/atac_counts_sel.mtx.gz",
            f"{INGEST}/rna_features_sel.tsv",
            f"{INGEST}/atac_features_sel.tsv",
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{CLUSTER}/obs_annotated.csv",
        output:
            manifest=f"{UTILITY}/sim/manifest.csv",
            folds=f"{UTILITY}/sim/folds.csv",
        conda:
            "envs/r_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            obs=f"{CLUSTER}/obs_annotated.csv",
            out=f"{UTILITY}/sim",
            seed=config["split"]["seed"],
            sizes=",".join(str(v) for v in config["utility"]["synthetic_sizes"]),
            sim_seeds=",".join(str(v) for v in config["utility"]["seeds"]),
            holdout_fraction=config["utility"]["donor_holdout_fraction"],
            folds=config["utility"].get("donor_folds", 1),
            n_cores=config["simulation"]["n_cores"],
            family_use=config["simulation"]["family_use"],
            mu_formula=config["simulation"]["mu_formula"],
            sigma_formula=config["simulation"]["sigma_formula"],
            corr_formula=config["simulation"]["corr_formula"],
        shell:
            """
            Rscript R/05_utility_sim.R \
                --ingest {params.ingest} \
                --features {params.features} \
                --obs {params.obs} \
                --out {params.out} \
                --seed {params.seed} \
                --sizes {params.sizes} \
                --sim-seeds {params.sim_seeds} \
                --holdout-fraction {params.holdout_fraction} \
                --folds {params.folds} \
                --n-cores {params.n_cores} \
                --family-use {params.family_use} \
                --mu-formula "{params.mu_formula}" \
                --sigma-formula "{params.sigma_formula}" \
                --corr-formula "{params.corr_formula}"
            """

    rule utility:
        """Train on synthetic cells, test on held-out real donors."""
        input:
            f"{UTILITY}/sim/manifest.csv",
            f"{UTILITY}/sim/folds.csv",
            INGEST_FILES,
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{SPLIT}/split.csv",
            f"{CLUSTER}/obs_annotated.csv",
        output:
            table=f"{UTILITY}/utility.csv",
            curve=f"{UTILITY}/learning_curve.csv",
            summary=f"{UTILITY}/summary.json",
        conda:
            "envs/py_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            obs=f"{CLUSTER}/obs_annotated.csv",
            sim=f"{UTILITY}/sim",
            out=UTILITY,
            label=config["utility"]["label"],
            folds=config["utility"]["folds"],
            sizes=",".join(str(v) for v in config["utility"]["synthetic_sizes"]),
            seed=config["split"]["seed"],
        shell:
            """
            python src/synthetic_utility.py \
                --ingest {params.ingest} \
                --features {params.features} \
                --obs {params.obs} \
                --sim {params.sim} \
                --out {params.out} \
                --label {params.label} \
                --folds {params.folds} \
                --sizes {params.sizes} \
                --seed {params.seed}
            """


# ---------------------------------------------------------------- dashboard --
_DASH_INPUT = [
    f"{FIDELITY}/summary.json",
    f"{HYPOTHESIS}/summary.json",
    "dashboard/app.py",
]
if REVERSAL_ENABLED:
    _DASH_INPUT.append(f"{REVERSAL}/summary.json")
if config.get("impact", {}).get("enabled", False):
    _DASH_INPUT.append(f"{IMPACT}/summary.json")
if config.get("utility", {}).get("enabled", False):
    _DASH_INPUT.append(f"{UTILITY}/summary.json")
if SUBTYPE_ENABLED:
    _DASH_INPUT.append(f"{SUBTYPE_ROOT}/summary.json")


rule dashboard:
    """Record which artifacts the Streamlit app can show (the app is served
    separately with `streamlit run dashboard/app.py`)."""
    input:
        _DASH_INPUT,
    output:
        manifest=f"{DASHBOARD}/manifest.json",
    conda:
        "envs/py_env.yaml"
    params:
        results=OUT,
    shell:
        """
        python src/dashboard_manifest.py --results {params.results} --out {output.manifest}
        """


# ----------------------------------------------------------- scaling series --
if config.get("reversal", {}).get("poc", False):
    rule reversal_poc:
        """Validate reversal: known-effect recovery, null specificity, donor confound."""
        input:
            INGEST_FILES,
            f"{INGEST}/rna_counts_sel.mtx.gz",
            f"{INGEST}/atac_counts_sel.mtx.gz",
            f"{INGEST}/rna_features_sel.tsv",
            f"{INGEST}/atac_features_sel.tsv",
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{SPLIT}/split.csv",
            f"{CLUSTER}/obs_annotated.csv",
            f"{SIM}/manifest.csv",
        output:
            summary=f"{REVERSAL}/poc/summary.json",
        conda:
            "envs/py_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            split=f"{SPLIT}/split.csv",
            obs=f"{CLUSTER}/obs_annotated.csv",
            sim=SIM,
            out=f"{REVERSAL}/poc",
            label=OUT.rstrip("/").split("/")[-1],
            seeds=",".join(str(v) for v in config["simulation"]["seeds"]),
            target=config.get("reversal", {}).get("source", "PD"),
            reference=config.get("reversal", {}).get("target", "Aged"),
            min_cells=config.get("reversal", {}).get("poc_min_cells", 30),
            top_n=config.get("reversal", {}).get("poc_top_n", 100),
            permutations=config.get("reversal", {}).get("poc_permutations", 10),
            seed=config["split"]["seed"],
        shell:
            """
            python src/reversal_poc.py \
                --ingest {params.ingest} \
                --features {params.features} \
                --split {params.split} \
                --obs {params.obs} \
                --sim {params.sim} \
                --out {params.out} \
                --label {params.label} \
                --seeds {params.seeds} \
                --target {params.target} \
                --reference {params.reference} \
                --min-cells {params.min_cells} \
                --top-n {params.top_n} \
                --permutations {params.permutations} \
                --seed {params.seed}
            """


rule scaling_compare:
    """Merge every step's metrics into one long table (run after the steps)."""
    output:
        comparison=f"{SCALING.get('root', 'results/scale')}/comparison.csv",
        summary=f"{SCALING.get('root', 'results/scale')}/summary.json",
    conda:
        "envs/py_env.yaml"
    params:
        root=SCALING.get("root", "results/scale"),
        steps=",".join(str(s) for s in SCALING.get("steps", [])),
        fixed_step=SCALING.get("fixed_feature_step", 24),
    shell:
        """
        python src/compare_scaling.py \
            --root {params.root} \
            --steps {params.steps} \
            --fixed-step {params.fixed_step}
        """


# ------------------------------------------------------- subtype transfer ----
if SUBTYPE_ENABLED:
    rule subtype_split:
        """Donor-disjoint train/test donors inside each target cell type."""
        input:
            f"{CLUSTER}/obs_annotated.csv",
        output:
            donors=f"{SUBTYPE_ROOT}/donors.csv",
            dropped=f"{SUBTYPE_ROOT}/dropped.csv",
            summary=f"{SUBTYPE_ROOT}/split_summary.csv",
        conda:
            "envs/py_env.yaml"
        params:
            obs=f"{CLUSTER}/obs_annotated.csv",
            out=SUBTYPE_ROOT,
            celltypes=",".join(SUBTYPE_CELLTYPES),
            fraction=SUBTYPE.get("fraction", 0.30),
            min_cells_per_donor=SUBTYPE.get("min_cells_per_donor", 20),
            seed=SUBTYPE.get("seed", 20261112),
        shell:
            """
            python src/subtype_split.py \
                --obs {params.obs} \
                --out {params.out} \
                --celltypes {params.celltypes} \
                --fraction {params.fraction} \
                --min-cells-per-donor {params.min_cells_per_donor} \
                --seed {params.seed}
            """

    rule subtype_cluster:
        """k=2 subclustering per target cell type and modality."""
        input:
            INGEST_FILES,
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
            f"{CLUSTER}/obs_annotated.csv",
            f"{SUBTYPE_ROOT}/donors.csv",
        output:
            labels=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/labels.csv",
            embedding=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/embedding.csv",
            markers=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/markers.csv",
            features=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/features.csv",
            filtered=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/filtered_cells.csv",
            qc=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/qc.json",
        conda:
            "envs/r_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            obs=f"{CLUSTER}/obs_annotated.csv",
            donors=f"{SUBTYPE_ROOT}/donors.csv",
            out=SUBTYPE_ROOT,
            n_features=SUBTYPE.get("n_features", 1000),
            dims=SUBTYPE.get("dims", 30),
            seed=SUBTYPE.get("seed", 20261112),
            max_doublet_score=SUBTYPE.get("max_doublet_score", 0.3),
            min_cells_per_donor=SUBTYPE.get("min_cells_per_donor", 10),
            min_cells_per_subtype=SUBTYPE.get("min_cells_per_subtype", 50),
        shell:
            """
            Rscript R/06_subcluster.R \
                --ingest {params.ingest} \
                --features {params.features} \
                --obs {params.obs} \
                --donors {params.donors} \
                --out {params.out} \
                --celltype {wildcards.celltype} \
                --modality {wildcards.modality} \
                --n-features {params.n_features} \
                --dims {params.dims} \
                --seed {params.seed} \
                --max-doublet-score {params.max_doublet_score} \
                --min-cells-per-donor {params.min_cells_per_donor} \
                --min-cells-per-subtype {params.min_cells_per_subtype}
            """

    rule subtype_simulate:
        """scDesign3 alt / null / null-zero simulations on the subtype labels."""
        input:
            f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/labels.csv",
            f"{INGEST}/rna_counts_sel.mtx.gz",
            f"{INGEST}/atac_counts_sel.mtx.gz",
            f"{INGEST}/rna_features_sel.tsv",
            f"{INGEST}/atac_features_sel.tsv",
            f"{FEATURES}/rna_features.tsv",
            f"{FEATURES}/atac_features.tsv",
        output:
            manifest=f"{SUBTYPE_ROOT}/{{celltype}}/{{modality}}/simulation/manifest.csv",
            obs=expand(
                f"{SUBTYPE_ROOT}/{{{{celltype}}}}/{{{{modality}}}}/simulation/{{mode}}/obs.csv",
                mode=SUBTYPE_MODES),
            counts=expand(
                f"{SUBTYPE_ROOT}/{{{{celltype}}}}/{{{{modality}}}}/simulation/{{mode}}/seed{{seed}}_counts.mtx.gz",
                mode=SUBTYPE_MODES, seed=SUBTYPE_SEEDS),
        conda:
            "envs/r_env.yaml"
        params:
            ingest=INGEST,
            features=FEATURES,
            out=SUBTYPE_ROOT,
            family_use=SUBTYPE.get("family_use", "nb"),
            mu_formula_alt=SUBTYPE.get("mu_formula_alt", "subtype + condition"),
            seeds=",".join(str(s) for s in SUBTYPE_SEEDS),
            modes=",".join(SUBTYPE_MODES),
            max_cells=SUBTYPE.get("max_cells", 20000),
        shell:
            """
            Rscript R/07_subtype_sim.R \
                --ingest {params.ingest} \
                --features {params.features} \
                --labels {SUBTYPE_ROOT}/{wildcards.celltype}/{wildcards.modality}/labels.csv \
                --out {params.out}/{wildcards.celltype}/{wildcards.modality}/simulation \
                --celltype {wildcards.celltype} \
                --modality {wildcards.modality} \
                --family-use {params.family_use} \
                --mu-formula-alt "{params.mu_formula_alt}" \
                --seeds {params.seeds} \
                --modes {params.modes} \
                --max-cells {params.max_cells}
            """

    rule subtype_transfer:
        """Synthetic <-> real subtype transfer matrix with controls."""
        input:
            SUBTYPE_CLUSTER_FILES,
            SUBTYPE_SIM_FILES,
            f"{INGEST}/obs.csv",
            f"{INGEST}/rna_features_sel.tsv",
            f"{INGEST}/atac_features_sel.tsv",
        output:
            transfer=f"{SUBTYPE_TRANSFER}/transfer.csv",
            transfer_summary=f"{SUBTYPE_TRANSFER}/transfer_summary.csv",
            jaccard=f"{SUBTYPE_TRANSFER}/jaccard.csv",
            summary=f"{SUBTYPE_TRANSFER}/summary.json",
        conda:
            "envs/py_env.yaml"
        params:
            root=SUBTYPE_ROOT,
            ingest=INGEST,
            out=SUBTYPE_TRANSFER,
            celltypes=",".join(SUBTYPE_CELLTYPES),
            modalities=",".join(SUBTYPE_MODALITIES),
            seeds=",".join(str(s) for s in SUBTYPE_SEEDS),
            models=",".join(SUBTYPE.get("models", [])),
            max_train=SUBTYPE.get("max_train", 10000),
            max_test=SUBTYPE.get("max_test", 10000),
            min_stratum=SUBTYPE.get("min_stratum", 20),
        shell:
            """
            python src/subtype_transfer.py \
                --root {params.root} \
                --ingest {params.ingest} \
                --out {params.out} \
                --celltypes {params.celltypes} \
                --modalities {params.modalities} \
                --seeds {params.seeds} \
                --models {params.models} \
                --max-train {params.max_train} \
                --max-test {params.max_test} \
                --min-stratum {params.min_stratum}
            """

    rule subtype_compare:
        """Merge the subtype QC, concordance and transfer matrix."""
        input:
            f"{SUBTYPE_TRANSFER}/transfer.csv",
            f"{SUBTYPE_TRANSFER}/jaccard.csv",
            SUBTYPE_CLUSTER_FILES,
        output:
            comparison=f"{SUBTYPE_ROOT}/comparison.csv",
            summary=f"{SUBTYPE_ROOT}/summary.json",
        conda:
            "envs/py_env.yaml"
        params:
            root=SUBTYPE_ROOT,
        shell:
            """
            python src/subtype_compare.py --root {params.root}
            """
