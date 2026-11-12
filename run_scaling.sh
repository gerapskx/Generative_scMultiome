#!/usr/bin/env bash
# Run the Aged/PD sample-size series after the one-time full ingest.
#
#   bash run_scaling.sh              # all steps, sequential
#   bash run_scaling.sh n24 n08      # selected steps, in the given order
#
# n24 must run first (or at least its `features` rule) because every smaller
# step reuses `results/scale/n24/features` as a fixed feature set.
#
# Each step is a separate Snakemake invocation with its own config, so a failed
# step can be resumed without touching the others.
set -euo pipefail

cd "$(dirname "$0")"

SM="$HOME/anaconda3/envs/multiome-py/bin/snakemake"
FULL="results/scale/full/ingest/obs.csv"
ROOT="results/scale"
JOBS="${JOBS:-1}"

if [[ ! -f "$FULL" ]]; then
    echo "full ingest is not finished: $FULL is missing" >&2
    echo "run: snakemake --configfile config/scaling/full.yaml -j 4 --use-conda" >&2
    exit 1
fi

steps=("$@")
if [[ ${#steps[@]} -eq 0 ]]; then
    steps=(24 08 16 02)
fi

for step in "${steps[@]}"; do
    tag="n${step}"
    cfg="config/scaling/${tag}.yaml"
    if [[ ! -f "$cfg" ]]; then
        echo "missing $cfg (run: python src/plan_scaling.py)" >&2
        exit 1
    fi
    echo "=== ${tag} $(date -Is) ==="
    # Clear a stale lock left behind by a killed run (safe: steps are sequential
    # and we only get here when no other Snakemake for this config is alive).
    if ! pgrep -f "snakemake --configfile ${cfg}" >/dev/null 2>&1; then
        "$SM" --configfile "$cfg" --unlock >/dev/null 2>&1 || true
    fi
    "$SM" --configfile "$cfg" -j "$JOBS" --use-conda --conda-frontend conda \
        --rerun-incomplete > "logs/scale_${tag}.log" 2>&1
    echo "=== ${tag} done $(date -Is) ==="
done

echo "=== comparison $(date -Is) ==="
~/anaconda3/envs/multiome-py/bin/python src/compare_scaling.py \
    --root "$ROOT" --steps 2,8,16,24 --fixed-step 24
echo "=== all done $(date -Is) ==="
