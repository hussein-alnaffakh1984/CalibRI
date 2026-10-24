#!/bin/bash
# Creates work/ and fills it with the provided estimator outputs and dataset metadata, so that the
# analysis scripts (code/03_analysis) and figures (code/04_figures) can be re-run without re-running the estimators.
# Run from the repository root:  bash setup_work.sh
set -e
mkdir -p work/ana figures
cp -r results/estimates/* work/
cp -r results/bundle_meta/* work/
echo "work/ is ready"
