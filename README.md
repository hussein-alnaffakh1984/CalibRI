# CalibRI

Code and summary results for the manuscript:

Interactive calculator: https://calibri-lab.netlify.app

> **Al-Naffakh HAH, et al.** Coverage audit and conformal recalibration of indirect reference-interval estimators (reflimR, RINet, refineR). *Submitted to Clinical Chemistry and Laboratory Medicine.*

The study measures how often the 95% intervals reported by three indirect reference-interval (RI) methods contain the true reference limit, and recalibrates them with split-conformal prediction. Benchmarks:
- RIbench-based simulated datasets.
- NHANES semi-synthetic datasets (pooled, five person-disjoint splits, healthy-only arm, n = 2,000 check).
- An external Korean benchmark (KNHANES 2013–2015).
- An application to TSH results from a laboratory in Al-Samawah, Iraq.

## Repository layout

| Folder | Contents |
|---|---|
| `code/01_data/` | Data preparation and dataset generators (simulation, NHANES, KNHANES, Samawah) |
| `code/02_estimators/` | refineR/reflimR worker (R), RINet inference with bootstrap (Python), run notebook, RINet model files |
| `code/03_analysis/` | Coverage, bias, widths, conformal recalibration, sensitivity and external analyses |
| `code/04_figures/` | Figures 1–4 (600 dpi) |
| `app/` | Browser calculator (reflimR and RINet ports in JavaScript, calibration factors, page source) |
| `code/05_checks/` | Check that the Python RINet implementation reproduces the rinet R package 0.1.1 |
| `data/nhanes/` | Download script for NHANES files; derived HbA1c file (NHANES 2011–2023) |
| `results/estimates/` | Per-dataset estimator outputs for the simulated and NHANES datasets |
| `results/bundle_meta/` | Dataset metadata, including the true limits, for the simulated and NHANES benchmarks |
| `results/tables/` | Summary results behind the tables and figures |

## Reproducing the results

All commands run from the repository root. Intermediate files go to `work/`, which is not tracked.

**A. Tables and figures from the provided estimator outputs (minutes)**

```bash
bash setup_work.sh
python code/03_analysis/disjoint.py      # person-disjoint coverage and conformal schemes (a)–(d)
python code/03_analysis/widths_d.py      # calibrated widths, scheme (d)
python code/03_analysis/bias.py
python code/03_analysis/sens.py          # five splits; RINet B = 50 vs 1,000
python code/03_analysis/asym.py          # signed and shape-group variants
python code/03_analysis/zero.py          # 0% contamination arm
python code/03_analysis/bin2000.py       # n = 2,000 calibration bin
python code/03_analysis/is_all.py        # interval scores
python code/03_analysis/supp_sim.py      # simulated benchmark
python code/03_analysis/knhanes_interval_scores.py   # KNHANES interval scores (S9)
python code/03_analysis/empirical_control.py         # empirical-percentile control (S7)
python code/03_analysis/refiner_bootstrap.py         # refineR bootstrap coverage
python code/04_figures/figs2.py          # writes figures/Figure1-4.png
```

We re-ran these scripts in this layout. The output tables matched the manuscript files, and the four figures were pixel-identical.

**B. Full pipeline (datasets → estimators → analyses)**

1. Simulated datasets: `Rscript code/01_data/gen_sim.R` (needs the RIbench package for its parameter grid).
2. NHANES thyroid data:
   - `bash data/nhanes/download_nhanes.sh`
   - `python code/01_data/prep_nhanes_thyroid.py`
3. NHANES datasets:
   - `python code/01_data/gen_nh.py`
   - `python code/01_data/merge_bundle1.py`
   - then `gen_nh2.py`, `gen_nh3.py` (split 1), `gen_nh4.py` (splits 2–5), `gen_nh5.py` (0% arm) and `gen_nh7.py` (n = 2,000).
   - Create each `work/ribundleN/` folder before running its generator.
4. Estimators:
   - `code/02_estimators/worker.R` (refineR + reflimR; `fast` and `boot` modes).
   - `rinet_runner.py` (RINet, B = 50).
   - `reflim_run.R`.
   - `run_all.sh` (RINet, B = 1,000).
   - `RUN_ri_audit.ipynb` runs the main parts on Kaggle or Colab.
5. Analyses and figures as in A.
6. Al-Samawah calibration (scheme d composition; needs the laboratory data): `python code/03_analysis/samawah_conf_d2.py d <samawah_dir> <out.csv>` (and `c` for the sensitivity analysis), then `gen_samawah_tables.py`. refineR on the n = 2,000 person-disjoint datasets: `Rscript code/02_estimators/worker.R fast work/ribundle7o work/res_b7o 1 1 <cores> <hours>` (bundle 7 restricted to total T4, total T3, thyroglobulin and HbA1c; `ribundle7t` = TSH).

Seeds are fixed in every generator (simulation and pooled NHANES 20261016; split 1 20261022). RINet bootstrap resamples use a per-dataset seed equal to `ds_id`. We regenerated the NHANES thyroid file and the split-1 datasets from the CDC files in this layout, and they were identical to those used in the paper.

## Data availability

- **NHANES** is public-use data from the U.S. CDC/NCHS. Thyroid files are downloaded by script. The HbA1c file is a derived NHANES 2011–2023 extract restricted to the variables used.
- **KNHANES** (Korea Disease Control and Prevention Agency) individual records may not be redistributed. Download them from the KNHANES website after registration (https://knhanes.kdca.go.kr), place `hn13_all.sas7bdat`, `hn14_all.sas7bdat` and `hn15_all.sas7bdat` in `work/knhanes/raw/`, then run `code/01_data/gen_kn.py` and `subset_kn_refiner.py`. Only the KNHANES true limits (`results/tables/knhanes_truth.csv`) and summary results are included here.
- **Samawah laboratory data** are not public. `code/01_data/prep_samawah.py` documents the cleaning steps, and only the summary estimates reported in the paper are included (`results/tables/samawah_*.csv`).

## Software

- **R 4.3.3:** refineR 2.0.0, reflimR 1.1.0, RIbench 1.0.2, data.table, future.
- **Python 3:** numpy, pandas, scipy, matplotlib, pyreadstat, keras with a TensorFlow backend.
- **RINet:** the rinet R package 0.1.1 was re-implemented in Python with the package's own model and scaler files (`code/02_estimators/rinet_models`, MIT licence, © 2025 Jack LeBien). In 60 test datasets the reference limits agreed with the R package to within a relative difference of 6 × 10⁻⁸ (`code/05_checks`).

## Licence

Code: MIT (see `LICENSE`). The RINet model files keep their original MIT licence. NHANES-derived files are public-domain U.S. government data.

## Contact

Hussein Ali Hussein Al-Naffakh, Al-Kafeel University, Najaf, Iraq (hussein.alnaffakh@alkafeel.edu.iq)
