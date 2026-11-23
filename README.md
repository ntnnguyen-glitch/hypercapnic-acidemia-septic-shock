# Hypercapnic acidemia and 28-day mortality in septic shock (MIMIC-IV)

R code for the retrospective cohort study "Hypercapnic Acidemia and 28-Day Mortality and Vasopressor Dependence in Septic Shock: A Retrospective Cohort Study Based on the MIMIC-IV Database".

## Data

The study uses MIMIC-IV version 3.1 (https://doi.org/10.13026/kpb9-mt58). Access requires credentialing on PhysioNet, completion of the required training and signature of the data use agreement. No patient data are included in this repository.

## Requirements

R (4.3 or later) with the packages DBI, duckdb, data.table, survival, cmprsk, splines, mice, pROC, EValue, WeightIt, MatchIt, ggplot2, patchwork, scales, jsonlite and ragg.

## How to run

1. Download the MIMIC-IV 3.1 files (`hosp/` and `icu/` folders with `.csv.gz` files).
2. Set the path to the data folder:
   ```sh
   export MIMIC_DIR=/path/to/mimic-iv-3.1
   ```
3. From the repository folder, run:
   ```sh
   Rscript run_all.R
   ```

## Scripts

| Script | Content |
|---|---|
| `R/00_setup.R` | Paths and packages |
| `R/01_extract.R` | Extraction of the study variables from the raw files into a local DuckDB file |
| `R/02_cohort.R` | Cohort selection, exposure definition (first arterial blood gas), covariates and outcomes |
| `R/03_analysis.R` | Logistic and Cox regression, restricted cubic splines, secondary outcomes, propensity score analyses, sensitivity and subgroup analyses |
| `R/04_tables_figures.R` | Table 1 and Figures 1-4 |

Outputs are written to `results/` (tables as CSV, `results.rds`, figures as TIFF, PDF and PNG). Intermediate files are written to `work/`.

## Exposure groups

Based on the first arterial blood gas between 6 h before and 24 h after ICU admission: hypocapnia (PaCO2 <35 mmHg), normocapnia (35-45 mmHg, reference), compensated hypercapnia (PaCO2 >45 mmHg, pH >=7.35) and hypercapnic acidemia (PaCO2 >45 mmHg, pH <7.35).

## License

MIT License (see `LICENSE`).

## Contact

Ngoc-Nguyen Thi Nguyen, Faculty of International Medicine, Nam Can Tho University, Can Tho City, Vietnam (ntnnguyen@nctu.edu.vn).
