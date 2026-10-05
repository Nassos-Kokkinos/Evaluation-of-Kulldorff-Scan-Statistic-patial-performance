# Evaluation of Kulldorff's Scan Statistic — spatial performance

R code accompanying my undergraduate thesis at the Department of Mathematics,
Aristotle University of Thessaloniki (supervisor: Georgios Tsaklidis, February 2026).

The thesis studies how precisely the Kulldorff spatial scan statistic
(as implemented in [SaTScan™](https://www.satscan.org)) *localises* clusters,
beyond simply detecting that they exist.

## Overview

| Script | Thesis section | What it does |
|---|---|---|
| [`R/01_nys_roc_pr_analysis.R`](R/01_nys_roc_pr_analysis.R) | Ch. 6.1, Appendix A | NYS breast cancer case study. Uses SaTScan's per-location relative risk as a score and evaluates each of the 5 detected clusters with ROC/AUC, Youden-optimal threshold and Precision–Recall/AUPRC. |
| [`R/02_simulation_spatial_overlap.R`](R/02_simulation_spatial_overlap.R) | Ch. 6.2, Appendix B | Synthetic 20 × 20 lattice with known high- and low-rate clusters. Compares detected vs. true regions with sensitivity, precision, Jaccard index and symmetric difference area. |

### SaTScan™ configuration

| Parameter | Setting |
|---|---|
| Version | SaTScan™ v10.3.3 |
| Model | Discrete Poisson |
| Analysis | Purely spatial |
| Maximum cluster size | 50% of population at risk |
| Monte Carlo replications | 999 |

## Key findings

- SaTScan detects the high-risk NYS clusters with very small Monte Carlo
  p-values, yet per-location discrimination is only **moderate**: relative risks
  inside and outside the clusters overlap substantially.
- In the simulation, a strong compact **high-rate** cluster is recovered almost
  perfectly, while the **low-rate** cluster is captured as a more diffuse region
  that covers the true core plus surrounding cells.
- Statistical significance does not guarantee precise spatial delineation of
  cluster boundaries.

## Repository structure

```
.
├── R/
│   ├── 01_nys_roc_pr_analysis.R
│   └── 02_simulation_spatial_overlap.R
├── CITATION.cff
├── LICENSE
└── README.md
```

## Input data

The data files are not included. Create the folders below in the repository
root and place the files there before running the scripts.

**`data/nys/`** — NYS breast cancer case study

| File | Description |
|---|---|
| `NYS_BreastCancer.geo` | SaTScan™ sample data — LOC_ID, LAT, LON |
| `NYS_BreastCancer_results.rr.dbf.xlsx` | SaTScan™ per-location relative risk (`*.rr.dbf`) exported to Excel; columns `LOC_ID`, `REL_RISK` |
| `cluster1_locations.kml` … `cluster5_locations.kml` | SaTScan™ KML output, one file per detected cluster |

**`data/simulation/`** — 20 × 20 lattice simulation

| File | Description |
|---|---|
| `sim_example_true_labels.csv` | `;`-separated; columns `loc_id`, `in_high_true` (0/1), `in_low_true` (0/1) |
| `sim_example_results.gis.dbf` | SaTScan™ GIS output; columns `LOC_ID`, `CLUSTER`, `LOC_RR` |

The NYS `.cas`, `.pop` and `.geo` inputs ship with SaTScan™ as sample data
(<https://www.satscan.org>).

## Usage

1. Install [R](https://cran.r-project.org/) (≥ 4.1).
2. Run SaTScan™ with the configuration above and place the files in `data/`
   as described in [Input data](#input-data).
3. From the repository root:

```bash
Rscript R/01_nys_roc_pr_analysis.R
Rscript R/02_simulation_spatial_overlap.R
```

Missing packages are installed automatically on first run
(`readxl`, `dplyr`, `stringr`, `ggplot2`, `pROC`, `PRROC`, `purrr`, `tibble`,
`readr`, `RANN`, `foreign`).

Results are written (folders are created automatically) to `outputs/nys/` (ROC/PR plots, CSV and LaTeX summary
table) and `outputs/simulation/` (overlap metrics CSV).

## Citation

> Kokkinos, A. (2026). *Evaluation of Kulldorff's Scan Statistic spatial
> performance*. Undergraduate thesis, Department of Mathematics, Aristotle
> University of Thessaloniki.

## References

- Kulldorff, M. (1997). A spatial scan statistic. *Communications in Statistics – Theory and Methods*, 26(6), 1481–1496.
- Kulldorff, M. (2022). *SaTScan™ User Guide*, v10.x.

## License

[MIT](LICENSE)
