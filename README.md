# Neighborhood Social Cohesion, Socioeconomic Disadvantage, and Incident Depression

Replication code for a prospective, multimodal cohort study in the **All of Us Research Program**
estimating the association between perceived neighborhood social cohesion and incident,
clinically diagnosed depression, its consistency across self-reported, wearable-measured, and
electronic-health-record (EHR) data, and its modification by socioeconomic disadvantage.

> Chen J, Larsson J. *Neighborhood Social Cohesion, Socioeconomic Disadvantage, and Incident Depression.* (manuscript)

---

## What this study does

- **Exposure:** perceived neighborhood social cohesion (5-item collective-efficacy scale, Social Determinants of Health survey), measured once at baseline.
- **Outcome:** incident depression — time from baseline to a first EHR diagnosis (155-concept SNOMED phenotype), among adults free of depression at baseline, right-censored at last EHR contact.
- **Design:** prospective Cox proportional-hazards models; primary estimate is the hazard ratio (HR) per 1-SD cohesion.
- **Multimodal triangulation:** the same exposure is related to self-reported affect (EHHWB survey), wearable-measured behavior (Fitbit), and the EHR diagnosis — data sources with independent measurement error.
- **Effect modification:** by individual income and area (ZIP3) deprivation.
- **Robustness:** 180-day lag period, E-values, representativeness raking, an antidepressant-initiation parallel outcome, and mediation analysis.

## Important: data access and environment

**No data are included in this repository**, and none can be. All analyses run **inside the
All of Us Researcher Workbench** (a secure cloud environment); individual-level All of Us data
cannot be downloaded or redistributed. To reproduce, you need:

1. An approved **All of Us Researcher Workbench** account with **Controlled Tier** access.
2. A **Controlled Tier workspace** (this study used curated data repository `C2025Q4R6`).
3. An **R** analysis environment (JupyterLab, R kernel) within that workspace.

The scripts read the workspace's BigQuery dataset via the environment variables
`WORKSPACE_CDR` and `GOOGLE_PROJECT`, which the Workbench sets automatically.

Only **aggregate** results (hazard ratios, coefficients, counts ≥ 20) leave the environment,
per the All of Us Data and Statistics Dissemination Policy.

## Requirements

R (≥ 4.0) with:

```r
install.packages(c("tidyverse", "bigrquery", "survival", "sandwich", "lmtest", "EValue", "splines", "ggplot2"))
```

`splines` and `ggplot2` ship with R / tidyverse. `bigrquery` authenticates automatically inside the Workbench.

## Repository structure

```
R/
  00_discovery/   Concept-ID and feasibility discovery (run first, informational)
  01_build/       Construct analysis datasets (BigQuery -> .rds)
  02_models/      Primary and multimodal models
  03_robustness/  Sensitivity analyses
  04_outputs/     Tables and figures
docs/
  TECHNICAL_DOCUMENTATION.md   Step-by-step explanation of every script
  DATA_DICTIONARY.md           Concept IDs, variables, and derived measures
```

## Pipeline (run order)

Datasets are cached as `.rds` files in the working directory; models read them. Run in this order:

| Step | Script | Produces | Needs |
|---|---|---|---|
| 1 | `01_build/build_cox.R` | `cox_dat.rds` (main at-risk cohort) | — |
| 2 | `01_build/build_cox_swb.R` | `swb_dat.rds` (self-report subcohort) | — |
| 3 | `01_build/build_wb_items.R` | `wb_item_dat.rds` (item-level affect) | — |
| 4 | `01_build/build_fitbit.R` | `fitbit_dat.rds` (wearable) | `cox_dat.rds` |
| 5 | `01_build/build_antidep.R` | `antidep_dat.rds` (parallel outcome) | `cox_dat.rds` |
| 7 | `02_models/*` | printed results | the relevant `.rds` |
| 8 | `03_robustness/*` | printed results | the relevant `.rds` |
| 9 | `04_outputs/figs_data.R` -> `figs_plot.R` | `fig_estimates.csv`, figures | `.rds` files |
| 10 | `04_outputs/table1.R`, `table2.R` | `table1.csv`, `table2.csv` | `.rds` files |

See `docs/TECHNICAL_DOCUMENTATION.md` for exactly what each script does, its inputs/outputs, and the methods it implements.

## Key result

Higher perceived cohesion was associated with lower hazard of incident depression
(HR, 0.88 per SD; 95% CI, 0.86-0.91), consistent across self-report and wearable channels,
and weaker among low-income individuals and in more-deprived areas.

## Acknowledgment

This research used data from the All of Us Research Program, supported by the National
Institutes of Health (see the manuscript for the full award list). The program would not be
possible without the partnership of its participants.

## License

Code released under the MIT License (see `LICENSE`). All of Us data are governed separately by
the program's Data Use and Registration Agreement.
