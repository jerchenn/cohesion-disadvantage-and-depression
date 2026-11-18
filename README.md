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

## Replication from a fresh notebook

Open a new R notebook in an All of Us Controlled Tier workspace and run one cell:

```r
if (!dir.exists("~/cohesion-disadvantage-and-depression"))
  system("git clone https://github.com/jerchenn/cohesion-disadvantage-and-depression.git ~/cohesion-disadvantage-and-depression")
setwd("~/cohesion-disadvantage-and-depression"); system("git pull")
source("R/run_all.R")
```

`R/run_all.R` installs any missing packages (survival, bigrquery, tidyverse, sandwich, lmtest, gridExtra, EValue),
rebuilds every analysis dataset from BigQuery, and runs every model, table, and figure in dependency order. Each
script runs in its own environment; a script that fails is logged and the run continues, and the failures are
listed at the end. To resume, set `START_AT <- "<path of the failed script>"` and source `R/run_all.R` again.
All printed output is aggregate (counts under 20 masked) and is copied to `replication_log.txt`. The full run
takes several hours.

## Repository structure

```
R/
  run_all.R       Replication script (the only entry point needed)
  00_discovery/   Concept discovery and cohort checks (verify_recodes.R and sanity_checks.R run in run_all.R)
  01_build/       Analysis datasets (BigQuery -> .rds): build_cox3.R is the primary cohort
  02_models/      Primary and subcohort models
  03_robustness/  Sensitivity analyses in the original submission
  04_outputs/     Tables and figures (fig_forest.R draws Figures 2 and 3)
  07_revision/    Sensitivity analyses added at revision (see R/07_revision/README.md)
docs/
  TECHNICAL_DOCUMENTATION.md   Step-by-step explanation of every script
  DATA_DICTIONARY.md           Concept IDs, variables, and derived measures
```

Scripts in 00_discovery other than the two checks, 05_recognition, 06_*, and the older build_cox.R and cox_model.R are
exploratory or superseded and are not used by any reported result.

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
