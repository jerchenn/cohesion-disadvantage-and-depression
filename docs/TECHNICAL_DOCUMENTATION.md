# Technical Documentation

Detailed, script-by-script explanation of the replication pipeline. For each script: its
**purpose**, **inputs**, **outputs**, and a **step-by-step** account of what the code does and
why. All scripts run inside the All of Us Researcher Workbench (Controlled Tier) against the
workspace BigQuery dataset (`WORKSPACE_CDR`, `GOOGLE_PROJECT`).

Conventions used throughout:
- `run_sql(q)` submits a BigQuery query and downloads the result with `bigint = "character"`
  (person_id is INT64 and would overflow R's 32-bit integer to NA otherwise).
- Survey items are pulled from `ds_survey`; diagnoses from `condition_occurrence`; drugs from
  `drug_exposure`; labs/vitals from `measurement`; demographics from `person`; concept metadata
  and hierarchies from `concept` and `concept_ancestor`.
- Datasets are cached as `.rds` in the working directory so models can be re-run without re-querying.

---

## 0. Discovery (`R/00_discovery/`)

These are informational: they identify concept IDs and check feasibility. Not required to run
the pipeline once the IDs (baked into the build scripts) are known, but included for transparency.

### `phenotype_cohort.R`
**Purpose:** derive the depression outcome phenotype and quantify the at-risk cohort.
**What it does:** starting from four SNOMED parent concepts (Depressive disorder, Major
depressive disorder, Recurrent depressive disorder, Dysthymia), it walks `concept_ancestor` to
collect all **standard, Condition-domain descendant concepts**, excluding any whose name
contains "bipolar". This yields the 155-concept depression phenotype. It then counts how many
SDOH respondents are EHR-linked, depression-free at baseline, and thus at risk, and how many
develop incident depression.

### `wb_items_discover.R`
**Purpose:** enumerate the Emotional Health History and Well-Being (EHHWB) survey items.
**What it does:** lists every `question_concept_id` in that survey with respondent counts and full
answer-option distributions, so the affect items and their response scales can be coded correctly.

### `fitbit_scout.R`
**Purpose:** check wearable data feasibility.
**What it does:** lists Fitbit-related tables, their per-person coverage and date spans, the overlap
with the cohesion cohort, and the columns available in the daily-summary tables.

### `biomarker_check.R`
**Purpose:** assess biomarker coverage and selection for the exploratory biology layer.
**What it does:** counts coverage of candidate biomarkers (e.g., HbA1c, CRP) and profiles how
biomarker-havers differ (age, prior depression, utilization) — establishing that lab availability
is utilization-driven and therefore selected.

---

## 1. Dataset construction (`R/01_build/`)

### `build_cox.R`  ->  `cox_dat.rds`  (the foundational dataset)
**Purpose:** assemble the prospective survival dataset: cohesion exposure + covariates +
time-to-incident-depression, on the SDOH∩EHR at-risk cohort.
**Steps:**
1. **Exposure + SES items** — pull the 5 cohesion items (help, get-along, trust, values,
   watch-out) plus income, education, employment from `ds_survey`; deduplicate on
   (person_id, question); reshape wide.
2. **Demographics** — age, sex, race, ethnicity from `person` joined to `concept`.
3. **Score cohesion** — recode each item to numeric with scale-aware maps (5-point agree scale;
   watch-out is 4-point); treat "Skip / Prefer not to answer / Don't know" as missing; z-score
   each item, average (requiring ≥60% items answered), re-standardize to give `z_cohesion`
   (mean 0, SD 1).
4. **SES** — map income and education to ordinal numerics; build missing-indicator factors
   (`income_f`, `educ_f`, `emp_f`) via `addNA` so non-response is retained as its own level.
5. **Survival outcome** — with three CTEs: `expo` (baseline = earliest date of the cohesion item
   per person), `dep` (first depression diagnosis, using the 155-concept phenotype), `ehr`
   (first/last condition date and total condition count). Keep participants who are EHR-observed
   around baseline (≥1 year prior lookback, ≥1 record after) and **depression-free before
   baseline**. Define `event` = 1 if a depression diagnosis occurs after baseline; `time_days` =
   baseline→first diagnosis (events) or baseline→last EHR record (censored).
6. **Utilization** — `log_util = log1p(n_cond)`, a proxy for ascertainment intensity.
7. Save `cox_dat.rds`.

### `build_cox_swb.R`  ->  `swb_dat.rds`
**Purpose:** the self-reported-wellbeing subcohort (SDOH∩EHHWB∩EHR), baseline = EHHWB survey
date, to test cohesion → subjective wellbeing → depression.
**Steps:** as `build_cox.R`, but additionally pulls the happiness (hedonic) and life-meaning
(eudaimonic) items, scores them into a standardized wellbeing index `z_wb`, and sets baseline to
the EHHWB date so the affect measure is contemporaneous with baseline.

### `build_wb_items.R`  ->  `wb_item_dat.rds`
**Purpose:** the item-level affect dataset — every individual environment item and a curated set
of EHHWB affect items, each oriented so higher = better mental state, with the survival outcome.
**Steps:** pull all environment items (cohesion, safety, disorder, walkability) and affect items
(happiness, meaning, social connection, plus PHQ-9 and GAD-7 symptom items tagged by
construct-overlap tier); recode with scale-aware maps; reverse-score symptoms so higher = fewer
symptoms; z-score each; attach the incident-depression survival outcome (baseline = happiness-item
date). Used for the item-resolution analyses.

### `build_fitbit.R`  ->  `fitbit_dat.rds`
**Purpose:** attach wearable-measured behavior to the cohesion cohort.
**Steps:** aggregate the Fitbit **daily-summary** tables (never the intraday tables — billions of
rows) to person-level means: mean daily **steps**, sedentary minutes, moderate-to-vigorous
activity (MVPA), over wear-days ≥ 10; mean **minutes asleep**, **sleep efficiency**
(asleep/in-bed), and **wake-after-sleep-onset (WASO)** over main-sleep nights ≥ 5. Merge onto
`cox_dat.rds`.

### `build_antidep.R`  ->  `antidep_dat.rds`
**Purpose:** an independent, prescribing-based parallel outcome.
**Steps:** define antidepressants as standard Drug descendants of ATC class `N06A`; among
participants **antidepressant-naive at baseline**, define time to first antidepressant
initiation; reuse cohesion + covariates from `cox_dat.rds`; fit the same Cox model. Corroborates
the diagnosis outcome against coding idiosyncrasies.

---

## 2. Models (`R/02_models/`)

### `cox_model.R`  (primary result)
Fits the primary Cox model `Surv(time_days, event) ~ z_cohesion + covariates` on `cox_dat.rds`,
adjusting for age, sex, race/ethnicity, income, education, employment, and utilization. Reports
the HR per SD; repeats with a **180-day lag period** (excluding early events to guard against
reverse causation); and tests the **proportional-hazards** assumption via scaled Schoenfeld
residuals.

### `cox_swb_model.R`
On `swb_dat.rds`: (Test 2) self-reported wellbeing → incident depression, full and 180-day lag;
(Test 1) cohesion → wellbeing → depression via the difference method (compare the cohesion HR with
and without the wellbeing mediator), reporting an approximate proportion mediated. Uses collapsed
covariate factors to avoid sparse-cell separation.

### `wb_items_cox.R`
On `wb_item_dat.rds`: for each affect item, a Cox HR per SD on incident depression, tagged by
construct-overlap tier, full and lag-period; plus a joint model of the low-overlap items
(happiness, meaning, social connection) to test independent contributions.

### `fitbit_model.R`
On `fitbit_dat.rds`, two directions: (1) cohesion → device behavior (linear models, HC3 SEs);
(2) device behavior → incident depression (Cox, per SD oriented so healthier = protective), full
and lag-period; plus a joint steps + sleep-efficiency model. Uses collapsed covariates.

### `sleep_headtohead.R`
Gating check for the sleep dissociation: within the intersection of the self-report and Fitbit
subcohorts, models self-reported sleep and device sleep **in the same people**, alone and jointly,
to test whether the self-reported signal is reducible to objective sleep.

### `env_wb_linkage.R`
On `wb_item_dat.rds`: (1) a matrix of standardized associations between each individual
environment item and each affect item (adjusted for age and income); (2) a "bridge" test —
whether the affect channels most coupled to cohesion are the same ones that most predict
depression (Spearman rank correlation of cohesion-coupling vs depression-protection across items).

---

## 3. Robustness (`R/03_robustness/`)

### `evalue.R`
Computes E-values (VanderWeele-Ding) for the key hazard ratios: the minimum association a single
unmeasured confounder would need with both exposure and outcome to explain away each estimate.
Requires the `EValue` package.

### `rake_ipw.R`
Representativeness sensitivity: rakes the analytic sample to approximate US adult margins (age,
sex, race/ethnicity, education) by iterative proportional fitting, then re-estimates the primary
HR weighted, with a bootstrap CI. (Note: `coxph` auto-enables robust SEs with non-integer weights
and can hit a residuals bug; the script uses `robust = FALSE` and a bootstrap instead.)

### `mediation.R`
Difference-method mediation of cohesion → incident depression through (a) objective physical
activity (steps) and (b) self-reported affect, on the log-HR scale, with the non-collapsibility and
cross-sectional-mediator caveats noted.

### `hetero_cox.R`
Effect modification and shape on the incident-depression outcome: cohesion quartile hazard ratios
(dose-response shape), a quadratic and spline nonlinearity test, and interaction tests /
subgroup HRs by income and age. Drops the collinear SES covariate within each modifier's model.

### `geo_hetero.R`
Area-level effect modification: joins ZIP3 socioeconomic data (`ds_zip_code_socioeconomic`,
`deprivation_index`), tests cohesion × area-deprivation interaction, and reports cohesion HRs by
deprivation tertile.

---

## 4. Outputs (`R/04_outputs/`)

### `table1.R`  ->  `table1.csv`, `etable1.csv`
Cohort characteristics by cohesion tertile (Table 1) and the subcohort comparison
(main / self-report / wearable; eTable 1), formatted as No. (%) and mean (SD) / median (IQR).

### `table2.R`  ->  `table2.csv`
Primary + sensitivity analyses (primary, 180-day lag, reweighted, antidepressant parallel outcome)
recomputed under the primary specification so the numbers match the manuscript; E-value footnote.

### `figs_data.R`  ->  `fig_estimates.csv`
Recomputes every figure estimate from the cached `.rds` files (avoiding stale hard-coded values)
into one tidy CSV: primary + lag, income and area-deprivation modification (with interaction P),
and the multimodal concordance estimates.

### `figs_plot.R`  ->  `figure1_modification.{png,pdf}`, `figure2_concordance.{png,pdf}`
Reads `fig_estimates.csv` and draws Figure 1 (effect modification, interaction P in-panel) and
Figure 2 (multimodal concordance), on a log-HR axis with a reference line at 1.0, grayscale-safe.
