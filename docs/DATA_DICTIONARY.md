# Data Dictionary

Concept IDs, source tables, and derived variables used in the pipeline. All IDs are OMOP
`concept_id`s in the All of Us Controlled Tier (CDR `C2025Q4R6`). Verify against the current CDR
if reproducing on a later release.

## Exposure — perceived neighborhood social cohesion (Social Determinants of Health survey)

| Item (label) | question_concept_id | Scale |
|---|---|---|
| Neighbors willing to help (`help`) | 40192463 | 5-pt agree |
| Neighbors get along (`getalong`) | 40192411 | 5-pt agree |
| Neighbors can be trusted (`trust`) | 40192499 | 5-pt agree |
| Neighbors share values (`values`) | 40192417 | 5-pt agree |
| Neighbors watch out for each other (`watchout`) | 40192400 | 4-pt agree |

Scored: recode to numeric, z-score each item, average (≥60% items required), re-standardize
→ `z_cohesion` (mean 0, SD 1). Cronbach α = 0.88.

## Other environment items (SDOH; used in `env_wb_linkage.R`)

Safety/danger: safe 40192384, crime 40192493, unsafe-night 40192492, unsafe-day 40192414.
Disorder: graffiti 40192420, noisy 40192522, vandalism 40192412, abandoned 40192469,
hang-around 40192500, drug 40192457, alcohol 40192476, trouble 40192404, clean 40192456,
upkeep 40192386. Walkability: shops 40192436, transit 40192440, sidewalks 40192437,
bike 40192431, recreation 40192410.

## Outcome — incident depression (EHR)

- Source: `condition_occurrence`.
- Phenotype: **155 standard SNOMED depressive-disorder concepts** = Condition-domain descendants
  (via `concept_ancestor`) of `Depressive disorder`, `Major depressive disorder`, `Recurrent
  depressive disorder`, `Dysthymia`; excluding any concept whose name contains "bipolar".
- Derived: `event` (1 = first diagnosis after baseline), `time_days` (baseline→first diagnosis, or
  baseline→last EHR record if censored).

## Parallel outcome — antidepressant initiation

- Source: `drug_exposure`; standard Drug descendants of ATC class `N06A`.
- Cohort: antidepressant-naive at baseline; `event`/`time_days` to first initiation.

## Affect items (EHHWB survey)

| Item | question_concept_id | Tier |
|---|---|---|
| Happiness (hedonic) | 1703980 | wellbeing (low overlap) |
| Life meaning (eudaimonic) | 1704001 | wellbeing (low overlap) |
| Felt distant/cut off (social connection) | 1703997 | social |
| PHQ-9 items (down, anhedonia, worthlessness, fatigue, sleep, appetite, concentration, psychomotor×2, suicidal) | 1704024, 1704026, 1704041, 1704039, 1703983, 1704004, 1704038, 1703996, 1704028, 1703977 | depressive (high overlap) |
| GAD-7 items (nervous, worry-stop, worry-much, relax, afraid, irritable) | 1703984, 1703995, 1704000, 1703987, 1704042, 1703920 | anxiety |

Orientation: positive-wellbeing items keep their scale; symptom items are reverse-scored so higher
= better mental state. Each z-scored.

## SES / demographics

Income 1585375, education 1585940, employment 1585952 (Basics survey); age/sex/race/ethnicity
from `person`. Missing-indicator factors `income_f`, `educ_f`, `emp_f` (via `addNA`) in the
primary spec; median-filled continuous `income_m`, `educ_m` in the parsimonious spec used for
subgroup/mediation models.

## Wearable metrics (Fitbit daily-summary tables)

- `activity_summary`: `steps`, `sedentary_minutes`, `very_active_minutes`+`fairly_active_minutes`
  (→ MVPA); person-level means over wear-days with 100 ≤ steps ≤ 45000, ≥10 days.
- `sleep_daily_summary` (main sleep, 120 ≤ minute_in_bed ≤ 900, ≥5 nights): `minute_asleep`,
  efficiency = `minute_asleep`/`minute_in_bed`, WASO = `minute_awake`+`minute_restless`.

## Biomarkers (allostatic load; `measurement`)

Concepts discovered at runtime by name (top standard concept by coverage): systolic/diastolic BP,
heart rate, HbA1c, total cholesterol, HDL, triglycerides, C-reactive protein, BMI. AL = count of
biomarkers in the high-risk quartile (bottom quartile for HDL), ≥6 systems required.

## Geography

`ds_zip_code_socioeconomic` (person-level, ZIP3): `deprivation_index` (ADI-style composite) used
for area-deprivation effect modification.

## Key derived analysis variables

| Variable | Meaning |
|---|---|
| `z_cohesion` | standardized cohesion exposure (per-SD unit of all HRs) |
| `event`, `time_days` | survival outcome (see above) |
| `log_util` | log1p(count of EHR condition records); ascertainment proxy |
| `z_wb` | standardized happiness+meaning wellbeing index |
| `z_steps`, `z_eff`, `z_waso` | standardized wearable metrics (oriented healthier = higher) |
| `al` | allostatic load score (0–9 systems dysregulated) |
