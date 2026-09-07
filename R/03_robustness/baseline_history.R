# baseline_history.R -- sensitivity: adjust the primary Cox for pre-baseline mental-health history.
#
# WHY: the exposure is a single self-reported perception, so a depressive/anxious disposition at
# baseline could lower perceived cohesion AND raise later diagnosis (confounding). The at-risk
# cohort is already free of DEPRESSION before baseline, so the informative pre-baseline markers are
# (a) a prior ANXIETY diagnosis and (b) a prior PSYCHOTROPIC prescription. Adding both to the
# primary model tests whether the cohesion association survives adjustment for baseline mental
# health. Report the z_cohesion HR before vs after adjustment.
#
# INPUT : cox_dat.rds (from build_cox.R)     OUTPUT: printed HRs
# RUN AFTER: build_cox.R.   Needs: tidyverse, survival, bigrquery.

library(tidyverse); library(survival); library(bigrquery)
cdr  <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

cox <- readRDS("cox_dat3.rds")

## 1. re-derive each person's baseline date = earliest cohesion-item survey date (as in build_cox.R).
##    (ds_survey date column is survey_datetime; adjust here if your CDR names it differently.)
bl <- run_sql(sprintf("
  SELECT person_id, MIN(DATE(survey_datetime)) AS bl_date
  FROM `%s.ds_survey`
  WHERE question_concept_id IN (40192463,40192411,40192499,40192417,40192400)
  GROUP BY person_id", cdr))
bl$bl_date <- as.Date(bl$bl_date)

## 2. pre-baseline ANXIETY diagnosis: standard SNOMED Condition descendants of 'Anxiety disorder'
##    (captures GAD, panic, phobias). Depression is already excluded from the cohort.
anx <- run_sql(sprintf("
  WITH parents AS (
    SELECT concept_id FROM `%s.concept`
    WHERE vocabulary_id='SNOMED' AND standard_concept='S' AND domain_id='Condition'
      AND LOWER(concept_name)='anxiety disorder'),
  ids AS (
    SELECT DISTINCT ca.descendant_concept_id AS concept_id
    FROM `%s.concept_ancestor` ca JOIN parents p ON ca.ancestor_concept_id=p.concept_id)
  SELECT co.person_id, MIN(DATE(co.condition_start_date)) AS first_anx
  FROM `%s.condition_occurrence` co JOIN ids i ON co.condition_concept_id=i.concept_id
  GROUP BY co.person_id", cdr, cdr, cdr))
anx$first_anx <- as.Date(anx$first_anx)

## 3. pre-baseline PSYCHOTROPIC prescription: standard Drug descendants of ATC N05A/N05B/N06A
##    (antipsychotics, anxiolytics, antidepressants) -- a proxy for treated baseline mental illness.
med <- run_sql(sprintf("
  WITH atc AS (
    SELECT concept_id FROM `%s.concept`
    WHERE vocabulary_id='ATC' AND concept_code IN ('N05A','N05B','N06A')),
  ids AS (
    SELECT DISTINCT ca.descendant_concept_id AS concept_id
    FROM `%s.concept_ancestor` ca JOIN atc a ON ca.ancestor_concept_id=a.concept_id)
  SELECT de.person_id, MIN(DATE(de.drug_exposure_start_date)) AS first_med
  FROM `%s.drug_exposure` de JOIN ids i ON de.drug_concept_id=i.concept_id
  GROUP BY de.person_id", cdr, cdr, cdr))
med$first_med <- as.Date(med$first_med)

## 4. build pre-baseline indicators on the analytic cohort (default 0 when no prior record exists)
cox <- cox %>% left_join(bl, by="person_id") %>%
               left_join(anx, by="person_id") %>%
               left_join(med, by="person_id")
cox$prior_anx       <- as.integer(!is.na(cox$first_anx) & cox$first_anx < cox$bl_date)
cox$prior_psych_med <- as.integer(!is.na(cox$first_med) & cox$first_med < cox$bl_date)
cat(sprintf("Pre-baseline anxiety dx: %d (%.1f%%) | pre-baseline psychotropic Rx: %d (%.1f%%)\n",
    sum(cox$prior_anx), 100*mean(cox$prior_anx),
    sum(cox$prior_psych_med), 100*mean(cox$prior_psych_med)))

## 5. primary model: unadjusted vs adjusted for baseline mental-health history
f0 <- Surv(time_days,event) ~ z_cohesion + age + sex + race + ethnicity +
        income_f + educ_f + emp_f + log_util
m0 <- coxph(f0, cox)
m1 <- coxph(update(f0, . ~ . + prior_anx + prior_psych_med), cox)
hr <- function(m) sprintf("%.3f (%.3f-%.3f)", exp(coef(m)[["z_cohesion"]]),
        exp(confint(m)["z_cohesion",1]), exp(confint(m)["z_cohesion",2]))

cat("\n=== Primary z_cohesion HR ===\n")
cat("  without baseline-MH adjustment : ", hr(m0), "\n")
cat("  + prior anxiety + psychotropic : ", hr(m1), "\n")
cat(sprintf("\nCovariate HRs: prior anxiety %.2f | prior psychotropic %.2f\n",
    exp(coef(m1)[["prior_anx"]]), exp(coef(m1)[["prior_psych_med"]])))
cat("Interpretation: if the cohesion HR is essentially unchanged, residual confounding by\n",
    "baseline mental health is unlikely to explain the association.\n", sep="")
