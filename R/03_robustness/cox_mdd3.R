# cox_mdd3.R -- referee additional-comment 1: restrict the outcome to major depressive disorder.
# Keeps the corrected primary cohort (cox_dat3.rds, free of ANY depression before baseline, so also free
# of MDD) and redefines the event as incident MDD (descendants of 'Major depressive disorder' and
# 'Recurrent depressive disorder' only); non-MDD depression cases are censored at last EHR record. Reuses
# cox_dat3 covariates (phenotype-independent). Compares with the broad-phenotype primary (HR 0.883).
# Aggregate only. Needs: survival, bigrquery.

library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat3.rds")

mdd_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Major depressive disorder','Recurrent depressive disorder') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)
sv <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
     mdd  AS (SELECT person_id, MIN(condition_start_date) AS first_mdd FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
     ehr  AS (SELECT person_id, MAX(condition_start_date) AS last_ehr FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT expo.person_id, expo.expo_date, mdd.first_mdd, ehr.last_ehr
FROM expo JOIN ehr USING(person_id) LEFT JOIN mdd USING(person_id)", cdr, cdr, mdd_ids, cdr))
sv$expo_date <- as.Date(sv$expo_date); sv$first_mdd <- as.Date(sv$first_mdd); sv$last_ehr <- as.Date(sv$last_ehr)

d <- merge(cox[, setdiff(names(cox), c("event","time_days"))], sv, by="person_id")
d$event <- as.integer(!is.na(d$first_mdd) & d$first_mdd > d$expo_date)
d$time_days <- ifelse(d$event==1, as.numeric(d$first_mdd - d$expo_date), as.numeric(d$last_ehr - d$expo_date))
d <- d[d$time_days > 0, ]

rhs <- "z_cohesion + age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
m <- coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", rhs)), d)
b <- summary(m)$coefficients["z_cohesion",]; ci <- exp(confint(m)["z_cohesion",])
cat(sprintf("MDD-only outcome: HR %.3f (%.3f-%.3f)  n=%d  events=%d\n",
            exp(b[["coef"]]), ci[1], ci[2], m$n, m$nevent))
cat("Compare broad-phenotype primary: HR 0.883 (0.858-0.909).\n")
