# build_antidep.R -- parallel outcome: incident ANTIDEPRESSANT INITIATION.
# A second, independent ascertainment channel (prescribing, not diagnosis). If cohesion
# predicts lower antidepressant initiation too, it corroborates the EHR-diagnosis result
# against coding/ascertainment idiosyncrasies. Reuses cohesion + covariates from cox_dat.
library(tidyverse); library(bigrquery); library(survival)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
cox <- readRDS("cox_dat.rds")

## antidepressant drug concept set = descendants of ATC 'N06A' (standard Drug concepts)
ad_ids <- sprintf("SELECT a.descendant_concept_id
FROM `%s.concept` p
JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id = p.concept_id
JOIN `%s.concept` c ON c.concept_id = a.descendant_concept_id
WHERE p.vocabulary_id='ATC' AND p.concept_code='N06A'
  AND c.standard_concept='S' AND c.domain_id='Drug'", cdr, cdr, cdr)

## baseline = cohesion-exposure date; require EHR overlap + antidepressant-naive at baseline
surv <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date
              FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
ad AS (SELECT person_id, MIN(drug_exposure_start_date) AS first_ad
       FROM `%s.drug_exposure` WHERE drug_concept_id IN (%s) GROUP BY person_id),
ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr,
               COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT expo.person_id, expo.expo_date, ad.first_ad, ehr.last_ehr, ehr.n_cond
FROM expo JOIN ehr USING(person_id) LEFT JOIN ad USING(person_id)
WHERE ehr.last_ehr >= expo.expo_date
  AND ehr.first_ehr <= DATE_SUB(expo.expo_date, INTERVAL 365 DAY)
  AND (ad.first_ad IS NULL OR ad.first_ad > expo.expo_date)     -- antidepressant-naive at baseline
", cdr, cdr, ad_ids, cdr))
surv$expo_date <- as.Date(surv$expo_date); surv$first_ad <- as.Date(surv$first_ad); surv$last_ehr <- as.Date(surv$last_ehr)
surv$event <- as.integer(!is.na(surv$first_ad))
surv$time_days <- ifelse(surv$event==1, as.numeric(surv$first_ad - surv$expo_date),
                                        as.numeric(surv$last_ehr - surv$expo_date))
surv <- surv[surv$time_days > 0, ]

ad_dat <- merge(cox[, c("person_id","z_cohesion","age","sex","race","ethnicity",
                        "income_f","educ_f","emp_f")],
                surv[, c("person_id","event","time_days","n_cond")], by="person_id")
ad_dat$log_util <- log1p(as.numeric(ad_dat$n_cond))
## collapsed covariates (avoid rare-level separation)
ad_dat$sex_c  <- factor(ifelse(ad_dat$sex %in% c("Female","Male"), ad_dat$sex, "Other"))
ad_dat$race_c <- factor(ifelse(ad_dat$race %in% c("White","Black or African American"), ad_dat$race, "Other"))
ad_dat$ethn_c <- factor(ifelse(ad_dat$ethnicity=="Hispanic or Latino","Hispanic",
                        ifelse(ad_dat$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
covs <- "age + sex_c + race_c + ethn_c + income_f + educ_f + log_util"

cat(sprintf("Antidepressant-initiation cohort: n=%d, events=%d, median FU=%.2f yrs\n",
    nrow(ad_dat), sum(ad_dat$event), median(ad_dat$time_days)/365.25))
m  <- coxph(as.formula(paste("Surv(time_days, event) ~ z_cohesion +", covs)), ad_dat)
m2 <- coxph(as.formula(paste("Surv(time_days, event) ~ z_cohesion +", covs)), ad_dat[ad_dat$time_days>180,])
hr <- function(mm){ ci<-confint(mm); sprintf("%.3f (%.3f-%.3f) p=%.2g",
     exp(coef(mm)[["z_cohesion"]]), exp(ci["z_cohesion",1]), exp(ci["z_cohesion",2]),
     summary(mm)$coefficients["z_cohesion","Pr(>|z|)"]) }
cat("Cohesion -> antidepressant initiation, full   :", hr(m),  "\n")
cat("Cohesion -> antidepressant initiation, lag period:", hr(m2), "\n")
saveRDS(ad_dat, "antidep_dat.rds")
