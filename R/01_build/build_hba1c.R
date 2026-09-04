# build_hba1c.R -- Study B (exploratory): does perceived cohesion register in
# cardiometabolic biology (HbA1c)?  CROSS-SECTIONAL, lab-selected, confounded ->
# NOT independent corroboration (lab availability is utilization-driven, same bias
# as EHR ascertainment).  Reported with selection + confounding caveats.
library(tidyverse); library(bigrquery); library(sandwich); library(lmtest)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
cox <- readRDS("cox_dat.rds")     # cohesion + covariates + survival (at-risk cohort)

## 1. discover HbA1c + BMI measurement concepts (coverage)
a1c_c <- run_sql(sprintf("SELECT m.measurement_concept_id AS cid, ANY_VALUE(c.concept_name) AS nm, COUNT(DISTINCT m.person_id) AS n
FROM `%s.measurement` m JOIN `%s.concept` c ON c.concept_id=m.measurement_concept_id
WHERE LOWER(c.concept_name) LIKE '%%hemoglobin a1c%%' OR LOWER(c.concept_name) LIKE '%%hba1c%%'
GROUP BY cid HAVING n > 500 ORDER BY n DESC", cdr, cdr))
cat("=== HbA1c concepts ===\n"); print(as.data.frame(a1c_c), row.names=FALSE)
bmi_c <- run_sql(sprintf("SELECT m.measurement_concept_id AS cid, ANY_VALUE(c.concept_name) AS nm, COUNT(DISTINCT m.person_id) AS n
FROM `%s.measurement` m JOIN `%s.concept` c ON c.concept_id=m.measurement_concept_id
WHERE LOWER(c.concept_name) LIKE '%%body mass index%%'
GROUP BY cid HAVING n > 500 ORDER BY n DESC", cdr, cdr))
cat("\n=== BMI concepts ===\n"); print(as.data.frame(bmi_c), row.names=FALSE)

## 2. pull person-level mean HbA1c (% scale) and mean BMI
a1c_ids <- paste(a1c_c$cid, collapse=",")
a1c <- run_sql(sprintf("SELECT person_id, AVG(value_as_number) AS hba1c, COUNT(*) AS n_a1c
FROM `%s.measurement` WHERE measurement_concept_id IN (%s)
AND value_as_number BETWEEN 3 AND 20 GROUP BY person_id", cdr, a1c_ids))
bmi_ids <- paste(bmi_c$cid, collapse=",")
bmi <- run_sql(sprintf("SELECT person_id, AVG(value_as_number) AS bmi FROM `%s.measurement`
WHERE measurement_concept_id IN (%s) AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id", cdr, bmi_ids))
a1c$hba1c <- as.numeric(a1c$hba1c); bmi$bmi <- as.numeric(bmi$bmi)

d <- merge(cox, a1c[c("person_id","hba1c")], by="person_id", all.x=TRUE)
d <- merge(d, bmi[c("person_id","bmi")], by="person_id", all.x=TRUE)
cat(sprintf("\nHbA1c coverage in cohort: %d / %d (%.1f%%);  BMI: %d\n",
    sum(!is.na(d$hba1c)), nrow(d), 100*mean(!is.na(d$hba1c)), sum(!is.na(d$bmi))))
cat("HbA1c distribution:\n"); print(round(quantile(d$hba1c, c(.05,.25,.5,.75,.95), na.rm=TRUE),2))

## 3. SELECTION check: who has an HbA1c?
d$has_a1c <- !is.na(d$hba1c)
cat("\n=== Selection: HbA1c-havers vs not ===\n")
for (v in c("z_cohesion","age","log_util")) cat(sprintf("  %-11s havers %.3f | non %.3f\n",
    v, mean(d[[v]][d$has_a1c], na.rm=TRUE), mean(d[[v]][!d$has_a1c], na.rm=TRUE)))

## 4. models (HC3): cohesion -> HbA1c, base -> +BMI -> subclinical (<6.5)
covs <- "age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
fitp <- function(m, lab){ ct <- coeftest(m, vcov=vcovHC(m,"HC3"))["z_cohesion",]
  cat(sprintf("  %-22s beta %+.4f  SE %.4f  p %.2g   (N=%d)\n", lab, ct[1], ct[2], ct[4], length(m$residuals))) }
cat("\n=== Cohesion -> HbA1c (percentage points per SD cohesion) ===\n")
fitp(lm(as.formula(paste("hba1c ~ z_cohesion +", covs)), d), "base (demog+SES+util)")
fitp(lm(as.formula(paste("hba1c ~ z_cohesion + bmi +", covs)), d), "+ BMI (block obesity path)")
fitp(lm(as.formula(paste("hba1c ~ z_cohesion +", covs)), d[which(d$hba1c < 6.5),]), "subclinical (<6.5%)")
cat("\nNote: base>+BMI gap = share of any cohesion-HbA1c link running through adiposity.\n")
saveRDS(d[c("person_id","hba1c","bmi")], "hba1c_dat.rds")
