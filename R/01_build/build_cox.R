# build_cox.R -- assemble the prospective Cox analysis dataset (primary analysis)
# exposure = perceived cohesion; outcome = time to incident depression; cohort = SDOH n EHR, at-risk
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## 1. cohesion items + SES + demographics (whole SDOH cohort, not just EHHWB completers)
item_ids <- c(40192463,40192411,40192499,40192417,40192400, 1585375,1585940,1585952)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey` WHERE question_concept_id IN (%s)",
                      cdr, paste(item_ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
labs <- c("40192463"="help","40192411"="getalong","40192499"="trust","40192417"="values",
          "40192400"="watchout","1585375"="income","1585940"="education","1585952"="employment")
sv$item <- labs[sv$question_concept_id]
wide <- tidyr::pivot_wider(sv[!is.na(sv$item), c("person_id","item","answer")], names_from=item, values_from=answer)

demo <- run_sql(sprintf(paste("SELECT p.person_id, DATE_DIFF(CURRENT_DATE, DATE(p.birth_datetime), YEAR) AS age,",
  "g.concept_name AS sex, r.concept_name AS race, e.concept_name AS ethnicity FROM `%s.person` p",
  "LEFT JOIN `%s.concept` g ON p.gender_concept_id=g.concept_id",
  "LEFT JOIN `%s.concept` r ON p.race_concept_id=r.concept_id",
  "LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id=e.concept_id"), cdr,cdr,cdr,cdr))
demo$age <- as.numeric(demo$age)
d <- merge(wide, demo, by="person_id")

## 2. score cohesion (watchout 4-pt, others 5-pt: z each item, average, re-z) + SES
map5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","Skip","Don't know","Prefer not to answer")
rc <- function(x,m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in c("help","getalong","trust","values")) d[[nm]] <- rc(d[[nm]], map5)
d$watchout <- rc(d$watchout, map4)
Z <- scale(as.matrix(d[c("help","getalong","trust","values","watchout")]))
coh <- rowMeans(Z, na.rm=TRUE); coh[rowMeans(!is.na(Z)) < 0.6] <- NA
d$z_cohesion <- scale(coh)[,1]
income_map <- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
educ_map <- c("Less than a high school degree or equivalent"=1,"Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"College graduate or advanced degree"=4)
d$income_n <- rc(d$income, income_map); d$educ_n <- rc(d$education, educ_map)
emp <- as.character(d$employment); emp[emp %in% na_vals] <- NA
d$income_f <- addNA(factor(d$income_n)); d$educ_f <- addNA(factor(d$educ_n)); d$emp_f <- addNA(factor(emp))

## 3. survival outcome (time to incident depression) + healthcare utilization
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)
surv <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
dep AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT expo.person_id, expo.expo_date, dep.first_dep, ehr.last_ehr, ehr.n_cond
FROM expo JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)
WHERE ehr.last_ehr >= expo.expo_date AND ehr.first_ehr <= DATE_SUB(expo.expo_date, INTERVAL 365 DAY)
  AND (dep.first_dep IS NULL OR dep.first_dep > expo.expo_date)
", cdr, cdr, dep_ids, cdr))
surv$expo_date <- as.Date(surv$expo_date); surv$first_dep <- as.Date(surv$first_dep); surv$last_ehr <- as.Date(surv$last_ehr)
surv$n_cond <- as.numeric(surv$n_cond)
surv$event <- as.integer(!is.na(surv$first_dep))
surv$time_days <- ifelse(surv$event==1, as.numeric(surv$first_dep - surv$expo_date),
                                         as.numeric(surv$last_ehr  - surv$expo_date))
surv <- surv[surv$time_days > 0, ]

cox_dat <- merge(d, surv[c("person_id","event","time_days","n_cond")], by="person_id")
cox_dat$log_util <- log1p(cox_dat$n_cond)
saveRDS(cox_dat, "cox_dat.rds")
cat("Cox dataset: n =", nrow(cox_dat), " events =", sum(cox_dat$event),
    " median FU (yrs) =", round(median(cox_dat$time_days)/365.25, 2), "\n")
