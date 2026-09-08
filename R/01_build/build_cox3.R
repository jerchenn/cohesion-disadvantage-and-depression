# build_cox3.R -- corrected primary Cox dataset addressing referee concerns 1-3, 8.
#  (3) healthcare utilization is now counted STRICTLY BEFORE baseline (condition_start_date < expo_date),
#      not over all history, so log_util cannot be a post-baseline collider/mediator.
#  (1) adds prior_ehr_days (expo_date - first EHR record) so washout can be stratified by length of
#      observable depression-free history; the depression washout itself is already lifetime-in-EHR.
#  (2) prints a full attrition ladder to reconcile cohort counts (e.g., 99,044 vs 98,795).
#  (8) carries deprivation_index (ZIP3-level) for cluster-robust variance and the interaction analyses.
# Primary design otherwise matches build_cox.R (SDOH cohesion baseline). Saves cox_dat3.rds. Aggregate only.
#
# INPUT: ds_survey, person, concept, concept_ancestor, condition_occurrence, ds_zip_code_socioeconomic.
# Needs: tidyverse, bigrquery.

library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## 1. cohesion items + SES + demographics (whole SDOH cohort) -- identical to build_cox.R
item_ids <- c(40192463,40192411,40192499,40192417,40192400, 1585375,1585940,1585952)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer, survey_datetime FROM `%s.ds_survey` WHERE question_concept_id IN (%s)",
                      cdr, paste(item_ids, collapse=",")))
sv <- sv[order(sv$person_id, sv$question_concept_id, sv$survey_datetime), ]
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
labs <- c("40192463"="help","40192411"="getalong","40192499"="trust","40192417"="values",
          "40192400"="watchout","1585375"="income","1585940"="education","1585952"="employment")
sv$item <- labs[sv$question_concept_id]
wide <- tidyr::pivot_wider(sv[!is.na(sv$item), c("person_id","item","answer")], names_from=item, values_from=answer)
demo <- run_sql(sprintf(paste("SELECT p.person_id, DATE(p.birth_datetime) AS birth_date,",
  "g.concept_name AS sex, r.concept_name AS race, e.concept_name AS ethnicity FROM `%s.person` p",
  "LEFT JOIN `%s.concept` g ON p.gender_concept_id=g.concept_id",
  "LEFT JOIN `%s.concept` r ON p.race_concept_id=r.concept_id",
  "LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id=e.concept_id"), cdr,cdr,cdr,cdr))
demo$birth_date <- as.Date(demo$birth_date)
d <- merge(wide, demo, by="person_id")
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

## 2. survival + PRE-BASELINE utilization + first/last EHR; pulled WITHOUT the cohort filters so we can
##    count attrition stepwise in R.
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)
raw <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
     dep  AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
     ehr  AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr FROM `%s.condition_occurrence` GROUP BY person_id),
     util AS (SELECT co.person_id, COUNT(*) AS n_cond FROM `%s.condition_occurrence` co JOIN expo USING(person_id)
              WHERE co.condition_start_date < expo.expo_date GROUP BY co.person_id)
SELECT expo.person_id, expo.expo_date, dep.first_dep, ehr.first_ehr, ehr.last_ehr, COALESCE(util.n_cond,0) AS n_cond
FROM expo LEFT JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id) LEFT JOIN util USING(person_id)
", cdr, cdr, dep_ids, cdr, cdr))
for (v in c("expo_date","first_dep","first_ehr","last_ehr")) raw[[v]] <- as.Date(raw[[v]])
raw$n_cond <- as.numeric(raw$n_cond)

## attrition ladder (referee concern 2)
cat("=== Attrition ladder ===\n")
cat(sprintf("Answered cohesion item (expo):                 %d\n", nrow(raw)))
s1 <- raw[!is.na(raw$first_ehr), ];                        cat(sprintf("  with any EHR condition record:               %d\n", nrow(s1)))
s2 <- s1[s1$last_ehr >= s1$expo_date, ];                   cat(sprintf("  with EHR record on/after baseline:          %d\n", nrow(s2)))
s3 <- s2[s2$first_ehr <= s2$expo_date - 365, ];            cat(sprintf("  with >=365 d prior EHR observation:         %d\n", nrow(s3)))
s4 <- s3[is.na(s3$first_dep) | s3$first_dep > s3$expo_date, ]; cat(sprintf("  free of depression before baseline:         %d\n", nrow(s4)))
s4$event <- as.integer(!is.na(s4$first_dep))
s4$time_days <- ifelse(s4$event==1, as.numeric(s4$first_dep - s4$expo_date), as.numeric(s4$last_ehr - s4$expo_date))
s5 <- s4[s4$time_days > 0, ];                              cat(sprintf("  follow-up time > 0 (final):                 %d  (events %d)\n", nrow(s5), sum(s5$event)))
s5$prior_ehr_days <- as.numeric(s5$expo_date - s5$first_ehr)

## 3. deprivation (ZIP3-level) for clustering + interaction
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]

cox_dat3 <- merge(d, s5[c("person_id","expo_date","event","time_days","n_cond","prior_ehr_days")], by="person_id")
cox_dat3 <- merge(cox_dat3, geo, by="person_id", all.x=TRUE)
cox_dat3$age <- as.numeric(difftime(cox_dat3$expo_date, cox_dat3$birth_date, units="days"))/365.25
cox_dat3$log_util <- log1p(cox_dat3$n_cond)   # PRE-baseline utilization
saveRDS(cox_dat3, "cox_dat3.rds")
cat(sprintf("\ncox_dat3 saved: n=%d, events=%d, median FU(yr)=%.2f\n",
            nrow(cox_dat3), sum(cox_dat3$event), median(cox_dat3$time_days)/365.25))
cat(sprintf("pre-baseline n_cond: median=%.0f (IQR %.0f-%.0f); with zero prior conditions: %.1f%%\n",
            median(cox_dat3$n_cond), quantile(cox_dat3$n_cond,.25), quantile(cox_dat3$n_cond,.75),
            100*mean(cox_dat3$n_cond==0)))
cat(sprintf("prior_ehr_days: median=%.0f (%.1f yr); with >=730 d: %.1f%%; >=1095 d: %.1f%%\n",
            median(cox_dat3$prior_ehr_days), median(cox_dat3$prior_ehr_days)/365.25,
            100*mean(cox_dat3$prior_ehr_days>=730), 100*mean(cox_dat3$prior_ehr_days>=1095)))
