# build_cox_swb.R -- SWB-nested prospective cohort (SDOH n EHHWB n EHR)
# baseline = EHHWB date (when wellbeing measured); exposure = cohesion (earlier);
# outcome = incident depression AFTER baseline.  Temporal order: cohesion -> SWB -> depression.
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## 1. items: cohesion + wellbeing + SES + demographics
item_ids <- c(40192463,40192411,40192499,40192417,40192400, 1703980,1704001, 1585375,1585940,1585952)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey` WHERE question_concept_id IN (%s)",
                      cdr, paste(item_ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
labs <- c("40192463"="help","40192411"="getalong","40192499"="trust","40192417"="values","40192400"="watchout",
          "1703980"="happy","1704001"="meaning","1585375"="income","1585940"="education","1585952"="employment")
sv$item <- labs[sv$question_concept_id]
wide <- tidyr::pivot_wider(sv[!is.na(sv$item), c("person_id","item","answer")], names_from=item, values_from=answer)

demo <- run_sql(sprintf(paste("SELECT p.person_id, DATE_DIFF(CURRENT_DATE, DATE(p.birth_datetime), YEAR) AS age,",
  "g.concept_name AS sex, r.concept_name AS race, e.concept_name AS ethnicity FROM `%s.person` p",
  "LEFT JOIN `%s.concept` g ON p.gender_concept_id=g.concept_id",
  "LEFT JOIN `%s.concept` r ON p.race_concept_id=r.concept_id",
  "LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id=e.concept_id"), cdr,cdr,cdr,cdr))
demo$age <- as.numeric(demo$age)
d <- merge(wide, demo, by="person_id")
.a3 <- readRDS("cox_dat3.rds"); d$age <- .a3$age[match(d$person_id, .a3$person_id)]

## 2. score cohesion + wellbeing + SES
map5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
happy_map   <- c("Extremely happy"=6,"Very happy"=5,"Moderately happy"=4,"Moderately unhappy"=3,"Very unhappy"=2,"Extremely unhappy"=1)
meaning_map <- c("Not at all"=1,"A little"=2,"A moderate amount"=3,"Very much"=4,"An extreme amount"=5)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","Skip","Don't know","Prefer not to answer")
rc <- function(x,m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in c("help","getalong","trust","values")) d[[nm]] <- rc(d[[nm]], map5)
d$watchout <- rc(d$watchout, map4)
Z <- scale(as.matrix(d[c("help","getalong","trust","values","watchout")]))
coh <- rowMeans(Z, na.rm=TRUE); coh[rowMeans(!is.na(Z)) < 0.6] <- NA
d$z_cohesion <- scale(coh)[,1]
d$happy_n <- rc(d$happy, happy_map); d$meaning_n <- rc(d$meaning, meaning_map)
d$wb_index <- rowMeans(cbind(scale(d$happy_n)[,1], scale(d$meaning_n)[,1]), na.rm=TRUE)
d$z_wb <- scale(d$wb_index)[,1]
income_map <- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
educ_map <- c("Less than a high school degree or equivalent"=1,"Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"College graduate or advanced degree"=4)
d$income_n <- rc(d$income, income_map); d$educ_n <- rc(d$education, educ_map)
emp <- as.character(d$employment); emp[emp %in% na_vals] <- NA
d$income_f <- addNA(factor(d$income_n)); d$educ_f <- addNA(factor(d$educ_n)); d$emp_f <- addNA(factor(emp))

## 3. baseline = EHHWB date (happiness item 1703980); survival to incident depression after baseline
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)
base <- run_sql(sprintf("
WITH bl AS (SELECT person_id, MIN(DATE(survey_datetime)) AS base_date FROM `%s.ds_survey` WHERE question_concept_id=1703980 GROUP BY person_id),
dep AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT bl.person_id, bl.base_date, dep.first_dep, ehr.last_ehr, ehr.n_cond
FROM bl JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)
WHERE ehr.last_ehr >= bl.base_date AND ehr.first_ehr <= DATE_SUB(bl.base_date, INTERVAL 365 DAY)
  AND (dep.first_dep IS NULL OR dep.first_dep > bl.base_date)
", cdr, cdr, dep_ids, cdr))
base$base_date <- as.Date(base$base_date); base$first_dep <- as.Date(base$first_dep); base$last_ehr <- as.Date(base$last_ehr)
base$n_cond <- as.numeric(base$n_cond)
base$event <- as.integer(!is.na(base$first_dep))
base$time_days <- ifelse(base$event==1, as.numeric(base$first_dep - base$base_date),
                                        as.numeric(base$last_ehr  - base$base_date))
base <- base[base$time_days > 0, ]

swb_dat <- merge(d, base[c("person_id","event","time_days","n_cond")], by="person_id")
swb_dat$log_util <- log1p(swb_dat$n_cond)
saveRDS(swb_dat, "swb_dat.rds")
cat("SWB cohort: n =", nrow(swb_dat), " events =", sum(swb_dat$event),
    " median FU (yrs) =", round(median(swb_dat$time_days)/365.25, 2), "\n")
