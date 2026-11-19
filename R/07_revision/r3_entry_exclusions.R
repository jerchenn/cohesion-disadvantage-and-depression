# r3_entry_exclusions.R -- Reviewer 2, point 3: who is excluded by the requirement of an EHR record on or after
# baseline? Compares participants who meet every other entry criterion (linked EHR, >=365 days of prior EHR,
# no recorded depression before baseline, computable cohesion score) but have no record on or after baseline
# with the at-risk cohort. Uses the queries of build_cox3.R. Aggregate output only; counts under 20 masked.
source("R/07_revision/_common.R", local = TRUE); bq()
mask2 <- function(x) ifelse(x > 0 & x < 20, "<20", format(x, big.mark = ","))
raw <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
     ehr  AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr FROM `%s.condition_occurrence` GROUP BY person_id),
     dep  AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (
               SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id
               JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S'
               AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia')
               AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%') GROUP BY person_id)
SELECT expo.person_id, expo.expo_date, ehr.first_ehr, ehr.last_ehr, dep.first_dep FROM expo JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)",
  cdr, cdr, cdr, cdr, cdr, cdr))
for (v in c("expo_date", "first_ehr", "last_ehr", "first_dep")) raw[[v]] <- as.Date(raw[[v]])
raw$person_id <- as.character(raw$person_id)
elig <- raw[raw$first_ehr <= raw$expo_date - 365 & (is.na(raw$first_dep) | raw$first_dep > raw$expo_date), ]
excl <- elig[elig$last_ehr < elig$expo_date, ]
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$person_id <- as.character(geo$person_id); geo <- geo[!duplicated(geo$person_id), ]
## covariates for the excluded group come from the same survey and person tables, recoded as in build_cox3.R
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer, survey_datetime FROM `%s.ds_survey` WHERE question_concept_id IN (40192463,40192411,40192499,40192417,40192400,1585375)", cdr))
sv$person_id <- as.character(sv$person_id); sv <- sv[sv$person_id %in% excl$person_id, ]
sv <- sv[order(sv$person_id, sv$question_concept_id, sv$survey_datetime, sv$answer), ]; sv <- sv[!duplicated(sv[c("person_id", "question_concept_id")]), ]
w <- tidyr::pivot_wider(sv[c("person_id", "question_concept_id", "answer")], names_from = question_concept_id, values_from = answer)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","PMI: Dont Know","Skip","Don't know","Prefer not to answer","Response removed due to invalid value")
num <- function(x, m) { x <- as.character(x); x[x %in% na_vals] <- NA; as.numeric(m[x]) }
m5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5); m4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
inc_map <- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
items <- cbind(num(w[["40192463"]], m5), num(w[["40192411"]], m5), num(w[["40192499"]], m5), num(w[["40192417"]], m5), num(w[["40192400"]], m4))
## cohesion on the at-risk cohort's scale: raw item mean, compared with the same raw mean in the cohort
cx <- cox; raw_coh <- function(M) ifelse(rowMeans(!is.na(M)) >= 0.6, rowMeans(M, na.rm = TRUE), NA)
excl_coh <- raw_coh(items[, 1:4]); coh_cohort <- raw_coh(as.matrix(cx[c("help", "getalong", "trust", "values")]))
cat(sprintf("otherwise eligible: %s | excluded for no EHR record on or after baseline: %s (%.1f%%)\n", mask2(nrow(elig)), mask2(nrow(excl)), 100 * nrow(excl) / nrow(elig)))
cat(sprintf("at-risk cohort (included): %s\n", mask2(nrow(cx))))
ex_inc <- num(w[["1585375"]], inc_map)
ex_dep <- geo$deprivation_index[match(excl$person_id, geo$person_id)]; ex_dep <- as.numeric(ex_dep)
cat(sprintf("\n%-46s %18s %18s\n", "", "excluded", "at-risk cohort"))
cat(sprintf("%-46s %18s %18s\n", "cohesion answered (>=3 of 4 agree-scale items), %", sprintf("%.1f", 100 * mean(!is.na(excl_coh))), sprintf("%.1f", 100 * mean(!is.na(coh_cohort)))))
cat(sprintf("%-46s %18s %18s\n", "cohesion, mean of 4 agree-scale items (1-5)", sprintf("%.2f", mean(excl_coh, na.rm = TRUE)), sprintf("%.2f", mean(coh_cohort, na.rm = TRUE))))
cat(sprintf("%-46s %18s %18s\n", "income <$35k, % of those reporting", sprintf("%.1f", 100 * mean(ex_inc <= 3, na.rm = TRUE)), sprintf("%.1f", 100 * mean(cx$income_n <= 3, na.rm = TRUE))))
cat(sprintf("%-46s %18s %18s\n", "income >=$75k, % of those reporting", sprintf("%.1f", 100 * mean(ex_inc >= 6, na.rm = TRUE)), sprintf("%.1f", 100 * mean(cx$income_n >= 6, na.rm = TRUE))))
cat(sprintf("%-46s %18s %18s\n", "income not reported, %", sprintf("%.1f", 100 * mean(is.na(ex_inc))), sprintf("%.1f", 100 * mean(is.na(cx$income_n)))))
cat(sprintf("%-46s %18s %18s\n", "area deprivation index, mean", sprintf("%.3f", mean(ex_dep, na.rm = TRUE)), sprintf("%.3f", mean(cx$deprivation_index, na.rm = TRUE))))
cat(sprintf("%-46s %18s %18s\n", "years from last EHR record to baseline, median", sprintf("%.1f", median(as.numeric(excl$expo_date - excl$last_ehr)) / 365.25), "-"))
cat(sprintf("%-46s %18s %18s\n", "baseline year 2023-2024, %", sprintf("%.1f", 100 * mean(format(excl$expo_date, "%Y") >= "2023")), sprintf("%.1f", 100 * mean(format(cx$expo_date, "%Y") >= "2023"))))
