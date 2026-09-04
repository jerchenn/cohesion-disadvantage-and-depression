# biomarker_check.R -- how selected are CRP vs HbA1c havers (decides the lead marker)
cdr <- Sys.getenv("WORKSPACE_CDR")

dep_ids <- "SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'"
dep_ids <- sprintf(dep_ids, cdr, cdr, cdr)

q <- sprintf("
WITH cohort AS (SELECT DISTINCT person_id FROM `%s.ds_survey` WHERE question_concept_id=40192463),
util AS (SELECT person_id, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id),
crp AS (SELECT DISTINCT person_id FROM `%s.measurement` WHERE measurement_concept_id IN (SELECT concept_id FROM `%s.concept` WHERE domain_id='Measurement' AND LOWER(concept_name) LIKE '%%c reactive protein%%')),
hb  AS (SELECT DISTINCT person_id FROM `%s.measurement` WHERE measurement_concept_id IN (SELECT concept_id FROM `%s.concept` WHERE domain_id='Measurement' AND LOWER(concept_name) LIKE '%%hemoglobin a1c%%')),
dep AS (SELECT DISTINCT person_id FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s))
SELECT c.person_id,
  DATE_DIFF(CURRENT_DATE, DATE(p.birth_datetime), YEAR) AS age,
  IF(crp.person_id IS NULL,0,1) AS has_crp,
  IF(hb.person_id  IS NULL,0,1) AS has_hba1c,
  IF(dep.person_id IS NULL,0,1) AS dep,
  IFNULL(util.n_cond,0)         AS n_cond
FROM cohort c
JOIN `%s.person` p USING(person_id)
LEFT JOIN util USING(person_id) LEFT JOIN crp USING(person_id)
LEFT JOIN hb USING(person_id)  LEFT JOIN dep USING(person_id)
", cdr, cdr, cdr, cdr, cdr, cdr, cdr, dep_ids, cdr)

bm <- run_sql(q)
for (v in c("age","has_crp","has_hba1c","dep","n_cond")) bm[[v]] <- as.numeric(bm[[v]])

summ <- function(flag){
  s <- split(bm, bm[[flag]])
  data.frame(group = c("without","with"),
             n          = sapply(s, nrow),
             mean_age   = sapply(s, function(x) round(mean(x$age, na.rm=TRUE),1)),
             pct_depr   = sapply(s, function(x) round(100*mean(x$dep),1)),
             mean_conds = sapply(s, function(x) round(mean(x$n_cond),0)))
}
cat("=== CRP: who gets tested (selection) ===\n");   print(summ("has_crp"),   row.names=FALSE)
cat("\n=== HbA1c: who gets tested (selection) ===\n"); print(summ("has_hba1c"), row.names=FALSE)
