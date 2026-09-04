# phenotype_cohort.R -- validated depression phenotype + prospective incident cohort
# needs run_sql + cdr.  Scans condition_occurrence -- allow a few minutes.
cdr <- Sys.getenv("WORKSPACE_CDR")

# 1. Depression phenotype: descendants of key SNOMED depressive-disorder parents, excluding bipolar
dep_set_sql <- sprintf("
SELECT DISTINCT a.descendant_concept_id AS concept_id, c.concept_name
FROM `%s.concept` p
JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id = p.concept_id
JOIN `%s.concept` c ON c.concept_id = a.descendant_concept_id
WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S'
  AND p.concept_name IN ('Depressive disorder','Major depressive disorder',
                         'Recurrent depressive disorder','Dysthymia')
  AND c.standard_concept='S' AND c.domain_id='Condition'
  AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'
", cdr, cdr, cdr)

concepts <- run_sql(dep_set_sql)
cat("Depression phenotype:", nrow(concepts), "concepts. Sample:\n")
print(head(concepts[order(concepts$concept_name), "concept_name", drop=FALSE], 25))
ids <- paste(concepts$concept_id, collapse = ",")

# 2. Prospective cohort with 1-year lookback (needed to call someone depression-free at baseline)
q <- sprintf("
WITH expo AS (
  SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date
  FROM `%s.ds_survey` WHERE question_concept_id = 40192463 GROUP BY person_id),
dep AS (
  SELECT person_id, MIN(condition_start_date) AS first_dep
  FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
ehr AS (
  SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr
  FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT
  COUNT(*)                                                                        AS n_expo_ehr,
  COUNTIF(ehr.first_ehr <= DATE_SUB(expo.expo_date, INTERVAL 365 DAY))            AS n_with_lookback,
  COUNTIF((dep.first_dep IS NULL OR dep.first_dep > expo.expo_date)
          AND ehr.first_ehr <= DATE_SUB(expo.expo_date, INTERVAL 365 DAY))        AS n_at_risk,
  COUNTIF(dep.first_dep > expo.expo_date
          AND ehr.first_ehr <= DATE_SUB(expo.expo_date, INTERVAL 365 DAY))        AS n_incident,
  ROUND(APPROX_QUANTILES(DATE_DIFF(ehr.last_ehr, expo.expo_date, DAY),2)[OFFSET(1)]/365.25,2) AS median_fu_yrs
FROM expo JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)
WHERE ehr.last_ehr >= expo.expo_date
", cdr, cdr, ids, cdr)
print(run_sql(q))
