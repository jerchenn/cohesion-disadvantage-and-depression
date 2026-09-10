# R/06_subjective_objective/06_build_analytic.R
# Build the person-level analytic file ONCE (one SQL query -> save to Workbench RDS). Downstream models
# read the RDS (fast). Person-level stays in the Workbench; only aggregate model output leaves (DUCC).
# Objective health = ENROLLMENT MEASUREMENTS ONLY (BMI, BP, HbA1c, glucose) -- universally measured, NOT
# coding-confounded. Subjective = 4 self-rated items. Social = race/income/education. Demographics = age/sex.

suppressMessages({library(bigrquery)})
CDR <- Sys.getenv("WORKSPACE_CDR")
sql <- gsub("__CDR__", CDR, "
WITH
srp AS (SELECT person_id, AVG(CASE WHEN LOWER(answer) LIKE '%excellent%' OR LOWER(answer) LIKE '%excllent%' THEN 5
   WHEN LOWER(answer) LIKE '%very good%' THEN 4 WHEN LOWER(answer) LIKE '%good%' THEN 3
   WHEN LOWER(answer) LIKE '%fair%' THEN 2 WHEN LOWER(answer) LIKE '%poor%' THEN 1 END) sr_phys
   FROM `__CDR__.ds_survey` WHERE question='Overall Health: General Physical Health' GROUP BY person_id),
srm AS (SELECT person_id, AVG(CASE WHEN LOWER(answer) LIKE '%excellent%' OR LOWER(answer) LIKE '%excllent%' THEN 5
   WHEN LOWER(answer) LIKE '%very good%' THEN 4 WHEN LOWER(answer) LIKE '%good%' THEN 3
   WHEN LOWER(answer) LIKE '%fair%' THEN 2 WHEN LOWER(answer) LIKE '%poor%' THEN 1 END) sr_ment
   FROM `__CDR__.ds_survey` WHERE question='Overall Health: General Mental Health' GROUP BY person_id),
srpain AS (SELECT person_id, AVG(SAFE_CAST(REGEXP_EXTRACT(answer,'([0-9]+)') AS INT64)) sr_pain
   FROM `__CDR__.ds_survey` WHERE question='Overall Health: Average Pain 7 Days'
     AND REGEXP_CONTAINS(answer,'[0-9]') GROUP BY person_id),
bmi AS (SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
   WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id),
sbp AS (SELECT person_id, AVG(value_as_number) sbp FROM `__CDR__.measurement`
   WHERE measurement_concept_id=3004249 AND value_as_number BETWEEN 70 AND 250 GROUP BY person_id),
hba1c AS (SELECT person_id, AVG(value_as_number) hba1c FROM `__CDR__.measurement`
   WHERE measurement_concept_id=3004410 AND value_as_number BETWEEN 3 AND 20 GROUP BY person_id),
glu AS (SELECT person_id, AVG(value_as_number) glucose FROM `__CDR__.measurement`
   WHERE measurement_concept_id=3004501 AND value_as_number BETWEEN 30 AND 500 GROUP BY person_id),
race AS (SELECT p.person_id, CASE WHEN c.concept_name='White' THEN 'White' WHEN c.concept_name LIKE 'Black%' THEN 'Black'
   WHEN c.concept_name='Asian' THEN 'Asian' WHEN c.concept_name LIKE 'American Indian%' THEN 'AIAN'
   WHEN c.concept_name LIKE 'More than one%' THEN 'Multi' ELSE 'Other' END race,
   p.year_of_birth, p.gender_concept_id FROM `__CDR__.person` p LEFT JOIN `__CDR__.concept` c ON c.concept_id=p.race_concept_id),
inc AS (SELECT person_id, ANY_VALUE(CASE
   WHEN answer LIKE '%less 10k%' OR answer LIKE '%10k 25k%' OR answer LIKE '%25k 35k%' THEN 1
   WHEN answer LIKE '%35k 50k%' OR answer LIKE '%50k 75k%' THEN 2
   WHEN answer LIKE '%75k 100k%' OR answer LIKE '%100k 150k%' OR answer LIKE '%150k 200k%' OR answer LIKE '%more 200k%' THEN 3 END) income
   FROM `__CDR__.ds_survey` WHERE question='Income: Annual Income' GROUP BY person_id),
edu AS (SELECT person_id, ANY_VALUE(CASE
   WHEN answer LIKE '%Never Attended%' OR answer LIKE '%One Through Four%' OR answer LIKE '%Five Through Eight%' OR answer LIKE '%Nine Through Eleven%' THEN 1
   WHEN answer LIKE '%Twelve Or GED%' THEN 2 WHEN answer LIKE '%College One to Three%' THEN 3
   WHEN answer LIKE '%College Graduate%' OR answer LIKE '%Advanced Degree%' THEN 4 END) educ
   FROM `__CDR__.ds_survey` WHERE question='Education Level: Highest Grade' GROUP BY person_id)
SELECT srp.person_id, srp.sr_phys, srm.sr_ment, srpain.sr_pain,
       bmi.bmi, sbp.sbp, hba1c.hba1c, glu.glucose,
       race.race, race.year_of_birth, race.gender_concept_id, inc.income, edu.educ
FROM srp
LEFT JOIN srm    ON srm.person_id=srp.person_id
LEFT JOIN srpain ON srpain.person_id=srp.person_id
LEFT JOIN bmi    ON bmi.person_id=srp.person_id
LEFT JOIN sbp    ON sbp.person_id=srp.person_id
LEFT JOIN hba1c  ON hba1c.person_id=srp.person_id
LEFT JOIN glu    ON glu.person_id=srp.person_id
LEFT JOIN race   ON race.person_id=srp.person_id
LEFT JOIN inc    ON inc.person_id=srp.person_id
LEFT JOIN edu    ON edu.person_id=srp.person_id
WHERE srp.sr_phys IS NOT NULL")

cat("running assembly query (person-level; downloads once)...\n")
d <- bq_table_download(bq_project_query(Sys.getenv("GOOGLE_PROJECT"), sql), bigint="character")
for (v in c("sr_phys","sr_ment","sr_pain","bmi","sbp","hba1c","glucose","income","educ","year_of_birth"))
  d[[v]] <- as.numeric(d[[v]])
d$age <- 2024 - d$year_of_birth
d$female <- ifelse(d$gender_concept_id=="8532",1L, ifelse(d$gender_concept_id=="8507",0L,NA_integer_))
saveRDS(d, "srh_analytic.rds")
cat("saved srh_analytic.rds  n=", nrow(d), "\n")
cat("coverage:\n")
for (v in c("sr_phys","sr_ment","sr_pain","bmi","sbp","hba1c","glucose","race","income","educ"))
  cat(sprintf("  %-10s %d\n", v, sum(!is.na(d[[v]]))))
cat("\nDONE\n")
