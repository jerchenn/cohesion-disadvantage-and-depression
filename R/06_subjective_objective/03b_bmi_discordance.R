# R/06_subjective_objective/03b_bmi_discordance.R
# Core B result (FAST, SQL-aggregated): at the SAME objective BMI, does self-rated PHYSICAL health differ
# across social groups (race, income, education)? Systematic shifts = socially-patterned reporting
# heterogeneity on a clean, coding-free anchor. All group-by in BigQuery; returns only small tables.
# BMI = enrollment measurement (universally measured, NOT coding-confounded). Aggregate only (HAVING n>=20).

suppressMessages({library(bigrquery)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

tmpl <- "
WITH
sr AS (SELECT person_id, AVG(CASE
         WHEN LOWER(answer) LIKE '%excellent%' OR LOWER(answer) LIKE '%excllent%' THEN 5
         WHEN LOWER(answer) LIKE '%very good%' THEN 4
         WHEN LOWER(answer) LIKE '%good%'      THEN 3
         WHEN LOWER(answer) LIKE '%fair%'      THEN 2
         WHEN LOWER(answer) LIKE '%poor%'      THEN 1 END) sr_phys
       FROM `__CDR__.ds_survey` WHERE question='Overall Health: General Physical Health' GROUP BY person_id),
bmi AS (SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
        WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id),
__MODCTE__
SELECT CASE WHEN b.bmi<25 THEN '1_<25' WHEN b.bmi<30 THEN '2_25-30' WHEN b.bmi<35 THEN '3_30-35' ELSE '4_35+' END bmi_cat,
       m.grp AS grp, COUNT(*) n, ROUND(AVG(s.sr_phys),3) sr_phys_mean
FROM sr s
JOIN bmi b ON b.person_id = s.person_id
JOIN m    ON m.person_id = s.person_id
WHERE s.sr_phys IS NOT NULL AND m.grp IS NOT NULL
GROUP BY bmi_cat, m.grp HAVING COUNT(*)>=20 ORDER BY m.grp, bmi_cat"

runmod <- function(modcte, name) {
  sql <- gsub("__MODCTE__", modcte, tmpl, fixed=TRUE)
  out <- tryCatch(q(sql), error=function(e){ cat("SQL ERROR:\n", conditionMessage(e), "\n"); NULL })
  if (!is.null(out)) { out$moderator <- name; as.data.frame(out) } else invisible(NULL)
}

race_cte <- "race_m AS (SELECT p.person_id, CASE
    WHEN c.concept_name='White' THEN 'White' WHEN c.concept_name LIKE 'Black%' THEN 'Black'
    WHEN c.concept_name='Asian' THEN 'Asian' WHEN c.concept_name LIKE 'American Indian%' THEN 'AIAN'
    WHEN c.concept_name LIKE 'More than one%' THEN 'Multi' ELSE NULL END grp
  FROM `__CDR__.person` p LEFT JOIN `__CDR__.concept` c ON c.concept_id=p.race_concept_id)"
inc_cte <- "inc_m AS (SELECT person_id, ANY_VALUE(CASE
    WHEN answer LIKE '%less 10k%' OR answer LIKE '%10k 25k%' OR answer LIKE '%25k 35k%' THEN '1_<35k'
    WHEN answer LIKE '%35k 50k%' OR answer LIKE '%50k 75k%' THEN '2_35-75k'
    WHEN answer LIKE '%75k 100k%' OR answer LIKE '%100k 150k%' OR answer LIKE '%150k 200k%' OR answer LIKE '%more 200k%' THEN '3_75k+'
    END) grp FROM `__CDR__.ds_survey` WHERE question='Income: Annual Income' GROUP BY person_id)"
edu_cte <- "edu_m AS (SELECT person_id, ANY_VALUE(CASE
    WHEN answer LIKE '%Never Attended%' OR answer LIKE '%One Through Four%' OR answer LIKE '%Five Through Eight%' OR answer LIKE '%Nine Through Eleven%' THEN '1_<HS'
    WHEN answer LIKE '%Twelve Or GED%' THEN '2_HS'
    WHEN answer LIKE '%College One to Three%' THEN '3_SomeCol'
    WHEN answer LIKE '%College Graduate%' OR answer LIKE '%Advanced Degree%' THEN '4_College+'
    END) grp FROM `__CDR__.ds_survey` WHERE question='Education Level: Highest Grade' GROUP BY person_id)"

# rename the moderator CTE alias to 'm' so the template's `JOIN m` works
fix <- function(cte, alias) sub(paste0(alias,"_m AS"), "m AS", cte, fixed=TRUE)

cat("== B: mean self-rated PHYSICAL health at fixed BMI, by RACE ==\n")
print(runmod(fix(race_cte,"race"), "race"))
cat("\n== ... by INCOME ==\n")
print(runmod(fix(inc_cte,"inc"), "income"))
cat("\n== ... by EDUCATION ==\n")
print(runmod(fix(edu_cte,"edu"), "educ"))

cat("\nREAD: within a BMI band, if sr_phys_mean differs across groups -> socially-patterned reporting\n",
    "(same objective body mass, different self-rated health). That is the B thesis.\n")
cat("\nDONE\n")
