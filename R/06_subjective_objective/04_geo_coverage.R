# R/06_subjective_objective/04_geo_coverage.R
# GATING DATA CHECK for the relative-comparison (Easterlin-in-health) design: is ZIP3 geography available
# for the FULL self-rated-health cohort (~670k), or only the 39k SDOH cohort? Geography source =
# ds_zip_code_socioeconomic (per-person zip3 + deprivation), per the Paper 1 build. If coverage is high,
# neighborhood reference groups are usable at scale; if not, we're limited to ~39k or to coarse state.
# Aggregate counts only.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

## columns available in the geography table
cat("== ds_zip_code_socioeconomic columns ==\n")
print(as.data.frame(q("SELECT column_name FROM `__CDR__`.INFORMATION_SCHEMA.COLUMNS
                       WHERE table_name='ds_zip_code_socioeconomic' ORDER BY ordinal_position")))

## total persons with a zip3 in that table
cat("\n== persons in ds_zip_code_socioeconomic ==\n")
print(as.data.frame(q("SELECT COUNT(DISTINCT person_id) n_persons, COUNT(DISTINCT zip3_as_string) n_zip3
                       FROM `__CDR__.ds_zip_code_socioeconomic`")))

## coverage among the self-rated-health cohort (Overall Health survey)
cat("\n== ZIP3 coverage of the self-rated-health cohort ==\n")
print(as.data.frame(q("
  WITH srh AS (SELECT DISTINCT person_id FROM `__CDR__.ds_survey` WHERE question='Overall Health: General Physical Health'),
       geo AS (SELECT DISTINCT person_id, zip3_as_string FROM `__CDR__.ds_zip_code_socioeconomic`)
  SELECT COUNT(DISTINCT s.person_id) srh_n,
         COUNT(DISTINCT g.person_id) srh_with_zip3
  FROM srh s LEFT JOIN geo g USING(person_id)")))

## reference-group feasibility: persons per zip3 among SRH cohort (need enough per zip3 for a stable ref mean)
cat("\n== persons-per-ZIP3 among SRH cohort (reference-group size distribution) ==\n")
print(as.data.frame(q("
  WITH srh AS (SELECT DISTINCT person_id FROM `__CDR__.ds_survey` WHERE question='Overall Health: General Physical Health'),
       g AS (SELECT person_id, zip3_as_string FROM `__CDR__.ds_zip_code_socioeconomic`)
  SELECT APPROX_QUANTILES(cnt, 4) q_persons_per_zip3 FROM (
    SELECT g.zip3_as_string, COUNT(DISTINCT s.person_id) cnt
    FROM srh s JOIN g USING(person_id) GROUP BY g.zip3_as_string)")))

cat("\nDONE -- if coverage high AND persons-per-zip3 healthy (median >~200), neighborhood reference at scale is GO.\n")
