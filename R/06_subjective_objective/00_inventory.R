# R/06_subjective_objective/00_inventory.R
# SSM-style study: systematic SUBJECTIVE-vs-OBJECTIVE concordance, multilevel moderation, individual-level.
# STEP 0 = DISCOVERY/INVENTORY (pre-specify the analysis matrix from what actually exists; no association
# mining yet). Catalog: (A) subjective measures + n, (B) matchable objective anchors + n, (C) multilevel
# covariates, (D) cohort overlaps. Aggregate counts only (DUCC). Nothing modelled here.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

## ---- (A) SUBJECTIVE measures: person-count for key constructs across surveys ----
cat("== (A) SUBJECTIVE measure availability (persons answering) ==\n")
subj <- q("
  SELECT survey, question, COUNT(DISTINCT person_id) n
  FROM `__CDR__.ds_survey`
  WHERE LOWER(question) LIKE '%how happy%'                 -- happiness
     OR LOWER(question) LIKE '%general%health%'            -- self-rated health
     OR LOWER(question) LIKE '%would you say%health%'
     OR LOWER(question) LIKE '%little interest or pleasure%'-- PHQ
     OR LOWER(question) LIKE '%feeling down, depressed%'
     OR LOWER(question) LIKE '%feeling nervous%'           -- GAD
     OR LOWER(question) LIKE '%worrying too much%'
     OR LOWER(question) LIKE '%stress%'                    -- perceived stress
     OR LOWER(question) LIKE '%lonely%' OR LOWER(question) LIKE '%left out%' OR LOWER(question) LIKE '%isolated%'
     OR LOWER(question) LIKE '%neighbor%'                  -- neighborhood perception
     OR LOWER(question) LIKE '%safe%'
     OR LOWER(question) LIKE '%discriminat%' OR LOWER(question) LIKE '%treated with less%'
     OR LOWER(question) LIKE '%meaning%' OR LOWER(question) LIKE '%life satisf%'
     OR LOWER(question) LIKE '%pain%'                      -- pain perception
  GROUP BY survey, question HAVING n >= 20 ORDER BY n DESC")
print(as.data.frame(subj), max = 400)

## ---- (B) OBJECTIVE anchors: most-populated physical measurements & labs ----
cat("\n== (B) OBJECTIVE measurements/labs (top by person-count) ==\n")
obj <- q("
  SELECT m.measurement_concept_id, c.concept_name, COUNT(DISTINCT m.person_id) n_persons
  FROM `__CDR__.measurement` m JOIN `__CDR__.concept` c ON c.concept_id = m.measurement_concept_id
  GROUP BY 1,2 HAVING n_persons >= 1000 ORDER BY n_persons DESC LIMIT 30")
print(as.data.frame(obj))

## objective EHR condition domains (phenotype anchors) -- broad prevalence via a few seed ancestors
cat("\n== (B2) objective EHR condition prevalence (seed ancestors: depression/anxiety/diabetes/htn/obesity/pain) ==\n")
seeds <- c(depression=440383, anxiety=442077, diabetes=201820, hypertension=316866, obesity=433736, pain=4329041)
for (nm in names(seeds)) {
  r <- q(sprintf("SELECT COUNT(DISTINCT co.person_id) n FROM `__CDR__.condition_occurrence` co
    JOIN `__CDR__.concept_ancestor` ca ON ca.descendant_concept_id = co.condition_concept_id
    WHERE ca.ancestor_concept_id = %d", seeds[[nm]]))
  cat(sprintf("  %-12s (%d): %s persons\n", nm, seeds[[nm]], r$n[1]))
}

## ---- (C) MULTILEVEL covariates already in hand ----
cat("\n== (C) covariates in cox_dat2 (individual + area) ==\n")
cx <- readRDS("cox_dat2.rds")
cat("  individual/area cols:", paste(intersect(c("age","sex","race","ethnicity","income_n","educ_n","emp_f",
    "deprivation_index","zip3_as_string","log_util"), names(cx)), collapse=", "), "\n")
cat("  + external ZIP3 linkages available: MH supply (zip3_mh_supply.csv), rurality (survey item)\n")

## ---- (D) cohort overlaps ----
cat("\n== (D) cohort sizes ==\n")
cat("  cox_dat2 (SDOH∩EHHWB∩EHR, re-anchored):", nrow(cx), "\n")
cat("  wearables (activity_summary):", q("SELECT COUNT(DISTINCT person_id) n FROM `__CDR__.activity_summary`")$n[1], "\n")

cat("\nDONE -- use this to PRE-SPECIFY the subjective x objective matrix + moderators before any association analysis.\n")
