# R/06_subjective_objective/03a_social_distributions.R
# Fast: just the answer distributions for the social moderators, so mappings can be locked before the
# (SQL-aggregated) B analysis. Small GROUP BY queries only -- seconds, tiny output. Aggregate only.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cat("== income question(s) + answer counts ==\n")
print(as.data.frame(q("SELECT question, answer, COUNT(DISTINCT person_id) n FROM `__CDR__.ds_survey`
   WHERE LOWER(question) LIKE '%income%' GROUP BY question, answer HAVING n>=20 ORDER BY question, n DESC")))

cat("\n== education question(s) + answer counts ==\n")
print(as.data.frame(q("SELECT question, answer, COUNT(DISTINCT person_id) n FROM `__CDR__.ds_survey`
   WHERE LOWER(question) LIKE '%education%' OR LOWER(question) LIKE '%grade%' OR LOWER(question) LIKE '%school%'
   GROUP BY question, answer HAVING n>=20 ORDER BY question, n DESC")))

cat("\n== race concept names + counts ==\n")
print(as.data.frame(q("SELECT c.concept_name race, COUNT(*) n FROM `__CDR__.person` p
   LEFT JOIN `__CDR__.concept` c ON c.concept_id=p.race_concept_id GROUP BY 1 ORDER BY n DESC")))

cat("\n== insurance question + answers ==\n")
print(as.data.frame(q("SELECT question, answer, COUNT(DISTINCT person_id) n FROM `__CDR__.ds_survey`
   WHERE LOWER(question) LIKE '%insurance%health%coverage%' OR question LIKE 'Insurance: Healthcare Coverage'
   GROUP BY question, answer HAVING n>=20 ORDER BY n DESC")))

cat("\nDONE\n")
