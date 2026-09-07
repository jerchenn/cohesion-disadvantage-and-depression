# household_discover2.R -- refine the household search:
#  (A) list the full SDOH housing block (untruncated) to find housing-insecurity items;
#  (B) print answer-option distributions for the confirmed candidate items (disability battery,
#      home ownership, and the two care items) so we can build scale-aware recodes.
# Informational only (aggregate).  INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

## ---- (A) full SDOH housing block (as.data.frame prints ALL rows, not just 10) ----
A <- run_sql(sprintf(
  "SELECT question_concept_id, question, COUNT(DISTINCT person_id) AS n
   FROM `%s.ds_survey`
   WHERE survey = 'Social Determinants of Health'
     AND (LOWER(question) LIKE '%%rent%%' OR LOWER(question) LIKE '%%hous%%'
       OR LOWER(question) LIKE '%%afford%%' OR LOWER(question) LIKE '%%worried%%'
       OR LOWER(question) LIKE '%%steady%%' OR LOWER(question) LIKE '%%place to live%%'
       OR LOWER(question) LIKE '%%utilit%%' OR LOWER(question) LIKE '%%evict%%'
       OR LOWER(question) LIKE '%%mov%%')
   GROUP BY question_concept_id, question ORDER BY n DESC", cdr))
A$n <- as.numeric(A$n); A$question <- gsub("\\s+"," ", A$question)
cat(sprintf("========== (A) SDOH housing block: %d rows (full text) ==========\n", nrow(A)))
for (i in seq_len(nrow(A)))
  cat(sprintf("%-9s | n=%6d | %s\n", A$question_concept_id[i], A$n[i], A$question[i]))

## ---- (B) answer distributions for confirmed candidate items ----
ids <- c(disab_deaf=903573, disab_blind=903574, disab_concen=903575, disab_walk=903576,
         disab_dress=903577, disab_errand=903578,           # Basics 6-item disability battery
         home_own=1585370,                                  # Basics housing tenure
         care_ehhwb=1704022, care_cope=1310136)             # care items (EHHWB, COPE)
dd <- run_sql(sprintf(
  "SELECT question_concept_id, answer, COUNT(*) AS n FROM `%s.ds_survey`
   WHERE question_concept_id IN (%s)
   GROUP BY question_concept_id, answer ORDER BY question_concept_id, n DESC",
  cdr, paste(ids, collapse=",")))
dd$n <- as.numeric(dd$n); dd$item <- names(ids)[match(as.numeric(dd$question_concept_id), ids)]
cat("\n========== (B) answer options for confirmed items (full) ==========\n")
for (nm in names(ids)){
  sub <- dd[dd$item==nm, c("answer","n")]
  if (nrow(sub)){ cat(sprintf("\n### %s (%d) ###\n", nm, ids[nm]))
    for (j in seq_len(nrow(sub))) cat(sprintf("  %6d  %s\n", sub$n[j], sub$answer[j])) }
}
