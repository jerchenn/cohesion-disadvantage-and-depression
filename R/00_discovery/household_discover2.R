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
A$n <- as.numeric(A$n); A$question <- substr(gsub("\\s+"," ", A$question), 1, 110)
cat("========== (A) SDOH housing block ==========\n")
print(as.data.frame(A[, c("n","question_concept_id","question")]), row.names=FALSE)

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
cat("\n========== (B) answer options for confirmed items ==========\n")
for (nm in names(ids)){
  sub <- dd[dd$item==nm, c("answer","n")]
  if (nrow(sub)){ cat(sprintf("\n### %-12s (%d) ###\n", nm, ids[nm]))
    print(as.data.frame(sub), row.names=FALSE) }
}
