# household_discover3.R -- full question text + answer options for the remaining household candidates,
# so we can finalize recodes: residential moves, housing type, the 2-item food-insecurity screener,
# and the ambiguous EHHWB "care" item (need its wording). Full output via cat(). Aggregate only.
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

ids <- c(moves=40192441, housing_type=40192458, food1=40192517, food2=40192426, care_ehhwb=1704022)

## full question text
qt <- run_sql(sprintf("SELECT DISTINCT question_concept_id, survey, question FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(ids, collapse=",")))
qt$question <- gsub("\\s+"," ", qt$question)
cat("========== full question text ==========\n")
for (i in seq_len(nrow(qt)))
  cat(sprintf("%s (%s) | %s\n", qt$question_concept_id[i], qt$survey[i], qt$question[i]))

## answer options per item
dd <- run_sql(sprintf("SELECT question_concept_id, answer, COUNT(*) AS n FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)
                       GROUP BY question_concept_id, answer ORDER BY question_concept_id, n DESC",
                      cdr, paste(ids, collapse=",")))
dd$n <- as.numeric(dd$n); dd$item <- names(ids)[match(as.numeric(dd$question_concept_id), ids)]
cat("\n========== answer options ==========\n")
for (nm in names(ids)){
  sub <- dd[dd$item==nm, c("answer","n")]
  if (nrow(sub)){ cat(sprintf("\n### %s (%d) ###\n", nm, ids[nm]))
    for (j in seq_len(nrow(sub))) cat(sprintf("  %7d  %s\n", sub$n[j], sub$answer[j])) }
}
