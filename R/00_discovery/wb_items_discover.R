# wb_items_discover.R -- enumerate the full EHHWB well-being / affect item set
# so the item-level SWB -> depression analysis is built on the real inventory.
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## 1. every question in the EHHWB survey: concept_id, text, n respondents
q_items <- run_sql(sprintf("
SELECT question_concept_id, ANY_VALUE(question) AS question,
       COUNT(DISTINCT person_id) AS n_resp
FROM `%s.ds_survey`
WHERE survey = 'Emotional Health History and Well-Being'
GROUP BY question_concept_id
ORDER BY n_resp DESC", cdr))
cat("=== EHHWB survey items (", nrow(q_items), "questions ) ===\n")
print(as.data.frame(q_items), row.names = FALSE)

## 2. answer distribution for each item (to see the response scale / recode map)
ans <- run_sql(sprintf("
SELECT question_concept_id, answer, COUNT(*) AS n
FROM `%s.ds_survey`
WHERE survey = 'Emotional Health History and Well-Being'
GROUP BY question_concept_id, answer", cdr))
cat("\n=== answer options per item ===\n")
for (qc in q_items$question_concept_id) {
  sub <- ans[ans$question_concept_id == qc, ]
  sub <- sub[order(-as.numeric(sub$n)), ]
  qtext <- substr(q_items$question[q_items$question_concept_id == qc], 1, 70)
  cat("\n[", qc, "] ", qtext, "\n", sep = "")
  print(data.frame(answer = sub$answer, n = sub$n), row.names = FALSE)
}
