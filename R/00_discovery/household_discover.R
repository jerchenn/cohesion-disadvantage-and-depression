# household_discover.R -- find survey items for broader HOUSEHOLD influences to add as predictors:
#   (1) housing insecurity, (2) unpaid/informal care duties, (3) long-term disability.
# Searches ds_survey question text by keyword and prints the matching questions with their
# concept_id, source survey, and respondent count, so we can pick the right items and then build
# scale-aware recodes. Informational only (aggregate).
#
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

groups <- list(
  housing    = c("rent","mortgage","housing","afford","evict","utilit","homeless",
                 "place to live","behind on","move"),
  care       = c("care for","caregiv","looking after","take care","care of","childcare",
                 "dependent","unpaid","help a family"),
  disability = c("disab","serious difficulty","difficulty","deaf","blind","dressing","bathing",
                 "errands","concentrat","remember","walking","climbing","limitation")
)

for (g in names(groups)){
  like <- paste(sprintf("LOWER(question) LIKE '%%%s%%'", groups[[g]]), collapse=" OR ")
  d <- run_sql(sprintf(
    "SELECT survey, question_concept_id, question, COUNT(DISTINCT person_id) AS n
     FROM `%s.ds_survey` WHERE %s
     GROUP BY survey, question_concept_id, question
     ORDER BY n DESC", cdr, like))
  d$n <- as.numeric(d$n)
  d$question <- substr(gsub("\\s+"," ", d$question), 1, 95)
  cat(sprintf("\n========== %s : %d matching questions ==========\n", toupper(g), nrow(d)))
  print(d[, c("n","survey","question_concept_id","question")], row.names=FALSE)
}
cat("\nPick the cleanest item(s) per domain and send them back; then we build recodes and add a\n",
    "'household' predictor family to item_symptom_map.R (and the two-pathway analysis).\n", sep="")
