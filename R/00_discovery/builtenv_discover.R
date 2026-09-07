# builtenv_discover.R -- print the answer-option distribution for each built-environment / cohesion
# item, so we can see the true response scale (and whether low coverage = different scale vs fewer
# respondents). Informational only.
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

ids <- c(safe=40192384, clean=40192456, upkeep=40192386, sidewalks=40192437, shops=40192436,
         transit=40192440, recreation=40192410, bike=40192431, crime=40192493,
         unsafe_night=40192492, unsafe_day=40192414, graffiti=40192420, vandalism=40192412,
         abandoned=40192469, drug=40192457, alcohol=40192476, trouble=40192404, noisy=40192522,
         hangaround=40192500, trust=40192499, help=40192463)

dd <- run_sql(sprintf("SELECT question_concept_id, answer, COUNT(*) AS n FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)
                       GROUP BY question_concept_id, answer
                       ORDER BY question_concept_id, n DESC", cdr, paste(ids, collapse=",")))
dd$n <- as.numeric(dd$n)
dd$item <- names(ids)[match(as.numeric(dd$question_concept_id), ids)]

for (nm in names(ids)){
  sub <- dd[dd$item==nm, c("answer","n")]
  cat(sprintf("\n### %-12s (%d)  total=%s ###\n", nm, ids[nm], format(sum(sub$n), big.mark=",")))
  print(sub, row.names=FALSE)
}
