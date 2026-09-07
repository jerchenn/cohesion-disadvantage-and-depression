# timeuse_discover.R -- search ds_survey question text for TIME-USE / WORK-HOURS / time-quantity items
# (hours worked, activity time and frequency, sedentary/sitting, sleep hours, commute, screen, caregiving
# time). Prints matching questions with concept_id, source survey, and respondent count so we can see
# what time-related measurement exists. Informational only (aggregate).
#
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

groups <- list(
  work_hours   = c("hours do you work","hours per week","hours a week","work per week",
                   "how many hours","overtime","shift","employed","employment","occupation","paid work"),
  activity_time= c("physical activit","how often","days per week","times per week","minutes",
                   "moderate","vigorous","exercise","brisk","how much time"),
  sedentary    = c("sitting","sedentary","sit","screen","television","computer"),
  sleep_time   = c("hours of sleep","how many hours","sleep","asleep","hours do you sleep"),
  other_time   = c("time spent","spend time","commute","travel to","leisure","volunteer","caring for","care per")
)

for (g in names(groups)){
  like <- paste(sprintf("LOWER(question) LIKE '%%%s%%'", groups[[g]]), collapse=" OR ")
  d <- run_sql(sprintf(
    "SELECT survey, question_concept_id, question, COUNT(DISTINCT person_id) AS n
     FROM `%s.ds_survey` WHERE %s
     GROUP BY survey, question_concept_id, question
     ORDER BY n DESC", cdr, like))
  d$n <- as.numeric(d$n); d$question <- substr(gsub("\\s+"," ", d$question), 1, 95)
  cat(sprintf("\n========== %s : %d matching questions ==========\n", toupper(g), nrow(d)))
  print(d[, c("n","survey","question_concept_id","question")], row.names=FALSE)
}
cat("\nThese show what time-quantity measurement exists (hours, frequency, minutes). Send back the\n",
    "useful ones and we can add time-use style predictors alongside the household and built-env items.\n", sep="")
