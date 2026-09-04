# fitbit_scout.R -- feasibility scout for the Fitbit layer (objective, device-measured,
# NOT self-report / NOT utilization-gated => genuinely independent measurement mode).
# Reports which Fitbit tables exist, their coverage, date span, and the analyzable
# intersection with the cohesion at-risk cohort.
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
cox <- readRDS("cox_dat.rds")
cohort_ids <- unique(cox$person_id)

## 1. which Fitbit-related tables exist?
tabs <- run_sql(sprintf("SELECT table_name FROM `%s.INFORMATION_SCHEMA.TABLES`", cdr))
fit <- tabs$table_name[grepl("activity|heart_rate|sleep|steps|fitbit|device", tabs$table_name, ignore.case=TRUE)]
cat("=== Candidate Fitbit tables ===\n"); print(fit)

## 2. per-table coverage + date span + cohort overlap
summ <- function(tbl, datecol){
  ok <- tbl %in% fit; if (!ok){ cat(sprintf("\n[%s] not present\n", tbl)); return(invisible()) }
  info <- tryCatch(run_sql(sprintf(
    "SELECT COUNT(DISTINCT person_id) AS n_people, COUNT(*) AS n_rows, MIN(%s) AS d0, MAX(%s) AS d1 FROM `%s.%s`",
    datecol, datecol, cdr, tbl)), error=function(e){cat("  (query failed:",conditionMessage(e),")\n");NULL})
  if (is.null(info)) return(invisible())
  ids <- tryCatch(run_sql(sprintf("SELECT DISTINCT person_id FROM `%s.%s`", cdr, tbl))$person_id,
                  error=function(e) character(0))
  cat(sprintf("\n[%s]  people=%s  rows=%s  span %s..%s  |  in cohesion cohort: %d\n",
      tbl, info$n_people, info$n_rows, as.character(info$d0), as.character(info$d1),
      length(intersect(ids, cohort_ids))))
}
summ("activity_summary",   "date")
summ("heart_rate_summary", "date")
summ("sleep_daily_summary","sleep_date")
summ("heart_rate_minute_level","datetime")
summ("steps_intraday","datetime")

## 3. peek at column names of the summary tables (so we know the metrics available)
peek <- function(tbl){ if(!(tbl %in% fit)) return(invisible())
  cols <- run_sql(sprintf("SELECT column_name FROM `%s.INFORMATION_SCHEMA.COLUMNS` WHERE table_name='%s'", cdr, tbl))
  cat(sprintf("\n--- columns: %s ---\n", tbl)); print(cols$column_name) }
peek("activity_summary"); peek("heart_rate_summary"); peek("sleep_daily_summary")
