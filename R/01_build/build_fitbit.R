# build_fitbit.R -- Fitbit layer: device-measured activity + sleep as an objective,
# non-self-report, non-utilization-gated measurement mode. Person-level means from the
# DAILY SUMMARY tables (never the intraday tables). Merged onto the cohesion at-risk cohort.
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
cox <- readRDS("cox_dat3.rds")

## 1. activity: person-level means over wear-days (steps in plausible range) ----
act <- run_sql(sprintf("
SELECT person_id,
       AVG(steps)                              AS steps,
       AVG(sedentary_minutes)                  AS sed_min,
       AVG(very_active_minutes+fairly_active_minutes) AS mvpa_min,
       COUNT(*)                                AS n_days
FROM `%s.activity_summary`
WHERE steps BETWEEN 100 AND 45000
GROUP BY person_id
HAVING n_days >= 10", cdr))

## 2. sleep: person-level means over main-sleep nights ----
slp <- run_sql(sprintf("
SELECT person_id,
       AVG(minute_asleep)                                  AS asleep_min,
       AVG(minute_asleep / NULLIF(minute_in_bed,0))        AS efficiency,
       AVG(minute_awake + minute_restless)                 AS waso_min,
       COUNT(*)                                            AS n_nights
FROM `%s.sleep_daily_summary`
WHERE LOWER(CAST(is_main_sleep AS STRING))='true'
  AND minute_in_bed BETWEEN 120 AND 900
GROUP BY person_id
HAVING n_nights >= 5", cdr))

for (v in c("steps","sed_min","mvpa_min","n_days")) act[[v]] <- as.numeric(act[[v]])
for (v in c("asleep_min","efficiency","waso_min","n_nights")) slp[[v]] <- as.numeric(slp[[v]])

fb <- merge(cox, act, by="person_id", all.x=TRUE)
fb <- merge(fb, slp, by="person_id", all.x=TRUE)
saveRDS(fb, "fitbit_dat.rds")

cat("Fitbit merge onto cohesion at-risk cohort (n =", nrow(cox), "):\n")
cat(sprintf("  activity: %d people (median %.0f wear-days)\n",
    sum(!is.na(fb$steps)), median(fb$n_days, na.rm=TRUE)))
cat(sprintf("  sleep:    %d people (median %.0f nights)\n",
    sum(!is.na(fb$asleep_min)), median(fb$n_nights, na.rm=TRUE)))
cat("  incident events among activity-havers:", sum(fb$event[!is.na(fb$steps)], na.rm=TRUE), "\n")
cat("\nMetric distributions (5/50/95 pct):\n")
for (v in c("steps","sed_min","mvpa_min","asleep_min","efficiency","waso_min"))
  cat(sprintf("  %-11s %s\n", v, paste(round(quantile(fb[[v]], c(.05,.5,.95), na.rm=TRUE),2), collapse=" / ")))
