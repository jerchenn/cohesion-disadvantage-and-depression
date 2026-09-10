# R/05_recognition/confound_check.R -- diagnostic before any modelling.
# Both surviving signals (discrimination-higher-recognition; low-supply-higher-recognition) run in the
# direction expected if those groups are simply SICKER / more care-ENGAGED. Test that directly:
# baseline utilization (log_util) and ordinal distress severity, by discrimination and by supply tertile.
# Self-contained; run from repo root: run("R/05_recognition/confound_check.R"). Aggregate output only.

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cx <- readRDS("cox_dat2.rds")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))

disc <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%Delayed Or No Care%'") |>
  transmute(person_id, discrim = case_when(
    grepl('always|most of the time|some of the time', ans) ~ 1L,
    grepl('none of the time', ans) ~ 0L, TRUE ~ NA_integer_))

dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all', ans) ~ 0L, grepl('a little', ans) ~ 1L, grepl('moderate', ans) ~ 2L,
    grepl('quite a bit', ans) ~ 3L, grepl('extremely', ans) ~ 4L, TRUE ~ NA_integer_))

d <- cx |> select(person_id, event, sup_t, log_util) |>
  left_join(disc, by = "person_id") |> left_join(dord, by = "person_id") |>
  filter(!is.na(sup_t), !is.na(discrim), !is.na(distress_ord))
cat("analytic n:", nrow(d), "\n")

cat("\n== baseline utilization (log_util) by discrimination ==\n")
print(as.data.frame(d |> group_by(discrim) |> summarise(n = n(),
  util_med = round(median(log_util, na.rm = TRUE), 2), util_mean = round(mean(log_util, na.rm = TRUE), 2))))

cat("\n== baseline utilization by supply tertile ==\n")
print(as.data.frame(d |> group_by(sup_t) |> summarise(n = n(),
  util_med = round(median(log_util, na.rm = TRUE), 2), util_mean = round(mean(log_util, na.rm = TRUE), 2))))

cat("\n== distress severity (ordinal 0-4) by discrimination ==\n")
print(as.data.frame(d |> group_by(discrim) |> summarise(n = n(), distress_mean = round(mean(distress_ord), 2))))

cat("\n== distress severity by supply tertile ==\n")
print(as.data.frame(d |> group_by(sup_t) |> summarise(n = n(), distress_mean = round(mean(distress_ord), 2))))
cat("\nDONE\n")
