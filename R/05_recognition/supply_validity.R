# R/05_recognition/supply_validity.R -- is the ZIP3 MH-supply tertile a real instrument or interpolation noise?
# Flat utilization AND flat distress across tertiles (confound_check.R) could mean supply is measuring nothing.
# Test: do supply tertiles differ on things they SHOULD -- area deprivation, rurality, the survey's own
# 'Delayed Medical Care: Rural Area' item -- and on the raw supply values themselves. If indistinguishable,
# report supply as a null instrument, not a finding. Self-contained; run from repo root. Aggregate only.

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
cxs <- cx |> filter(!is.na(sup_t))

## 1. sanity: the supply values themselves should differ sharply across their own tertiles
cat("== raw mh_prov_per_100k by tertile (must be well separated) ==\n")
print(as.data.frame(cxs |> group_by(sup_t) |> summarise(n = n(),
  supply_med = round(median(mh_prov_per_100k), 0),
  supply_p10 = round(quantile(mh_prov_per_100k, .1), 0),
  supply_p90 = round(quantile(mh_prov_per_100k, .9), 0))))

## 2. area deprivation by tertile (higher supply should track LOWER deprivation, if it's real)
if ("deprivation_index" %in% names(cxs)) {
  cat("\n== area deprivation_index by supply tertile ==\n")
  print(as.data.frame(cxs |> group_by(sup_t) |> summarise(n = n(),
    dep_mean = round(mean(deprivation_index, na.rm = TRUE), 3),
    dep_med  = round(median(deprivation_index, na.rm = TRUE), 3))))
} else cat("\n(deprivation_index not in cox_dat2)\n")

## 3. rurality: survey 'Delayed Medical Care: Rural Area' endorsement by supply tertile
rural <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%Rural Area%'") |>
  transmute(person_id, rural_delay = case_when(grepl('yes', ans) ~ 1L, grepl('no', ans) ~ 0L, TRUE ~ NA_integer_))
cxr <- left_join(cxs, rural, by = "person_id")
cat("\n== 'Delayed care: rural area' rate by supply tertile (should be HIGHER where supply is LOW) ==\n")
print(as.data.frame(cxr |> filter(!is.na(rural_delay)) |> group_by(sup_t) |>
  summarise(n = n(), rural_delay_rate = round(mean(rural_delay), 4))))

cat("\nDONE\n")
