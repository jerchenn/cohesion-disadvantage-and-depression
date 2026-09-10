# R/05_recognition/denominator_check.R
# The last threat to the "higher recognition in low-supply/rural areas" pattern: a DENOMINATOR effect.
# Outcome = first depression code AFTER baseline, in a cohort with pre-baseline depression washed out.
# Low-supply areas have fewer records + shorter observation -> fewer depression dx captured pre-baseline
# -> fewer washed out -> more remain as "incident". If the pre-baseline depression EXCLUSION rate is lower
# and prior-EHR observation is SHORTER in low-supply tertiles, the reversal is mechanical, not recognition.
# Aggregate only. Self-contained; run from repo root.

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
cx <- cx |> filter(!is.na(sup_t))

## 1. prior-EHR observation length by tertile (shorter in low-supply => less chance to capture/washout)
cat("== prior-EHR observation by supply tertile ==\n")
pe_col <- intersect(c("prior_ehr_days", "prior_obs_days", "prior_ehr_years"), names(cx))
if (length(pe_col)) {
  cx$prior <- as.numeric(cx[[pe_col[1]]])
  print(as.data.frame(cx |> group_by(sup_t) |> summarise(n = n(),
    prior_med = round(median(prior, na.rm = TRUE), 1), prior_mean = round(mean(prior, na.rm = TRUE), 1))))
  cat("(column used:", pe_col[1], ")\n")
} else cat("(no prior-EHR column in cox_dat2 -- names:", paste(names(cx), collapse=", "), ")\n")

## 2. full ordinal distress distribution by tertile (advisor: not the truncated within-distressed mean)
dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all', ans) ~ 0L, grepl('a little', ans) ~ 1L, grepl('moderate', ans) ~ 2L,
    grepl('quite a bit', ans) ~ 3L, grepl('extremely', ans) ~ 4L, TRUE ~ NA_integer_))
dd <- cx |> left_join(dord, by = "person_id") |> filter(!is.na(distress_ord))
cat("\n== full distress_ord distribution by supply tertile (row %) ==\n")
tab <- table(dd$sup_t, dd$distress_ord)
print(round(100 * prop.table(tab, 1), 1))

## 3. THE denominator test: pre-baseline depression EXCLUSION rate by tertile.
## Build the at-risk denominator (SDOH-anchored, cox_dat3) BEFORE the depression washout, flag who had a
## pre-baseline depression code, and compare exclusion rate across supply tertiles.
d3 <- readRDS("cox_dat3.rds")
d3$zip3 <- substr(as.character(d3$zip3_as_string), 1, 3)
d3 <- left_join(d3, sup, by = "zip3")
d3$sup_t <- cut(d3$mh_prov_per_100k, quantile(d3$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))
cat("\n== proxy denominator check: incident-event RATE and prior-obs among the ANALYSIS cohort by tertile ==\n")
print(as.data.frame(d3 |> filter(!is.na(sup_t)) |> group_by(sup_t) |> summarise(
  n = n(), event_rate = round(mean(event), 4),
  prior_med = if ("prior_ehr_days" %in% names(d3)) round(median(as.numeric(prior_ehr_days), na.rm=TRUE),0) else NA)))
cat("\nNOTE: if a pre-washout roster with a pre-baseline-depression flag exists in the build, report the\n",
    "EXCLUSION rate by tertile directly. This proxy shows event-rate + prior-obs by tertile as a first look.\n")
cat("\nDONE\n")
