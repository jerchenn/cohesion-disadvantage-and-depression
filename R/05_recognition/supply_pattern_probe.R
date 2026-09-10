# R/05_recognition/supply_pattern_probe.R
# The supply/rurality recognition reversal (higher recognition where supply is LOW) -- is it real, or
# (a) the distressed subgroup being sicker in low-supply areas, or (b) an EHR ASCERTAINMENT artifact
# (rural records more complete/concentrated)? Two checks before any substantive model. Aggregate only.
#
# NOTE construct: at ZIP3, "low MH supply" == "more rural" (rural-delay 3.4% vs 1.4% across tertiles).
# Treat this arm as RURALITY / provider-scarcity, not pure MH-provider availability.

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

dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all', ans) ~ 0L, grepl('a little', ans) ~ 1L, grepl('moderate', ans) ~ 2L,
    grepl('quite a bit', ans) ~ 3L, grepl('extremely', ans) ~ 4L, TRUE ~ NA_integer_))
d <- cx |> left_join(dord, by = "person_id") |> filter(!is.na(sup_t), !is.na(distress_ord))
d$distress_bin <- as.integer(d$distress_ord >= 2)

## CHECK A -- WITHIN the distressed (where the recognition contrast lives): are low-supply distressed
## people simply sicker / more engaged? (advisor: flat marginals don't rule this out)
cat("== within distress_bin==1: severity + utilization by supply tertile ==\n")
print(as.data.frame(d |> filter(distress_bin == 1) |> group_by(sup_t) |> summarise(
  n = n(), distress_ord_mean = round(mean(distress_ord), 2),
  log_util_med = round(median(log_util, na.rm = TRUE), 2))))

## CHECK B -- ASCERTAINMENT: is EHR capture just denser/more concentrated in low-supply(rural) areas?
## Total condition-record count per person by tertile (proxy for EHR completeness/capture).
ehr <- q("SELECT person_id, COUNT(*) ehr_cond FROM `__CDR__.condition_occurrence` GROUP BY person_id")
ehr$ehr_cond <- as.numeric(ehr$ehr_cond)
de <- d |> select(person_id, sup_t) |> left_join(ehr, by = "person_id")   # keep only needed cols -> no name clash
de$ehr_cond[is.na(de$ehr_cond)] <- 0
cat("\n== EHR record density (total condition rows) by supply tertile ==\n")
print(as.data.frame(de |> group_by(sup_t) |> summarise(
  n = n(), cond_med = round(median(ehr_cond), 0), cond_mean = round(mean(ehr_cond), 0))))

## site concentration: distinct EHR sites per tertile (fewer sites = more concentrated capture)
site <- tryCatch(
  q("SELECT COUNT(DISTINCT src_id) n_site FROM `__CDR__.condition_occurrence`"), error = function(e) NULL)
if (!is.null(site)) { cat("\n== distinct EHR source sites (whole cohort) ==\n"); print(as.data.frame(site)) } else
  cat("\n(site/src_id column not available in condition_occurrence)\n")

cat("\nDONE\n")
