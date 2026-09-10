# R/05_recognition/salvage_model.R
# Decisive test for the supply/rurality recognition reversal: does it survive the DENOMINATOR confound
# (differential pre-baseline EHR observation)? Among the DISTRESSED (distress_ord>=2):
#  (A) Cox time-to-recognition ~ supply, UNADJUSTED vs +prior_ehr_days +log_util +core covariates.
#      If the low-vs-high supply HR attenuates to ~1 after prior_ehr_days, the reversal is denominator artifact.
#  (B) washout-adequacy restriction: require >=3y prior EHR (comparable observation), re-check event rate by tertile.
# Aggregate only. Self-contained; run from repo root.

suppressMessages({library(bigrquery); library(dplyr); library(readr); library(survival)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cx <- readRDS("cox_dat2.rds")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low", "mid", "high"))
cx$sup_t <- relevel(factor(cx$sup_t), ref = "high")   # HR vs high-supply reference

dord <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%distant or cut off%'") |>
  transmute(person_id, distress_ord = case_when(
    grepl('not at all', ans) ~ 0L, grepl('a little', ans) ~ 1L, grepl('moderate', ans) ~ 2L,
    grepl('quite a bit', ans) ~ 3L, grepl('extremely', ans) ~ 4L, TRUE ~ NA_integer_))

d <- cx |> left_join(dord, by = "person_id") |>
  filter(!is.na(sup_t), !is.na(distress_ord), distress_ord >= 2,          # the distressed
         !is.na(time_days), time_days > 0, !is.na(event))
cat("distressed analytic n:", nrow(d), " events:", sum(d$event), "\n")
print(as.data.frame(d |> group_by(sup_t) |> summarise(n = n(), events = sum(event), rate = round(mean(event), 4))))

hr <- function(m, terms) {
  s <- summary(m)$coefficients
  for (t in terms) if (t %in% rownames(s))
    cat(sprintf("  %-14s HR %.2f (%.2f-%.2f)  p=%.3f\n", t, exp(s[t,"coef"]),
        exp(s[t,"coef"]-1.96*s[t,"se(coef)"]), exp(s[t,"coef"]+1.96*s[t,"se(coef)"]), s[t,"Pr(>|z|)"]))
}
tt <- c("sup_tlow","sup_tmid")

cat("\n== (A1) UNADJUSTED: recognition ~ supply ==\n")
hr(coxph(Surv(time_days, event) ~ sup_t, d), tt)

cat("\n== (A2) + prior_ehr_days (the denominator confound) ==\n")
d$prior_yr <- as.numeric(d$prior_ehr_days) / 365.25
hr(coxph(Surv(time_days, event) ~ sup_t + prior_yr, d), tt)

cat("\n== (A3) + prior_ehr_days + log_util + age + sex (full) ==\n")
f <- "sup_t + prior_yr + log_util + age + sex"
hr(coxph(as.formula(paste("Surv(time_days, event) ~", f)), d), tt)

cat("\n== (B) washout-adequacy: restrict to >=3y prior EHR, event rate by tertile ==\n")
db <- d |> filter(as.numeric(prior_ehr_days) >= 3*365.25)
cat("restricted n:", nrow(db), " events:", sum(db$event), "\n")
print(as.data.frame(db |> group_by(sup_t) |> summarise(n = n(), events = sum(event), rate = round(mean(event), 4))))
cat("\n-- supply HR within the >=3y-observation subset --\n")
hr(coxph(Surv(time_days, event) ~ sup_t, db), tt)

cat("\nDONE\n")
