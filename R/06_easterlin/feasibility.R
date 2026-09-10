# R/06_easterlin/feasibility.R
# Iterated thesis: EASTERLIN DECOUPLING. Counter-consensus: subjective wellbeing does NOT track objective
# circumstance the way consensus assumes; RELATIVE position governs SWB beyond ABSOLUTE health/income.
# Objective inputs (EHR morbidity, BMI) vs SUBJECTIVE outcome (happiness) -> no common-method bias.
#
# THIS SCRIPT = feasibility + well-posedness kill tests (cheap, run first):
#  KA: does SWB vary with objective health at all? (if flat -> decoupling trivial/ill-posed)
#  KB: does relative position (own income WITHIN area-deprivation strata) associate with SWB beyond absolute?
# Uses cox_dat2 (income_n, deprivation_index, log_util) + happiness item. Aggregate only. Dedupe per person.

suppressMessages({library(bigrquery); library(dplyr); library(readr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cx <- readRDS("cox_dat2.rds")
cat("base cohort n:", nrow(cx), " | has:", paste(intersect(c("income_n","deprivation_index","log_util","age","sex"), names(cx)), collapse=", "), "\n")

## --- SWB: 'In general, how happy are you?' (dedupe one row/person; map to numeric, higher = happier) ---
h <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%how happy are you%'")
cat("\n== happiness answer distribution (to lock the mapping) ==\n"); print(sort(table(h$ans), decreasing=TRUE))
h <- h |> mutate(happy = case_when(
    grepl('very happy|extremely happy', ans) ~ 3L,
    grepl('pretty happy|quite happy|fairly happy', ans) ~ 2L,
    grepl('not too happy|not very happy', ans) ~ 1L,
    grepl('not at all happy', ans) ~ 0L, TRUE ~ NA_integer_)) |>
  filter(!is.na(happy)) |> group_by(person_id) |> summarise(happy = mean(happy), .groups="drop")   # one row/person
cat("mapping check -- dedup rows==persons:", nrow(h)==n_distinct(h$person_id), " n=", nrow(h), "\n")

d <- cx |> left_join(h, by="person_id") |> filter(!is.na(happy))
cat("\nanalytic n (cohort with happiness):", nrow(d), "\n")

## objective-health proxy (this probe): pre-baseline morbidity via log_util quartiles.
## (If it clears, replace with a proper Elixhauser index + BMI/BP.)
d <- d |> filter(!is.na(log_util))
d$health_q <- cut(d$log_util, quantile(d$log_util, 0:4/4, na.rm=TRUE), include.lowest=TRUE, labels=c("Q1","Q2","Q3","Q4"))

cat("\n== KA (well-posedness): mean happiness by objective-morbidity quartile (log_util) ==\n")
print(as.data.frame(d |> group_by(health_q) |> summarise(n=n(), happy_mean=round(mean(happy),3))))
cat("  (steep decline => SWB tracks objective health, consensus; FLAT => decoupling, the thesis)\n")

## --- absolute vs relative income ---
di <- d |> filter(!is.na(income_n), !is.na(deprivation_index))
di$inc_q  <- cut(di$income_n, quantile(di$income_n, 0:4/4, na.rm=TRUE), include.lowest=TRUE, labels=c("I1","I2","I3","I4"))
di$dep_t  <- cut(di$deprivation_index, quantile(di$deprivation_index, 0:3/3, na.rm=TRUE), include.lowest=TRUE, labels=c("dep_low","dep_mid","dep_high"))
cat("\n== KB1: mean happiness by ABSOLUTE income quartile ==\n")
print(as.data.frame(di |> group_by(inc_q) |> summarise(n=n(), happy_mean=round(mean(happy),3))))
cat("\n== KB2: mean happiness by income quartile WITHIN area-deprivation tertile (relative position) ==\n")
print(as.data.frame(di |> group_by(dep_t, inc_q) |> summarise(n=n(), happy_mean=round(mean(happy),3), .groups="drop")))
cat("  (relative-deprivation signal: same absolute income buys MORE happiness in a MORE-deprived area,\n",
    "   i.e. happiness at fixed inc_q rises from dep_low -> dep_high => rank matters beyond absolute)\n")

cat("\nDONE\n")
