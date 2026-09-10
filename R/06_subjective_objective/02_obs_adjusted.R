# R/06_subjective_objective/02_obs_adjusted.R
# Fix the D1 confound: the distinct-condition COUNT accumulates with EHR observation time, and age drives
# both -> "same rating, 2.7x morbidity by age" is partly artifact (advisor). Rebuild with observation held
# constant and an accumulation-free objective anchor (BMI), plus a bounded morbidity measure.
# TESTS:
#  T1: BMI (does NOT accumulate) by self-rated physical health WITHIN age band. If BMI mapping is STABLE
#      across age but the count moved, the age "finding" was observation artifact.
#  T2: condition count PER OBSERVATION-YEAR by rating x age (crude accumulation adjustment).
#  T3: same, restricted to comparable EHR observation (>=3y and <=10y) -> observation held roughly constant.
# Aggregate only. Individual-level.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

lik5 <- function(a) case_when(grepl('excellent|excllent',a)~5L, grepl('very good',a)~4L,
                              grepl('\\bgood\\b',a)~3L, grepl('fair',a)~2L, grepl('poor',a)~1L, TRUE~NA_integer_)
phys <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey`
           WHERE question LIKE 'Overall Health: General Physical Health'") |>
  transmute(person_id, v=lik5(ans)) |> filter(!is.na(v)) |> group_by(person_id) |> summarise(sr_phys=mean(v),.groups="drop")

## objective anchors + observation window
como <- q("SELECT person_id, COUNT(DISTINCT condition_concept_id) n_cond,
                  DATE_DIFF(MAX(condition_start_date), MIN(condition_start_date), DAY) obs_days
           FROM `__CDR__.condition_occurrence` GROUP BY person_id")
como <- como |> mutate(n_cond=as.numeric(n_cond), obs_yr=as.numeric(obs_days)/365.25)
bmi <- q("SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
          WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id")
dem <- q("SELECT person_id, year_of_birth FROM `__CDR__.person`") |>
  mutate(age=2024-as.numeric(year_of_birth), age_band=cut(age,c(0,45,65,120),labels=c("<45","45-64","65+")))

d <- dem |> inner_join(phys,by="person_id") |> left_join(como,by="person_id") |> left_join(bmi,by="person_id")

cat("== T1: BMI by self-rated physical health WITHIN age band (accumulation-free anchor) ==\n")
cat("   (if BMI mapping is STABLE across age bands, the count-based age gradient was observation artifact)\n")
print(as.data.frame(d |> filter(!is.na(bmi)) |> group_by(age_band,sr_phys) |>
  summarise(n=n(), bmi_mean=round(mean(bmi),1), .groups="drop") |> filter(n>=20)))

cat("\n== T2: condition count PER OBSERVATION-YEAR by rating x age band ==\n")
print(as.data.frame(d |> filter(!is.na(n_cond),!is.na(obs_yr),obs_yr>=1) |> mutate(cpy=n_cond/obs_yr) |>
  group_by(age_band,sr_phys) |> summarise(n=n(), cond_per_yr_med=round(median(cpy),1), .groups="drop") |> filter(n>=20)))

cat("\n== T3: raw condition count by rating x age, RESTRICTED to comparable observation (3-10y) ==\n")
print(as.data.frame(d |> filter(!is.na(n_cond),!is.na(obs_yr),obs_yr>=3,obs_yr<=10) |>
  group_by(age_band,sr_phys) |> summarise(n=n(), comorbid_med=round(median(n_cond),0), obs_med=round(median(obs_yr),1), .groups="drop") |> filter(n>=20)))

cat("\n== reference: observation years by age band (shows the confound) ==\n")
print(as.data.frame(d |> filter(!is.na(obs_yr)) |> group_by(age_band) |> summarise(n=n(), obs_yr_med=round(median(obs_yr),1))))

cat("\nDONE\n")
