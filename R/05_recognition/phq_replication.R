# R/05_recognition/phq_replication.R
# Replicate the modality finding with a validated, larger distress measure: PHQ-2 (cardinal items
# 'little interest' + 'feeling down', screen-positive at sum >= 3) instead of the placeholder PTSD item.
# Re-run the two things that matter: (1) supply -> EHR depression-coding gradient among the distressed
# (event rate + Cox HR, unadjusted vs +prior_ehr_days), and (2) the modality mechanism (frequent
# MH-professional contact -> EHR-code% by supply). If both hold at n>>3173, the finding is not a
# small-sample fluke. Aggregate only. Self-contained; run from repo root.

suppressMessages({library(bigrquery); library(dplyr); library(readr); library(survival)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

cx <- readRDS("cox_dat2.rds")
cx$zip3 <- substr(as.character(cx$zip3_as_string), 1, 3)
sup <- read_csv("zip3_mh_supply.csv", col_types = cols(zip3 = "c")) |> filter(!is.na(mh_prov_per_100k))
cx  <- left_join(cx, sup, by = "zip3")
cx$sup_t <- relevel(factor(cut(cx$mh_prov_per_100k, quantile(cx$mh_prov_per_100k, 0:3/3, na.rm = TRUE),
                include.lowest = TRUE, labels = c("low","mid","high"))), ref = "high")

## PHQ-2 cardinal items; frequency -> 0..3 (not at all / several days / more than half / nearly every day)
freq <- function(ans) case_when(grepl('not at all', ans) ~ 0L, grepl('several days', ans) ~ 1L,
                                grepl('more than half', ans) ~ 2L, grepl('nearly every day', ans) ~ 3L,
                                TRUE ~ NA_integer_)
p1 <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%little interest or pleasure%'") |>
  transmute(person_id, phq_interest = freq(ans))
p2 <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%feeling down, depressed%'") |>
  transmute(person_id, phq_down = freq(ans))
cat("== PHQ item answer mapping check (interest) ==\n"); print(table(p1$phq_interest, useNA="ifany"))

phq <- full_join(p1, p2, by = "person_id") |>
  mutate(phq2 = phq_interest + phq_down, distress = as.integer(phq2 >= 3))       # standard PHQ-2 screen +

d <- cx |> left_join(phq, by="person_id") |>
  filter(!is.na(sup_t), !is.na(distress), !is.na(time_days), time_days > 0, !is.na(event))
cat("\nPHQ cohort n:", nrow(d), " | distressed (PHQ2>=3):", sum(d$distress), "\n")

## sanity: PHQ distress -> recognition (positive control, should replicate 2.8 vs 9.9-ish)
cat("\n== recognition rate by PHQ distress ==\n")
print(as.data.frame(d |> group_by(distress) |> summarise(n=n(), events=sum(event), rate=round(mean(event),4))))

dd <- d |> filter(distress == 1)
cat("\nPHQ-distressed n:", nrow(dd), " events:", sum(dd$event), "\n")
cat("\n== (1) supply -> EHR-coding gradient among PHQ-distressed: event rate by tertile ==\n")
print(as.data.frame(dd |> group_by(sup_t) |> summarise(n=n(), events=sum(event), rate=round(mean(event),4))))
hr <- function(m,t){s<-summary(m)$coefficients; for(x in t) if(x %in% rownames(s))
  cat(sprintf("  %-10s HR %.2f (%.2f-%.2f) p=%.3f\n",x,exp(s[x,1]),exp(s[x,1]-1.96*s[x,3]),exp(s[x,1]+1.96*s[x,3]),s[x,5]))}
dd$prior_yr <- as.numeric(dd$prior_ehr_days)/365.25
cat("-- unadjusted --\n");         hr(coxph(Surv(time_days,event)~sup_t, dd), c("sup_tlow","sup_tmid"))
cat("-- + prior_ehr_days --\n");   hr(coxph(Surv(time_days,event)~sup_t+prior_yr, dd), c("sup_tlow","sup_tmid"))

## (2) modality mechanism among PHQ-distressed
vis <- q("SELECT person_id, LOWER(answer) ans FROM `__CDR__.ds_survey` WHERE question LIKE '%Mental Health Professional Visits%'") |>
  transmute(person_id, mh_hi = case_when(grepl('16 or more|10 to 12|8 to 9|6 to 7',ans)~1L, grepl('visits:',ans)~0L, TRUE~NA_integer_))
dm <- dd |> left_join(vis, by="person_id") |> filter(!is.na(mh_hi))
cat("\n== (2) modality: among frequent-MH-contact PHQ-distressed, EHR-code% by supply ==\n")
print(as.data.frame(dm |> filter(mh_hi==1) |> group_by(sup_t) |> summarise(n=n(), got_code_pct=round(100*mean(event),1))))
cat("\nDONE\n")
