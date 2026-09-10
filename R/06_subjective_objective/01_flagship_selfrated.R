# R/06_subjective_objective/01_flagship_selfrated.R
# Flagship pair for the SSM discordance study: SELF-RATED HEALTH (Overall Health survey, n~680k) vs
# OBJECTIVE clinical anchors. Establish (1) real overlap n, (2) the answer-scale mappings (PRINTED -- verify
# before trusting), (3) core concordance for each self-rating x objective pair, (4) ONE discordance-patterning
# look (by age + race) as proof of concept. Individual-level. Aggregate only. Dedupe per person.

suppressMessages({library(bigrquery); library(dplyr)})
CDR <- Sys.getenv("WORKSPACE_CDR")
q   <- function(sql) bq_table_download(
  bq_project_query(Sys.getenv("GOOGLE_PROJECT"), gsub("__CDR__", CDR, sql, fixed = TRUE)), bigint = "character")

## ---- SUBJECTIVE: 4 self-rated items; print raw distributions, then map (higher = better health) ----
sr <- q("SELECT person_id, question, LOWER(answer) ans FROM `__CDR__.ds_survey`
         WHERE question LIKE 'Overall Health: General Physical Health'
            OR question LIKE 'Overall Health: General Mental Health'
            OR question LIKE 'Overall Health: General Health'
            OR question LIKE 'Overall Health: Average Pain 7 Days'")
cat("== raw answer distributions (LOCK the scales) ==\n")
for (qq in unique(sr$question)) { cat("\n--", qq, "--\n"); print(sort(table(sr$ans[sr$question==qq]), decreasing=TRUE)) }

lik5 <- function(a) case_when(grepl('excellent',a)~5L, grepl('very good',a)~4L, grepl('^good| good$|\\bgood\\b',a)~3L,
                              grepl('fair',a)~2L, grepl('poor',a)~1L, TRUE~NA_integer_)   # higher=better
mk <- function(qname, fn) sr |> filter(question==qname) |> transmute(person_id, v=fn(ans)) |>
        filter(!is.na(v)) |> group_by(person_id) |> summarise(v=mean(v), .groups="drop")
phys <- mk("Overall Health: General Physical Health", lik5) |> rename(sr_phys=v)
ment <- mk("Overall Health: General Mental Health",   lik5) |> rename(sr_ment=v)
genh <- mk("Overall Health: General Health",          lik5) |> rename(sr_gen =v)
# pain is 0-10 (numeric-ish); map digits, higher = MORE pain (worse)
painfn <- function(a){ x<-suppressWarnings(as.integer(gsub("[^0-9]","",substr(a,1,3)))); ifelse(is.na(x)|x>10,NA_integer_,x) }
pain <- mk("Overall Health: Average Pain 7 Days", painfn) |> rename(sr_pain=v)
cat("\nmapping check phys:\n"); print(round(quantile(phys$sr_phys,c(0,.25,.5,.75,1)),2))
cat("mapping check pain (0-10):\n"); print(round(quantile(pain$sr_pain,c(0,.25,.5,.75,1)),2))

## ---- OBJECTIVE anchors ----
como <- q("SELECT person_id, COUNT(DISTINCT condition_concept_id) n_cond FROM `__CDR__.condition_occurrence` GROUP BY person_id")
como$n_cond <- as.numeric(como$n_cond)
bmi <- q("SELECT person_id, AVG(value_as_number) bmi FROM `__CDR__.measurement`
          WHERE measurement_concept_id=3038553 AND value_as_number BETWEEN 12 AND 80 GROUP BY person_id")
flag <- function(anc) q(sprintf("SELECT DISTINCT co.person_id FROM `__CDR__.condition_occurrence` co
  JOIN `__CDR__.concept_ancestor` ca ON ca.descendant_concept_id=co.condition_concept_id
  WHERE ca.ancestor_concept_id=%d", anc))
dep <- flag(440383)  |> mutate(dx_dep=1L)
anx <- flag(442077)  |> mutate(dx_anx=1L)
pnx <- flag(4329041) |> mutate(dx_pain=1L)

## demographics (whole cohort) from person
dem <- q("SELECT person_id, year_of_birth, race_concept_id, gender_concept_id FROM `__CDR__.person`")
dem$age <- 2024 - as.numeric(dem$year_of_birth)

## ---- assemble ----
d <- dem |>
  left_join(phys,by="person_id") |> left_join(ment,by="person_id") |> left_join(genh,by="person_id") |>
  left_join(pain,by="person_id") |> left_join(como,by="person_id") |> left_join(bmi,by="person_id") |>
  left_join(dep,by="person_id") |> left_join(anx,by="person_id") |> left_join(pnx,by="person_id") |>
  mutate(across(c(dx_dep,dx_anx,dx_pain), ~ifelse(is.na(.),0L,.)))
cat("\n== overlap n: self-rated PHYS x comorbidity ==\n")
cat("  sr_phys present:", sum(!is.na(d$sr_phys)), " | with comorbidity:", sum(!is.na(d$sr_phys)&!is.na(d$n_cond)),
    " | with BMI:", sum(!is.na(d$sr_phys)&!is.na(d$bmi)), "\n")

## ---- CONCORDANCE (core pairs) ----
cat("\n== C1: objective comorbidity count by self-rated PHYSICAL health (expect decline if concordant) ==\n")
print(as.data.frame(d |> filter(!is.na(sr_phys),!is.na(n_cond)) |> group_by(sr_phys) |>
  summarise(n=n(), comorbid_med=round(median(n_cond),0), bmi_mean=round(mean(bmi,na.rm=TRUE),1))))
cat("\n== C2: EHR depression/anxiety dx rate by self-rated MENTAL health ==\n")
print(as.data.frame(d |> filter(!is.na(sr_ment)) |> group_by(sr_ment) |>
  summarise(n=n(), dep_rate=round(mean(dx_dep),3), anx_rate=round(mean(dx_anx),3))))
cat("\n== C3: EHR pain dx rate by self-rated PAIN (0-10) ==\n")
print(as.data.frame(d |> filter(!is.na(sr_pain)) |> mutate(pain_g=cut(sr_pain,c(-1,0,3,6,10),labels=c("0","1-3","4-6","7-10"))) |>
  group_by(pain_g) |> summarise(n=n(), pain_dx_rate=round(mean(dx_pain),3))))

## ---- DISCORDANCE PATTERNING (proof of concept): comorbidity by self-rated phys WITHIN age band ----
cat("\n== D1: comorbidity by self-rated PHYSICAL health, WITHIN age band (does the mapping differ by age?) ==\n")
dd <- d |> filter(!is.na(sr_phys),!is.na(n_cond)) |>
  mutate(age_band=cut(age,c(0,45,65,120),labels=c("<45","45-64","65+")))
print(as.data.frame(dd |> group_by(age_band,sr_phys) |> summarise(n=n(), comorbid_med=round(median(n_cond),0), .groups="drop") |> filter(n>=20)))

cat("\nDONE\n")
