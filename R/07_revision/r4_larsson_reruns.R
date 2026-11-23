# r4_larsson_reruns.R -- Johan Larsson's rerun items 2 to 6 of 23 Nov 2026 (item 1 is r2_baseline_mood.R).
#  2. survey-response weights with education and income from The Basics in the response model (eTable 5A; letter R1.1)
#  3. censoring weights with area deprivation and income in the censoring model; both interaction models refitted with the weights
#  4. affect subcohort with delayed entry: risk time starts at the wellbeing-survey date; both subcohort HRs printed
#  5. outcome codes: sensitivity excluding events whose first code is a remission or recurrent code; PMDD and schizoaffective dropped
#  6. area-mean cohesion per SD and per IQR of the area means
# Primary specification unchanged. Aggregate output only; counts below 20 masked. Needs: survival, bigrquery.
source("R/07_revision/_common.R", local = TRUE); bq()
## Only item 3 answers a reviewer request (R2.3: detection by income and area under censoring weights). Items 2, 4, 5 and 6 run only with LARSSON_ALL=1.
ALL <- Sys.getenv("LARSSON_ALL") == "1"

## ============================================================ 2. response weights with education and income
if (ALL) {
hdr("2. Survey-response weights: response model with education and income from The Basics")
pop <- run_sql(sprintf("
WITH ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id),
     resp AS (SELECT DISTINCT person_id FROM `%s.ds_survey` WHERE question_concept_id = 40192463),
     inc AS (SELECT person_id, ARRAY_AGG(answer ORDER BY survey_datetime, answer LIMIT 1)[OFFSET(0)] AS income FROM `%s.ds_survey` WHERE question_concept_id = 1585375 GROUP BY person_id),
     edu AS (SELECT person_id, ARRAY_AGG(answer ORDER BY survey_datetime, answer LIMIT 1)[OFFSET(0)] AS education FROM `%s.ds_survey` WHERE question_concept_id = 1585940 GROUP BY person_id)
SELECT p.person_id, EXTRACT(YEAR FROM p.birth_datetime) AS birth_year, g.concept_name AS gender, r.concept_name AS race, e.concept_name AS ethnicity,
       ehr.first_ehr, ehr.last_ehr, ehr.n_cond, IF(resp.person_id IS NULL, 0, 1) AS responder, inc.income, edu.education
FROM `%s.person` p JOIN ehr USING(person_id)
LEFT JOIN resp USING(person_id) LEFT JOIN inc USING(person_id) LEFT JOIN edu USING(person_id)
LEFT JOIN `%s.concept` g ON p.gender_concept_id = g.concept_id
LEFT JOIN `%s.concept` r ON p.race_concept_id = r.concept_id
LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id = e.concept_id", cdr, cdr, cdr, cdr, cdr, cdr, cdr, cdr))
pop$person_id <- as.character(pop$person_id); pop$responder <- as.integer(pop$responder)
pop$age <- 2023 - as.numeric(pop$birth_year)
pop$age_band <- cut(pop$age, c(-Inf, 29, 44, 64, Inf), labels = c("18-29","30-44","45-64","65+"))
pop$gender_c <- ifelse(pop$gender %in% c("Female","Male"), pop$gender, "Other/not reported")
pop$race_eth <- with(pop, ifelse(ethnicity == "Hispanic or Latino", "Hispanic", ifelse(race == "White", "White", ifelse(race == "Black or African American", "Black", ifelse(race == "Asian", "Asian", "Other/not reported")))))
pop$ehr_years <- as.numeric(as.Date(pop$last_ehr) - as.Date(pop$first_ehr)) / 365.25; pop$log_ncond <- log1p(as.numeric(pop$n_cond))
income_map <- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
educ_map <- c("Highest Grade: Never Attended"=1,"Highest Grade: One Through Four"=1,"Highest Grade: Five Through Eight"=1,"Highest Grade: Nine Through Eleven"=1,"Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"Highest Grade: College Graduate"=4,"Highest Grade: Advanced Degree"=4)
pop$income_c <- factor(ifelse(is.na(income_map[pop$income]), "not reported", ifelse(income_map[pop$income] <= 3, "<35k", ifelse(income_map[pop$income] <= 5, "35-75k", ">=75k"))), levels = c("<35k","35-75k",">=75k","not reported"))
pop$educ_c <- factor(ifelse(is.na(educ_map[pop$education]), "not reported", c("less than high school","high school","some college","college or more")[educ_map[pop$education]]), levels = c("less than high school","high school","some college","college or more","not reported"))
cmp <- function(v){ t <- prop.table(table(pop[[v]], pop$responder), 2) * 100; out <- data.frame(level = rownames(t), non_responders_pct = sprintf("%.1f", t[, "0"]), responders_pct = sprintf("%.1f", t[, "1"])); cat(sprintf("\n-- %s (for eTable 5, panel A) --\n", v)); print(out, row.names = FALSE) }
cmp("income_c"); cmp("educ_c")
cat(sprintf("Basics income answered: %.1f%% of non-responders, %.1f%% of responders; education: %.1f%% and %.1f%%\n", 100*mean(pop$income_c[pop$responder==0] != "not reported"), 100*mean(pop$income_c[pop$responder==1] != "not reported"), 100*mean(pop$educ_c[pop$responder==0] != "not reported"), 100*mean(pop$educ_c[pop$responder==1] != "not reported")))
fit_w <- function(form, lab){
  pm <- glm(form, data = pop, family = binomial); pop$p <- fitted(pm)
  d <- cox; d$w <- 1 / pop$p[match(d$person_id, pop$person_id)]; d <- d[!is.na(d$w), ]; d$w <- d$w / mean(d$w); cap <- quantile(d$w, .99); d$w[d$w > cap] <- cap
  cat(sprintf("%s: pseudo-R2 %.3f; weights median %.2f, 95th pct %.2f, cap %.2f; effective n %s\n", lab, 1 - pm$deviance / pm$null.deviance, median(d$w), quantile(d$w, .95), cap, mask(round(sum(d$w)^2 / sum(d$w^2)))))
  hr(coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), d, weights = w, robust = TRUE), lab = paste("  primary weighted,", lab)) }
hr(fit(RHS, cox), lab = "primary, unweighted")
fit_w(responder ~ age_band + gender_c + race_eth + ehr_years + log_ncond, "age, gender, race/ethnicity, EHR depth (as submitted)")
fit_w(responder ~ age_band + gender_c + race_eth + ehr_years + log_ncond + income_c + educ_c, "plus Basics income and education (Larsson item 2)")
}
## ============================================================ 3. censoring weights with deprivation and income; interactions refitted
hdr("3. IPCW with area deprivation and income in the censoring model; interaction models with the weights")
d <- cox[!is.na(cox$z_cohesion) & !is.na(cox$deprivation_index), ]
d$end_date <- d$expo_date + d$time_days; admin_end <- max(d$end_date)
d$lost <- as.integer(d$event == 0 & d$end_date < admin_end - 180)
d$dep_z <- as.numeric(scale(d$deprivation_index)); d$income_z <- as.numeric(scale(d$income_n)); d$income_z[is.na(d$income_z)] <- 0; d$income_miss <- as.integer(is.na(d$income_n))
ipcw <- function(rhs_c){
  cm <- coxph(as.formula(paste("Surv(time_days, lost) ~", rhs_c)), d); cm0 <- coxph(Surv(time_days, lost) ~ 1, d)
  H <- basehaz(cm, centered = TRUE); H0 <- basehaz(cm0, centered = TRUE)
  Hi <- approx(H$time, H$hazard, xout = d$time_days, method = "constant", rule = 2, f = 0)$y; H0i <- approx(H0$time, H0$hazard, xout = d$time_days, method = "constant", rule = 2, f = 0)$y
  w <- exp(-H0i) / exp(-Hi * exp(predict(cm, type = "lp"))); cap <- quantile(w, .99); w[w > cap] <- cap; list(w = w, cm = cm) }
cat(sprintf("ZIP3-available rows: %s; lost to follow-up (>180 d before data end): %s (%.1f%%)\n", mask(nrow(d)), mask(sum(d$lost)), 100*mean(d$lost)))
I0 <- ipcw(RHS); I1 <- ipcw(paste(RHS, "+ dep_z + income_z + income_miss"))
hr(I1$cm, "dep_z", "censoring model: area deprivation -> loss (per SD)"); hr(I1$cm, "income_z", "censoring model: income -> loss (per SD)")
d$w0 <- I0$w; d$w1 <- I1$w
cat(sprintf("weights with deprivation and income: median %.3f, IQR %.3f-%.3f, cap %.2f\n", median(d$w1), quantile(d$w1, .25), quantile(d$w1, .75), max(d$w1)))
hr(fit(RHS, d), lab = "primary, unweighted (ZIP3 rows)")
hr(coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), d, weights = w0, robust = TRUE), lab = "primary, IPCW as submitted")
hr(coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), d, weights = w1, robust = TRUE), lab = "primary, IPCW with deprivation and income")
fdep <- as.formula(paste("Surv(time_days, event) ~ z_cohesion * dep_z +", RHS_NOCOH)); finc <- as.formula(paste("Surv(time_days, event) ~ z_cohesion * income_z + income_miss +", sub("income_f \\+ ", "", RHS_NOCOH)))
hr(coxph(fdep, d), "z_cohesion:dep_z", "cohesion x deprivation, unweighted"); hr(coxph(fdep, d, weights = w1, robust = TRUE), "z_cohesion:dep_z", "cohesion x deprivation, IPCW-weighted (robust SE)")
hr(coxph(finc, d), "z_cohesion:income_z", "cohesion x income, unweighted"); hr(coxph(finc, d, weights = w1, robust = TRUE), "z_cohesion:income_z", "cohesion x income, IPCW-weighted (robust SE)")
## ============================================================ 4. affect subcohort with delayed entry
if (ALL) {
hdr("4. Subcohorts: affect subcohort with delayed entry at the wellbeing-survey date; wearable subcohort as before")
wb <- readRDS("wb_item_dat.rds"); wb$person_id <- as.character(wb$person_id)
wbd <- run_sql(sprintf("SELECT person_id, MIN(DATE(survey_datetime)) AS wb_date FROM `%s.ds_survey` WHERE question_concept_id = 1703980 GROUP BY person_id", cdr))
wbd$person_id <- as.character(wbd$person_id); wbd$wb_date <- as.Date(wbd$wb_date)
s <- cox[cox$person_id %in% wb$person_id[!is.na(wb$event)], ]; s$wb_date <- wbd$wb_date[match(s$person_id, wbd$person_id)]
s$entry <- as.numeric(s$wb_date - s$expo_date)
cat(sprintf("affect subcohort in the at-risk cohort: %s; wellbeing survey before the cohesion baseline: %s; same day: %s; after: %s; after the end of follow-up (excluded under delayed entry): %s\n",
            mask(nrow(s)), mask(sum(s$entry < 0, na.rm = TRUE)), mask(sum(s$entry == 0, na.rm = TRUE)), mask(sum(s$entry > 0, na.rm = TRUE)), mask(sum(s$entry >= s$time_days, na.rm = TRUE))))
hr(fit(RHS, s), lab = "affect subcohort, time from cohesion baseline (as submitted)")
s2 <- s[!is.na(s$entry) & s$entry < s$time_days, ]; s2$entry <- pmax(s2$entry, 0)
m_de <- coxph(as.formula(paste("Surv(entry, time_days, event) ~", RHS)), s2); hr(m_de, lab = "affect subcohort, delayed entry at the wellbeing-survey date")
cat(sprintf("  events after entry: %s; person-years after entry: %.0f\n", mask(sum(s2$event)), sum(s2$time_days - s2$entry) / 365.25))
fb <- readRDS("fitbit_dat.rds"); fb_ids <- as.character(fb$person_id[!is.na(fb$steps) & !is.na(fb$event)])
hr(fit(RHS, cox[cox$person_id %in% fb_ids, ]), lab = "wearable subcohort (time from cohesion baseline)")
}
## ============================================================ 5. outcome codes
if (ALL) {
hdr("5. Outcome codes: first qualifying concept per event; remission and recurrent codes; PMDD and schizoaffective")
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr, cdr, cdr)
fc <- run_sql(sprintf("WITH f AS (SELECT person_id, MIN(condition_start_date) AS d FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id)
SELECT co.person_id, STRING_AGG(DISTINCT c.concept_name, ' | ' ORDER BY c.concept_name) AS first_names FROM `%s.condition_occurrence` co JOIN f ON co.person_id = f.person_id AND co.condition_start_date = f.d JOIN `%s.concept` c ON c.concept_id = co.condition_concept_id WHERE co.condition_concept_id IN (%s) GROUP BY co.person_id", cdr, dep_ids, cdr, cdr, dep_ids))
fc$person_id <- as.character(fc$person_id)
d <- cox; d$first_names <- fc$first_names[match(d$person_id, fc$person_id)]
ev <- d[d$event == 1, ]
cat(sprintf("events with a first-code name: %s of %s\n", mask(sum(!is.na(ev$first_names))), mask(nrow(ev))))
flag <- function(pat) grepl(pat, ev$first_names, ignore.case = TRUE)
rem <- flag("remission"); rec <- flag("recurrent"); pmdd <- flag("premenstrual"); sza <- flag("schizoaffective")
cat(sprintf("first code includes: remission %s; recurrent %s; remission or recurrent %s; premenstrual dysphoric %s; schizoaffective %s\n", mask(sum(rem)), mask(sum(rec)), mask(sum(rem | rec)), mask(sum(pmdd)), mask(sum(sza))))
cat("  most frequent first-code names among remission/recurrent events:\n"); tb <- sort(table(unlist(strsplit(ev$first_names[rem | rec], " \\| "))), decreasing = TRUE); for (k in head(names(tb), 12)) cat(sprintf("    %-70s %s\n", k, mask(tb[[k]])))
recode <- function(drop_ids, lab, censor = TRUE){
  dd <- d; hit <- dd$person_id %in% drop_ids
  if (censor) dd$event[hit] <- 0L else dd <- dd[!hit, ]
  hr(fit(RHS, dd), lab = sprintf("%s (%s)", lab, if (censor) "events recoded as censored at the event date" else "participants dropped")) }
hr(fit(RHS, d), lab = "primary (Table 2)")
recode(ev$person_id[rem | rec], "excluding first codes of remission or recurrent depression"); recode(ev$person_id[rem | rec], "excluding first codes of remission or recurrent depression", censor = FALSE)
recode(ev$person_id[pmdd | sza], "dropping PMDD and schizoaffective first codes")
recode(ev$person_id[rem | rec | pmdd | sza], "all of the above")
}
## ============================================================ 6. area-mean cohesion per SD and per IQR of the area means
if (ALL) {
hdr("6. Area-mean (contextual) cohesion rescaled per SD and per IQR of the area means")
d <- cox[!is.na(cox$zip3_as_string) & !is.na(cox$z_cohesion), ]; d$clus <- factor(d$zip3_as_string)
S <- tapply(d$z_cohesion, d$clus, sum); N <- tapply(d$z_cohesion, d$clus, length)
d$area_mean <- as.numeric((S[d$clus] - d$z_cohesion) / pmax(N[d$clus] - 1, 1)); d$area_mean[N[d$clus] == 1] <- NA
dd <- d[!is.na(d$area_mean), ]
sd_a <- sd(dd$area_mean); iqr_a <- IQR(dd$area_mean)
cat(sprintf("area means: SD across participants %.3f; IQR %.3f (individual score has SD 1)\n", sd_a, iqr_a))
dd$area_sd <- dd$area_mean / sd_a; dd$area_iqr <- dd$area_mean / iqr_a
mb <- fit(paste(RHS, '+ area_mean + frailty(clus, distribution = "gaussian")'), dd); hr(mb, "area_mean", "area mean per 1 SD of the individual scale (as submitted)")
mb_sd <- fit(paste(RHS, '+ area_sd + frailty(clus, distribution = "gaussian")'), dd); hr(mb_sd, "area_sd", "area mean per 1 SD of the area means (for the Results paragraph)")
hr(mb_sd, "z_cohesion", "  individual cohesion in the same model")
mb_iqr <- fit(paste(RHS, '+ area_iqr + frailty(clus, distribution = "gaussian")'), dd); hr(mb_iqr, "area_iqr", "area mean per IQR of the area means")
cat("\nSentence for the Results: replace '(HR, 0.77; 95% CI, 0.64-0.93)' with the per-SD estimate above, stated as 'per SD of the area means'.\n")
}
