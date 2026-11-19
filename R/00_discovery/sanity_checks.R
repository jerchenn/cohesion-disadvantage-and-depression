# sanity_checks.R -- one pass over cox_dat3.rds looking for errors of the education-recode kind:
# a variable that is silently degenerate, a level in the wrong direction, a sample that differs
# between scripts without being said, a time variable that cannot be what it claims.
# Prints aggregates only (any cell below 20 is shown as "<20"). Reads cox_dat3.rds; one BigQuery
# query for the depression concept list. Needs: survival, bigrquery.

library(survival)
d <- readRDS("cox_dat3.rds")
mask <- function(x){ x <- as.numeric(x); ifelse(!is.na(x) & x > 0 & x < 20, "<20", format(x, big.mark=",", trim=TRUE)) }
tab  <- function(x, lab){ t <- table(x, useNA="ifany"); cat(sprintf("\n-- %s --\n", lab))
  print(data.frame(level=names(t), n=mask(t), pct=sprintf("%.1f", 100*as.numeric(t)/length(x))), row.names=FALSE) }
hdr  <- function(s) cat("\n==========", s, "==========\n")

## 1. integrity ---------------------------------------------------------------------------------
hdr("1. Integrity of the at-risk file")
cat(sprintf("rows %s | unique person_id %s | events %s (%.2f%%)\n", mask(nrow(d)), mask(length(unique(d$person_id))),
            mask(sum(d$event)), 100*mean(d$event)))
cat(sprintf("event values: %s | time_days <= 0: %s | time_days NA: %s\n",
            paste(sort(unique(d$event)), collapse=","), mask(sum(d$time_days <= 0, na.rm=TRUE)), mask(sum(is.na(d$time_days)))))
# percentiles only: a minimum or maximum is a value from one participant and is not printed
qf <- function(x, p) quantile(x, p, na.rm=TRUE)
cat(sprintf("follow-up years: 1st pct %.2f, median %.2f, 99th pct %.2f | over 3.5 y: %s\n",
            qf(d$time_days,.01)/365.25, median(d$time_days)/365.25, qf(d$time_days,.99)/365.25, mask(sum(d$time_days > 3.5*365.25))))
cat(sprintf("baseline (expo_date): 1st pct %s, 99th pct %s | before 2021-11: %s | after 2024-12: %s\n",
            qf(as.numeric(d$expo_date),.01) |> as.Date(origin="1970-01-01"), qf(as.numeric(d$expo_date),.99) |> as.Date(origin="1970-01-01"),
            mask(sum(d$expo_date < as.Date("2021-11-01"))), mask(sum(d$expo_date > as.Date("2024-12-31")))))
cat(sprintf("age: 1st pct %.1f, median %.1f, 99th pct %.1f | under 18: %s | over 100: %s | NA: %s\n",
            qf(d$age,.01), median(d$age,na.rm=TRUE), qf(d$age,.99),
            mask(sum(d$age < 18, na.rm=TRUE)), mask(sum(d$age > 100, na.rm=TRUE)), mask(sum(is.na(d$age)))))
cat(sprintf("prior_ehr_days: 1st pct %.0f, median %.0f | below 365 (should be none): %s\n",
            qf(d$prior_ehr_days,.01), median(d$prior_ehr_days), mask(sum(d$prior_ehr_days < 365))))
cat(sprintf("n_cond: zero %s | median %.0f | log_util 1st-99th pct %.2f to %.2f\n",
            mask(sum(d$n_cond == 0)), median(d$n_cond), qf(d$log_util,.01), qf(d$log_util,.99)))
cat("\nevents by year of baseline (an error in date handling shows up as a year with no events):\n")
yr <- format(d$expo_date, "%Y"); print(data.frame(year=names(table(yr)), n=mask(table(yr)),
     events=mask(tapply(d$event, yr, sum)), fu_median_y=sprintf("%.2f", tapply(d$time_days, yr, median)/365.25)), row.names=FALSE)

## 2. every covariate: levels, counts, missingness -------------------------------------------------
hdr("2. Covariate levels in the at-risk cohort (a degenerate variable has one level or mostly NA)")
for (v in c("sex","race","ethnicity","income_f","educ_f","emp_f")) tab(d[[v]], v)
cat("\nNA count per model variable (this is what drops 98,795 to the complete-case n):\n")
mv <- c("z_cohesion","age","sex","race","ethnicity","income_f","educ_f","emp_f","log_util","time_days","event")
print(data.frame(variable=mv, n_NA=mask(sapply(mv, function(v) sum(is.na(d[[v]]))))), row.names=FALSE)
cc <- complete.cases(d[, mv]); cat(sprintf("complete cases: %s of %s\n", mask(sum(cc)), mask(nrow(d))))

## 3. cohesion score ------------------------------------------------------------------------------
hdr("3. Cohesion score")
items <- c("help","getalong","trust","values","watchout")
cat(sprintf("z_cohesion in at-risk cohort: mean %.3f, SD %.3f (standardised on all responders, so not 0/1 here)\n",
            mean(d$z_cohesion, na.rm=TRUE), sd(d$z_cohesion, na.rm=TRUE)))
cat("item correlation matrix (should be all positive, roughly 0.4 to 0.8):\n")
print(round(cor(d[, items], use="pairwise"), 2))
X <- scale(d[, items]); k <- ncol(X); a <- k/(k-1) * (1 - sum(apply(X,2,var,na.rm=TRUE)) / var(rowSums(X, na.rm=TRUE)))
cat(sprintf("Cronbach alpha (standardised items, complete rows): %.3f (manuscript says 0.88)\n", a))
cat("crude event rate by cohesion quintile (should fall monotonically):\n")
q5 <- cut(d$z_cohesion, quantile(d$z_cohesion, 0:5/5, na.rm=TRUE), include.lowest=TRUE, labels=paste0("Q",1:5))
print(data.frame(quintile=levels(q5), n=mask(table(q5)), event_pct=sprintf("%.2f", 100*tapply(d$event, q5, mean))), row.names=FALSE)
cat(sprintf("z_cohesion NA: %s (respondents with fewer than 3 of 5 items; this alone is the 98,795 -> 96,581 drop)\n", mask(sum(is.na(d$z_cohesion)))))
na_c <- is.na(d$z_cohesion)
cat("participants without a cohesion score vs with one (is the exclusion ignorable?):\n")
print(data.frame(group=c("no score","score"), n=mask(c(sum(na_c), sum(!na_c))),
      event_pct=sprintf("%.2f", 100*c(mean(d$event[na_c]), mean(d$event[!na_c]))),
      age_mean=sprintf("%.1f", c(mean(d$age[na_c]), mean(d$age[!na_c]))),
      female_pct=sprintf("%.1f", 100*c(mean(d$sex[na_c]=="Female"), mean(d$sex[!na_c]=="Female"))),
      income_NA_pct=sprintf("%.1f", 100*c(mean(is.na(d$income_n[na_c])), mean(is.na(d$income_n[!na_c])))),
      BAplus_pct=sprintf("%.1f", 100*c(mean(d$educ_n[na_c]==4, na.rm=TRUE), mean(d$educ_n[!na_c]==4, na.rm=TRUE))),
      n_cond_median=c(median(d$n_cond[na_c]), median(d$n_cond[!na_c]))), row.names=FALSE)

## 4. crude gradients that must run the expected way ----------------------------------------------
hdr("4. Crude event rates by covariate level (direction check; a reversed recode shows here)")
crude <- function(v, lab){ t <- table(d[[v]], useNA="ifany"); e <- sapply(names(t), function(l) if (is.na(l)) mean(d$event[is.na(d[[v]])]) else mean(d$event[!is.na(d[[v]]) & d[[v]] == l]))
  cat(sprintf("\n-- %s --\n", lab)); print(data.frame(level=names(t), n=mask(t), event_pct=sprintf("%.2f", 100*e)), row.names=FALSE) }
crude("educ_f", "education 1=<HS 2=HS/GED 3=some college 4=BA+ (rate should fall with level)")
crude("income_f", "income 1=<10k ... 9=>200k (rate should fall with level)")
crude("emp_f", "employment (Unable To Work should be highest)")
crude("sex", "sex (female higher)")
ad <- cut(d$age, c(0,30,40,50,60,70,80,120)); crude_ad <- tapply(d$event, ad, mean)
cat("\n-- age band --\n"); print(data.frame(band=names(crude_ad), n=mask(table(ad)), event_pct=sprintf("%.2f", 100*crude_ad)), row.names=FALSE)

## 5. the primary model, coefficient by coefficient ------------------------------------------------
hdr("5. Primary model: every coefficient (a level with HR far from 1 and a huge CI is a sparse or broken level)")
rhs <- "z_cohesion + age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
m <- coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", rhs)), d)
s <- summary(m)$coefficients; ci <- suppressMessages(exp(confint(m)))
print(data.frame(term=rownames(s), HR=sprintf("%.3f", exp(s[,"coef"])), lo=sprintf("%.3f", ci[,1]), hi=sprintf("%.3f", ci[,2]),
                 p=signif(s[,"Pr(>|z|)"], 2)), row.names=FALSE)
cat(sprintf("n=%s events=%s\n", mask(m$n), mask(m$nevent)))

## 6. the three "primary" numbers printed by three scripts ----------------------------------------
hdr("6. Why cox_model3.R, mi_sensitivity.R and etable_robustness.R print different primary HRs")
m0 <- coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", sub(" \\+ log_util","",rhs))), d)
cat(sprintf("factor spec, with log_util:    HR %.5f (%.5f-%.5f)  coef on log_util = %.4f (p=%.2g)\n",
            exp(coef(m)["z_cohesion"]), exp(confint(m)["z_cohesion",1]), exp(confint(m)["z_cohesion",2]),
            coef(m)["log_util"], s["log_util","Pr(>|z|)"]))
cat(sprintf("factor spec, without log_util: HR %.5f (%.5f-%.5f)\n",
            exp(coef(m0)["z_cohesion"]), exp(confint(m0)["z_cohesion",1]), exp(confint(m0)["z_cohesion",2])))
r <- cor(d$z_cohesion, d$log_util, use="complete"); dd <- d[!is.na(d$z_cohesion), ]; ra <- cor(resid(lm(z_cohesion ~ age, dd)), resid(lm(log_util ~ age, dd)))
cat(sprintf("corr(z_cohesion, log_util) = %.3f; partialling out age = %.3f  (near zero explains identical HRs)\n", r, ra))
dc <- d[!is.na(d$zip3_as_string), ]
mc <- coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", rhs)), dc)
cat(sprintf("ZIP3-available subset: n=%s of %s, HR %.4f (the cluster-robust line in cox_model3.R is this sample, not the full cohort)\n",
            mask(mc$n), mask(m$n), exp(coef(mc)["z_cohesion"])))
d$sex_c <- factor(ifelse(d$sex %in% c("Female","Male"), as.character(d$sex), "Other"))
d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), as.character(d$race), "Other"))
d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic","NotHispanic"))
d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
mh <- coxph(Surv(time_days,event) ~ z_cohesion + age + sex_c + race_c + ethn_c + income_m + educ_m + log_util, d)
cat(sprintf("harmonised/collapsed spec (figs_data.R, etable_robustness.R footnote, Figure 3): HR %.4f (%.4f-%.4f) n=%s\n",
            exp(coef(mh)["z_cohesion"]), exp(confint(mh)["z_cohesion",1]), exp(confint(mh)["z_cohesion",2]), mask(mh$n)))
cat("mi_sensitivity.R's 'missing-indicator primary HR 0.883' is a hard-coded string in that script, not a computed value.\n")

## 7. the outcome: which depression concepts actually produce the events --------------------------
hdr("7. Depression concept set: codes that produce the incident events (Reviewer 1 asked for this list)")
library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr, cdr, cdr)   # identical to the outcome set in build_cox3.R
cat("concepts in the outcome set:", nrow(run_sql(sprintf("SELECT DISTINCT descendant_concept_id FROM (%s)", dep_ids))), "\n")
fc <- run_sql(sprintf("
WITH dep AS (SELECT person_id, condition_start_date, condition_concept_id FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s)),
     first AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM dep GROUP BY person_id)
SELECT DISTINCT dep.person_id, c.concept_name, c.concept_code
FROM dep JOIN first USING(person_id) JOIN `%s.concept` c ON c.concept_id = dep.condition_concept_id
WHERE dep.condition_start_date = first.first_dep", cdr, dep_ids, cdr))
ev_ids <- d$person_id[d$event == 1]
fc$person_id <- as.character(fc$person_id); ev_ids <- as.character(ev_ids)
inc <- fc[fc$person_id %in% ev_ids, ]
tab <- aggregate(person_id ~ concept_name + concept_code, inc, function(x) length(unique(x)))
tab <- tab[order(-tab$person_id), ]; names(tab)[3] <- "n_events"
cat(sprintf("Codes at the FIRST depression record among the %s incident events in the cohort (a person can carry more than one code that day); cells <20 masked:\n", mask(length(ev_ids))))
tab$n_events <- mask(tab$n_events); print(as.data.frame(tab), row.names = FALSE, right = FALSE)
cat("\nAll concepts in the outcome set that appear anywhere in the CDR at a first record (top 60):\n")
allt <- aggregate(person_id ~ concept_name + concept_code, fc, function(x) length(unique(x)))
allt <- allt[order(-allt$person_id), ][1:min(60, nrow(allt)), ]; allt$person_id <- mask(allt$person_id)
print(as.data.frame(allt), row.names = FALSE, right = FALSE)
cat("\nLook for concepts that are not depression (bipolar episodes, adjustment disorder, grief, perinatal screening, mixed anxiety-depression).\n")
cat("The washout excludes anyone with ANY of these codes before baseline; that step removes ~36% of the EHR-observed sample.\n")

## 8. Cronbach alpha on ALL cohesion respondents (the manuscript's 0.88 may be this sample) -----------
hdr("8. Cronbach alpha, all respondents vs at-risk cohort")
tryCatch({
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer, survey_datetime FROM `%s.ds_survey` WHERE question_concept_id IN (40192463,40192411,40192499,40192417,40192400)", cdr))
sv$qid <- as.character(sv$question_concept_id); sv$answer <- as.character(sv$answer)
sv <- sv[order(sv$person_id, sv$qid, sv$survey_datetime), ]
sv <- sv[!duplicated(sv[c("person_id","qid")]), ]
labs8 <- c("40192463"="help","40192411"="getalong","40192499"="trust","40192417"="values","40192400"="watchout")
sv$item <- unname(labs8[sv$qid])
tb <- sort(table(sv$answer), decreasing=TRUE); cat("answer strings seen (counts, <20 masked):\n"); print(data.frame(answer=names(tb), n=mask(as.numeric(tb)), row.names=NULL))
m5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
m4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
sv$v <- ifelse(sv$item == "watchout", unname(m4[sv$answer]), unname(m5[sv$answer]))
W <- tidyr::pivot_wider(sv[!is.na(sv$item), c("person_id","item","v")], names_from=item, values_from=v)
Xa <- as.matrix(W[, items])
alpha_of <- function(X){ X <- scale(X); k <- ncol(X); k/(k-1) * (1 - sum(apply(X,2,var)) / var(rowSums(X))) }
alpha_raw <- function(X){ k <- ncol(X); k/(k-1) * (1 - sum(apply(X,2,var)) / var(rowSums(X))) }
cc5 <- complete.cases(Xa)
cat(sprintf("all respondents with all 5 items: n=%s, alpha=%.3f standardised, %.3f raw\n", mask(sum(cc5)),
            if (sum(cc5) > 1) alpha_of(Xa[cc5, , drop=FALSE]) else NA, if (sum(cc5) > 1) alpha_raw(Xa[cc5, , drop=FALSE]) else NA))
Xc <- as.matrix(d[, items]); ccc <- complete.cases(Xc)
cat(sprintf("at-risk cohort with all 5 items: n=%s, alpha=%.3f standardised, %.3f raw\n", mask(sum(ccc)), alpha_of(Xc[ccc, ]), alpha_raw(Xc[ccc, ])))
}, error=function(e) cat("section 8 failed:", conditionMessage(e), "\n"))
