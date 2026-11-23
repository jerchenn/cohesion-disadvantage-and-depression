# r6_mood_pathway.R -- four analyses on the attenuation by baseline mood (r2_baseline_mood.R, r5_movers.R), 23 Nov 2026:
#  1. does cohesion predict LATER affect (EHHWB survey, answered after the cohesion survey) given the enrolment mood items?
#  2. cohesion -> incident depression among participants who were well at enrolment (mental health excellent/very good, emotional problems never/rarely)
#  3. attenuation by the age of the mood rating (same day, under 1 year, 1 to 2 years, 2+ years before the cohesion survey)
#  4. pooled estimate adjusted for mood rated 1+ years before baseline, movers and stayers together, with and without the 180-day lag
#  5. antidepressant initiation (parallel outcome) adjusted for prior mood, if antidep_dat.rds is present
# Aggregate output only. Needs: survival, bigrquery; files cox_dat3.rds and wb_item_dat.rds.
source("R/07_revision/_common.R", local = TRUE); bq()

## ---- enrolment mood items (as in r2/r5) ----
q_mh <- "Overall Health: General Mental Health"; q_em <- "Overall Health: Emotional Problem 7 Days"
sv <- run_sql(sprintf("SELECT person_id, question, LOWER(answer) AS ans, DATE(survey_datetime) AS dt FROM `%s.ds_survey` WHERE question IN ('%s','%s')", cdr, q_mh, q_em))
sv$person_id <- as.character(sv$person_id); sv$dt <- as.Date(sv$dt)
map_mh <- function(a) ifelse(grepl("excellent|excllent", a), 5L, ifelse(grepl("very good", a), 4L, ifelse(grepl("\\bgood\\b", a), 3L, ifelse(grepl("fair", a), 2L, ifelse(grepl("poor", a), 1L, NA_integer_)))))
map_em <- function(a) ifelse(grepl("never", a), 0L, ifelse(grepl("rarely", a), 1L, ifelse(grepl("sometimes", a), 2L, ifelse(grepl("often", a), 3L, ifelse(grepl("always", a), 4L, NA_integer_)))))
sv$v <- ifelse(sv$question == q_mh, map_mh(sv$ans), map_em(sv$ans)); sv <- sv[!is.na(sv$v), ]
d <- cox
pick <- function(qq, nm){ x <- sv[sv$question == qq & sv$person_id %in% d$person_id, ]; x$expo <- d$expo_date[match(x$person_id, d$person_id)]; x <- x[x$dt <= x$expo, ]
  x <- x[order(x$person_id, x$dt, decreasing = TRUE), ]; x <- x[!duplicated(x$person_id), ]; d[[nm]] <<- x$v[match(d$person_id, x$person_id)]; d[[paste0(nm, "_dt")]] <<- x$dt[match(d$person_id, x$person_id)] }
pick(q_mh, "sr_mental"); pick(q_em, "emo_prob")
d$mood_years <- as.numeric(d$expo_date - d$sr_mental_dt) / 365.25
db <- d[!is.na(d$sr_mental) & !is.na(d$emo_prob), ]
db$sr_mental_f <- relevel(factor(db$sr_mental), ref = "5"); db$emo_prob_f <- relevel(factor(db$emo_prob), ref = "0")
MOOD <- "+ sr_mental_f + emo_prob_f"
cat(sprintf("cohort members with both enrolment mood items on/before baseline: %s\n", mask(nrow(db))))

## ---- 1. later affect given earlier mood ----
hdr("1. Does cohesion predict later affect (EHHWB) given the enrolment mood items?")
wb <- readRDS("wb_item_dat.rds"); wb$person_id <- as.character(wb$person_id); wb_names <- attr(wb, "wb_names")
cat("affect items (standardised, higher = better):", paste(wb_names, collapse = ", "), "\n")
wbd <- run_sql(sprintf("SELECT person_id, MIN(DATE(survey_datetime)) AS wb_date FROM `%s.ds_survey` WHERE question_concept_id = 1703980 GROUP BY person_id", cdr))
wbd$person_id <- as.character(wbd$person_id); wbd$wb_date <- as.Date(wbd$wb_date)
wb$later_affect <- rowMeans(wb[, wb_names, drop = FALSE], na.rm = TRUE); wb$wb_date <- wbd$wb_date[match(wb$person_id, wbd$person_id)]
a <- merge(db, wb[, c("person_id", "later_affect", "wb_date")], by = "person_id")
a$affect_after <- as.numeric(a$wb_date - a$expo_date)
cat(sprintf("members with both: %s; affect survey after the cohesion survey: %s (median %.0f days after); before: %s\n", mask(nrow(a)), mask(sum(a$affect_after > 0, na.rm = TRUE)), median(a$affect_after[a$affect_after > 0], na.rm = TRUE), mask(sum(a$affect_after <= 0, na.rm = TRUE))))
a2 <- a[!is.na(a$affect_after) & a$affect_after > 0 & !is.na(a$later_affect), ]; a2$later_z <- as.numeric(scale(a2$later_affect))
COV <- sub("z_cohesion \\+ ", "", RHS)
lmco <- function(form, lab){ m <- lm(as.formula(form), a2); b <- coef(summary(m))["z_cohesion", ]; cat(sprintf("%-78s cohesion coef %.3f (SE %.3f) P=%.2g; R2 %.3f; n=%s\n", lab, b[1], b[2], b[4], summary(m)$r.squared, mask(nobs(m)))) ; m }
m_a <- lmco(paste("later_z ~ z_cohesion +", COV), "later affect (SD units) on cohesion, covariates only")
m_b <- lmco(paste("later_z ~ z_cohesion", MOOD, "+", COV), "  plus enrolment mood items")
m_c <- lm(as.formula(paste("later_z ~", sub("^\\+ ", "", MOOD), "+", COV)), a2)
cat(sprintf("  R2 with enrolment mood but without cohesion: %.3f; added by cohesion: %.4f\n", summary(m_c)$r.squared, summary(m_b)$r.squared - summary(m_c)$r.squared))
a3 <- a2[a2$mood_years >= 1, ]
if (nrow(a3) > 500) { m_d <- lm(as.formula(paste("later_z ~ z_cohesion", MOOD, "+", COV)), a3); b <- coef(summary(m_d))["z_cohesion", ]; cat(sprintf("%-78s cohesion coef %.3f (SE %.3f) P=%.2g; n=%s\n", "  same, enrolment mood rated 1+ years before the cohesion survey", b[1], b[2], b[4], mask(nobs(m_d)))) }
cat("Reading: a cohesion coefficient that stays positive given prior mood means cohesion predicts where mood goes next, which a pathway requires; a coefficient near zero means prior mood carries the association.\n")

## ---- 2. well at enrolment ----
hdr("2. Cohesion -> incident depression among participants well at enrolment")
well <- db$sr_mental >= 4 & db$emo_prob <= 1
cat(sprintf("well at enrolment (mental health excellent/very good AND emotional problems never/rarely): %s (%.1f%%), events %s; not well: %s, events %s\n", mask(sum(well)), 100*mean(well), mask(sum(db$event[well])), mask(sum(!well)), mask(sum(db$event[!well]))))
hr(fit(RHS, db[well, ]), lab = "well at enrolment: primary specification")
hr(fit(paste(RHS, MOOD), db[well, ]), lab = "well at enrolment: + mood items (within-stratum residual)")
hr(fit(RHS, db[!well, ]), lab = "not well at enrolment: primary specification")
hr(fit(paste(RHS, MOOD), db[!well, ]), lab = "not well at enrolment: + mood items")
hr(fit(RHS, db[db$sr_mental == 5 & db$emo_prob == 0, ]), lab = "excellent mental health and never bothered (strictest stratum)")

## ---- 3. attenuation by age of the mood rating ----
hdr("3. Attenuation by the age of the enrolment mood rating")
db$band <- cut(db$mood_years, c(-Inf, 0, 1, 2, Inf), labels = c("same day", "under 1 year", "1 to 2 years", "2+ years"), right = TRUE)
for (b in levels(db$band)) { dd <- db[db$band == b, ]; if (nrow(dd) < 1000) next
  m0 <- fit(RHS, dd); m1 <- fit(paste(RHS, MOOD), dd); hr(m0, lab = sprintf("%-14s primary", b)); hr(m1, lab = sprintf("%-14s + mood items", b))
  cat(sprintf("    share of the log HR removed: %.0f%%\n", 100 * (1 - coef(m1)[["z_cohesion"]] / coef(m0)[["z_cohesion"]]))) }

## ---- 4. pooled prior-mood estimate ----
hdr("4. Pooled estimate adjusted for mood rated 1+ years before the cohesion survey (movers and stayers)")
p1 <- db[db$mood_years >= 1, ]
hr(fit(RHS, p1), lab = "mood rated 1+ years before: primary")
m_p <- fit(paste(RHS, MOOD), p1); hr(m_p, lab = "mood rated 1+ years before: + mood items  <- the number for the letter")
hr(fit(paste(RHS, MOOD), p1[p1$time_days > 180, ]), lab = "  same, with the 180-day lag period")
hr(fit(paste(RHS, MOOD), p1[p1$time_days > 365, ]), lab = "  same, with a 365-day lag period")
cat(sprintf("E-value for the prior-mood-adjusted HR: %.2f (CI limit %.2f)\n", { h <- 1 / exp(coef(m_p)[["z_cohesion"]]); h + sqrt(h * (h - 1)) }, { u <- 1 / exp(coef(m_p)[["z_cohesion"]] + 1.96 * sqrt(vcov(m_p)["z_cohesion", "z_cohesion"])); if (u > 1) u + sqrt(u * (u - 1)) else 1 }))

## ---- 5. antidepressant outcome ----
hdr("5. Antidepressant initiation adjusted for enrolment mood (if the file is present)")
if (file.exists("antidep_dat.rds")) {
  ad <- readRDS("antidep_dat.rds"); ad$person_id <- as.character(ad$person_id)
  ad <- merge(ad, db[, c("person_id", "sr_mental_f", "emo_prob_f", "mood_years")], by = "person_id")
  hr(fit(RHS, ad), lab = "antidepressant initiation: primary specification")
  hr(fit(paste(RHS, MOOD), ad), lab = "antidepressant initiation: + mood items")
  hr(fit(paste(RHS, MOOD), ad[ad$mood_years >= 1, ]), lab = "antidepressant initiation: + mood items rated 1+ years before")
} else cat("antidep_dat.rds not present; run R/01_build/build_antidep.R first if this is wanted\n")
