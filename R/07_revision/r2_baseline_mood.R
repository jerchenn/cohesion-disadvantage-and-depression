# r2_baseline_mood.R -- Reviewer 2, item 2: adjust for markers of baseline mood collected BEFORE the cohesion survey.
# The Overall Health survey (offered at enrolment, before the SDOH survey for most participants) carries the two
# PROMIS Global Health mental items: self-rated mental health (excellent to poor) and frequency of being bothered by
# emotional problems such as feeling anxious, depressed or irritable in the past 7 days (never to always).
# This script (a) lists the Overall Health question labels in the CDR, (b) locks the answer scales, (c) keeps each
# cohort member's latest answer dated on or before the cohesion baseline (expo_date), (d) reports coverage and the
# interval between the two surveys, (e) refits the primary model with the items, alone and with the prior anxiety and
# psychotropic markers of baseline_history.R. Primary specification unchanged. Aggregate output only.
# The letter sentence "All of Us has no measure of mood or personality collected before the cohesion survey" is wrong
# and must be replaced with the result printed at the end.
source("R/07_revision/_common.R", local = TRUE); bq()

hdr("(a) Overall Health question labels in this CDR (labels only; persons masked below 20)")
labs <- run_sql(sprintf("SELECT question, COUNT(DISTINCT person_id) AS n FROM `%s.ds_survey` WHERE question LIKE 'Overall Health:%%' GROUP BY question ORDER BY question", cdr))
for (i in seq_len(nrow(labs))) cat(sprintf("  %-70s persons %s\n", labs$question[i], mask(labs$n[i])))
q_mh <- "Overall Health: General Mental Health"
q_em <- labs$question[grepl("emotional", labs$question, ignore.case = TRUE)]
if (length(q_em) != 1) { cat("\n!! emotional-problems item not identified uniquely; candidates:", paste(q_em, collapse = " | "), "\n"); q_em <- q_em[1] }
cat(sprintf("\nItems used: [%s] and [%s]\n", q_mh, q_em))

hdr("(b) Raw answer distributions (lock the scales before trusting the mapping)")
sv <- run_sql(sprintf("SELECT person_id, question, LOWER(answer) AS ans, DATE(survey_datetime) AS dt FROM `%s.ds_survey` WHERE question IN ('%s','%s')", cdr, q_mh, q_em))
sv$person_id <- as.character(sv$person_id); sv$dt <- as.Date(sv$dt)
for (qq in c(q_mh, q_em)) { cat("\n--", qq, "--\n"); tb <- sort(table(sv$ans[sv$question == qq]), decreasing = TRUE); for (k in names(tb)) cat(sprintf("  %-45s %s\n", k, mask(tb[[k]]))) }
map_mh <- function(a) ifelse(grepl("excellent|excllent", a), 5L, ifelse(grepl("very good", a), 4L, ifelse(grepl("\\bgood\\b", a), 3L, ifelse(grepl("fair", a), 2L, ifelse(grepl("poor", a), 1L, NA_integer_)))))   # higher = better
map_em <- function(a) ifelse(grepl("never", a), 0L, ifelse(grepl("rarely", a), 1L, ifelse(grepl("sometimes", a), 2L, ifelse(grepl("often", a), 3L, ifelse(grepl("always", a), 4L, NA_integer_)))))   # higher = more bothered
sv$v <- ifelse(sv$question == q_mh, map_mh(sv$ans), map_em(sv$ans))
cat(sprintf("\nunmapped answers (skip, prefer not to answer, other): mental health %s of %s rows; emotional problems %s of %s rows\n",
            mask(sum(is.na(sv$v) & sv$question == q_mh)), mask(sum(sv$question == q_mh)), mask(sum(is.na(sv$v) & sv$question == q_em)), mask(sum(sv$question == q_em))))
sv <- sv[!is.na(sv$v), ]

hdr("(c) Latest answer on or before the cohesion baseline, cohort members")
d <- cox
pick <- function(qq, nm){
  x <- sv[sv$question == qq, ]; m <- match(x$person_id, d$person_id); x <- x[!is.na(m), ]; x$expo <- d$expo_date[match(x$person_id, d$person_id)]
  x <- x[x$dt <= x$expo, ]; x <- x[order(x$person_id, x$dt, decreasing = TRUE), ]; x <- x[!duplicated(x$person_id), ]
  d[[nm]] <<- x$v[match(d$person_id, x$person_id)]; d[[paste0(nm, "_dt")]] <<- x$dt[match(d$person_id, x$person_id)]
  after <- sv[sv$question == qq & sv$person_id %in% d$person_id, ]; after$expo <- d$expo_date[match(after$person_id, d$person_id)]; after <- after[after$dt > after$expo, ]
  cat(sprintf("  %-45s on/before baseline: %s (%.1f%%); strictly before: %s; same day: %s; only after baseline (not used): %s persons\n", nm,
              mask(sum(!is.na(d[[nm]]))), 100*mean(!is.na(d[[nm]])), mask(sum(!is.na(d[[nm]]) & d[[paste0(nm, "_dt")]] < d$expo_date)),
              mask(sum(!is.na(d[[nm]]) & d[[paste0(nm, "_dt")]] == d$expo_date)), mask(length(setdiff(unique(after$person_id), d$person_id[!is.na(d[[nm]])])))))
}
pick(q_mh, "sr_mental"); pick(q_em, "emo_prob")
d$has_both <- !is.na(d$sr_mental) & !is.na(d$emo_prob)
gap <- as.numeric(d$expo_date - d$sr_mental_dt) / 365.25
cat(sprintf("  both items: %s (%.1f%%). Years from Overall Health answer to cohesion baseline: median %.1f, IQR %.1f-%.1f; share answered 1+ years before: %.1f%%\n",
            mask(sum(d$has_both)), 100*mean(d$has_both), median(gap, na.rm = TRUE), quantile(gap, .25, na.rm = TRUE), quantile(gap, .75, na.rm = TRUE), 100*mean(gap >= 1, na.rm = TRUE)))
cat("  distribution in the cohort (persons): self-rated mental health 1 poor .. 5 excellent | emotional problems 0 never .. 4 always\n")
for (k in 1:5) cat(sprintf("    mental health = %d: %s\n", k, mask(sum(d$sr_mental == k, na.rm = TRUE))))
for (k in 0:4) cat(sprintf("    emotional problems = %d: %s\n", k, mask(sum(d$emo_prob == k, na.rm = TRUE))))
cat(sprintf("  correlation of the two items in the cohort: r = %.2f; correlation with z_cohesion: mental health %.2f, emotional problems %.2f\n",
            cor(d$sr_mental, d$emo_prob, use = "complete.obs"), cor(d$sr_mental, d$z_cohesion, use = "complete.obs"), cor(d$emo_prob, d$z_cohesion, use = "complete.obs")))

WITH_PRIOR <- Sys.getenv("MOOD_PRIOR") == "1"   # the anxiety and drug pulls are the heavy step; the reviewer's request is answered without them
hdr("(d) Prior anxiety diagnosis and psychotropic prescription (as in baseline_history.R), for the combined model")
if (WITH_PRIOR) {
anx <- first_date("Anxiety disorder"); d$prior_anx <- as.integer(pre_baseline(anx, d) & anx$first_dt[match(d$person_id, anx$person_id)] < d$expo_date)
med <- run_sql(sprintf("WITH atc AS (SELECT concept_id FROM `%s.concept` WHERE vocabulary_id='ATC' AND concept_code IN ('N05A','N05B','N06A')), ids AS (SELECT DISTINCT ca.descendant_concept_id AS concept_id FROM `%s.concept_ancestor` ca JOIN atc a ON ca.ancestor_concept_id=a.concept_id) SELECT de.person_id, MIN(DATE(de.drug_exposure_start_date)) AS first_dt FROM `%s.drug_exposure` de JOIN ids i ON de.drug_concept_id=i.concept_id GROUP BY de.person_id", cdr, cdr, cdr))
med$person_id <- as.character(med$person_id); med$first_dt <- as.Date(med$first_dt)
d$prior_med <- as.integer(!is.na(match(d$person_id, med$person_id)) & med$first_dt[match(d$person_id, med$person_id)] < d$expo_date); d$prior_med[is.na(d$prior_med)] <- 0L
d$prior_anx[is.na(d$prior_anx)] <- 0L
cat(sprintf("  prior anxiety: %s (%.1f%%); prior psychotropic: %s (%.1f%%)\n", mask(sum(d$prior_anx)), 100*mean(d$prior_anx), mask(sum(d$prior_med)), 100*mean(d$prior_med)))
} else { d$prior_anx <- 0L; d$prior_med <- 0L; cat("  skipped (set MOOD_PRIOR=1 to include the EHR markers)\n") }

hdr("(e) Primary model with baseline mood items")
hr(fit(RHS, d), lab = "primary, full cohort (Table 2)")
db <- d[d$has_both, ]
hr(fit(RHS, db), lab = "primary, participants with both items on/before baseline (reference)")
db$sr_mental_f <- relevel(factor(db$sr_mental), ref = "5"); db$emo_prob_f <- relevel(factor(db$emo_prob), ref = "0")
hr(fit(paste(RHS, "+ sr_mental_f"), db),              lab = "  + self-rated mental health (categorical, ref excellent)")
hr(fit(paste(RHS, "+ emo_prob_f"), db),               lab = "  + emotional problems 7 days (categorical, ref never)")
m_both <- fit(paste(RHS, "+ sr_mental_f + emo_prob_f"), db); hr(m_both, lab = "  + both items (categorical)")
hr(fit(paste(RHS, "+ sr_mental + emo_prob"), db),     lab = "  + both items (linear)")
if (WITH_PRIOR) { m_all <- fit(paste(RHS, "+ sr_mental_f + emo_prob_f + prior_anx + prior_med"), db); hr(m_all, lab = "  + both items + prior anxiety + prior psychotropic") } else m_all <- m_both
hr(fit(paste(RHS, "+ sr_mental_f + emo_prob_f"), db[db$sr_mental_dt < db$expo_date, ]), lab = "  both items, answered strictly before baseline only")
cat("\n  Item HRs in the model with both items (validity: worse mood should predict depression):\n")
cf <- summary(m_both)$coefficients; for (r in grep("sr_mental_f|emo_prob_f", rownames(cf), value = TRUE)) cat(sprintf("    %-16s HR %.2f (%.2f-%.2f)\n", r, cf[r, "exp(coef)"], exp(cf[r, "coef"] - 1.96*cf[r, "se(coef)"]), exp(cf[r, "coef"] + 1.96*cf[r, "se(coef)"])))

hdr("(f) Interaction check: does adjustment change the income and deprivation modification?")
## continuous interactions as in r2_reri_subsamples_vif.R (dep_z from deprivation_index; income_z from income_n, income factor dropped)
dd <- db[!is.na(db$deprivation_index), ]; dd$dep_z <- as.numeric(scale(dd$deprivation_index))
hr(fit(paste("z_cohesion * dep_z +", RHS_NOCOH), dd), "z_cohesion:dep_z", "  cohesion x deprivation, without mood items")
hr(fit(paste("z_cohesion * dep_z + sr_mental_f + emo_prob_f +", RHS_NOCOH), dd), "z_cohesion:dep_z", "  cohesion x deprivation, with mood items")
di <- db[!is.na(db$income_n), ]; di$income_z <- as.numeric(scale(di$income_n))
hr(fit(paste("z_cohesion * income_z +", sub("income_f \\+ ", "", RHS_NOCOH)), di), "z_cohesion:income_z", "  cohesion x income, without mood items")
hr(fit(paste("z_cohesion * income_z + sr_mental_f + emo_prob_f +", sub("income_f \\+ ", "", RHS_NOCOH)), di), "z_cohesion:income_z", "  cohesion x income, with mood items")

hdr("Sentence for the letter and the Limitations (fill from above)")
b <- coef(m_all)[["z_cohesion"]]; se <- sqrt(vcov(m_all)["z_cohesion", "z_cohesion"])
cat(sprintf("The Overall Health survey, completed at enrolment, includes the PROMIS items on self-rated mental health and on how often the participant was bothered by emotional problems such as feeling anxious, depressed, or irritable. Among %s participants with both items recorded on or before the cohesion survey (median %.1f years earlier), adjusting for them%s gave a hazard ratio of %.2f (95%% CI, %.2f-%.2f).\n",
            mask(m_all$n), median(gap, na.rm = TRUE), if (WITH_PRIOR) " together with prior anxiety and psychotropic prescription" else "", exp(b), exp(b - 1.96*se), exp(b + 1.96*se)))
