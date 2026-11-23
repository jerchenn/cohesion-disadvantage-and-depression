# r5_movers.R -- does the attenuation by baseline mood (r2_baseline_mood.R: 0.88 -> 0.97) look like confounding or mediation?
# Test: among participants who moved home in the 12 months before the cohesion survey, a mood rating given before that
# year predates the current neighborhood and cannot mediate its cohesion. If the attenuation is as large among recent
# movers as among long-term residents, prior mood is acting as a confounder. The SDOH survey asks how many times the
# respondent moved in the last 12 months; the question label is discovered from the CDR. Aggregate output only.
source("R/07_revision/_common.R", local = TRUE); bq()

hdr("(a) Moves item: label and answers")
labs <- run_sql(sprintf("SELECT question, COUNT(DISTINCT person_id) AS n FROM `%s.ds_survey` WHERE LOWER(question) LIKE '%%moved%%' OR LOWER(question) LIKE '%%one home to another%%' GROUP BY question ORDER BY n DESC", cdr))
for (i in seq_len(nrow(labs))) cat(sprintf("  %-90s persons %s\n", labs$question[i], mask(labs$n[i])))
q_mv <- labs$question[grepl("moved", labs$question, ignore.case = TRUE)][1]; cat("Item used:", q_mv, "\n")
mv <- run_sql(sprintf("SELECT person_id, LOWER(answer) AS ans, DATE(survey_datetime) AS dt FROM `%s.ds_survey` WHERE question = '%s'", cdr, q_mv))
mv$person_id <- as.character(mv$person_id); mv$dt <- as.Date(mv$dt)
tb <- sort(table(mv$ans), decreasing = TRUE); for (k in names(tb)) cat(sprintf("  %-50s %s\n", k, mask(tb[[k]])))
mv$moved <- ifelse(grepl("skip|prefer|know", mv$ans), NA, ifelse(grepl("^[^0-9]*0|none|never|not moved|did not", mv$ans), 0L, 1L))
cat(sprintf("mapping: moved = 0 for answers reading 0/none/never, 1 for any other numeric answer; unmapped %s rows\n", mask(sum(is.na(mv$moved)))))

hdr("(b) Mood items on or before baseline (as in r2_baseline_mood.R)")
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
## moves answer closest to (at or before) the cohesion baseline
m2 <- mv[!is.na(mv$moved) & mv$person_id %in% d$person_id, ]; m2$expo <- d$expo_date[match(m2$person_id, d$person_id)]; m2 <- m2[m2$dt <= m2$expo + 30, ]
m2 <- m2[order(m2$person_id, abs(as.numeric(m2$dt - m2$expo))), ]; m2 <- m2[!duplicated(m2$person_id), ]
d$moved <- m2$moved[match(d$person_id, m2$person_id)]
d$mood_years_before <- as.numeric(d$expo_date - d$sr_mental_dt) / 365.25
db <- d[!is.na(d$sr_mental) & !is.na(d$emo_prob) & !is.na(d$moved), ]
db$sr_mental_f <- relevel(factor(db$sr_mental), ref = "5"); db$emo_prob_f <- relevel(factor(db$emo_prob), ref = "0")
cat(sprintf("cohort members with both mood items and a moves answer: %s; moved in the past 12 months: %s (%.1f%%)\n", mask(nrow(db)), mask(sum(db$moved)), 100*mean(db$moved)))

hdr("(c) Attenuation by baseline mood, among recent movers and long-term residents")
att <- function(dd, lab){
  m0 <- fit(RHS, dd); m1 <- fit(paste(RHS, "+ sr_mental_f + emo_prob_f"), dd)
  hr(m0, lab = paste(lab, ": primary")); hr(m1, lab = paste(lab, ": + both mood items"))
  cat(sprintf("    share of the log HR removed by adjustment: %.0f%%\n", 100 * (1 - coef(m1)[["z_cohesion"]] / coef(m0)[["z_cohesion"]]))) }
att(db[db$moved == 0, ], "did not move in the past 12 months")
att(db[db$moved == 1, ], "moved in the past 12 months")
att(db[db$moved == 1 & db$mood_years_before >= 1, ], "moved in the past 12 months, mood rated 1+ years before baseline (mood predates the current neighborhood)")
att(db[db$moved == 0 & db$mood_years_before >= 1, ], "did not move, mood rated 1+ years before baseline")
cat("\nReading: under mediation, adjustment should remove less of the association among recent movers whose mood predates the move; under confounding, the same share.\n")
