# r2_ascertainment.R -- Reviewer 2, third point: does health-care contact differ by cohesion, income and
# deprivation; a more specific engagement measure (outpatient visits); negative control outcomes.
#  (a) post-baseline contact: condition records and outpatient visits per person-year of follow-up, and the
#      share with any contact in the first year, by cohesion tertile, income group, deprivation tertile
#  (b) pre-baseline outpatient visit count as the utilisation covariate (replacing / alongside log_util)
#  (c) post-baseline outpatient visit rate as an additional covariate (post-exposure; sensitivity only)
#  (d) negative control outcomes: incident cataract and incident inguinal hernia, same model as the primary
# Primary specification unchanged. Aggregate output only. Needs: survival, bigrquery.
source("R/07_revision/_common.R", local = TRUE); bq()
d <- cox[!is.na(cox$z_cohesion), ]
d$end_date <- d$expo_date + d$time_days

## ---- contact counts ---------------------------------------------------------------------------------
expo_sql <- sprintf("(SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id = 40192463 GROUP BY person_id)", cdr)
cnd <- run_sql(sprintf("
SELECT co.person_id,
       COUNTIF(co.condition_start_date >= e.expo_date) AS cond_post,
       COUNTIF(co.condition_start_date >= e.expo_date AND co.condition_start_date < DATE_ADD(e.expo_date, INTERVAL 365 DAY)) AS cond_y1
FROM `%s.condition_occurrence` co JOIN %s e USING(person_id) GROUP BY co.person_id", cdr, expo_sql))
vis <- run_sql(sprintf("
SELECT v.person_id,
       COUNTIF(v.visit_start_date <  e.expo_date) AS visit_pre,
       COUNTIF(v.visit_start_date >= e.expo_date) AS visit_post,
       COUNTIF(v.visit_start_date >= e.expo_date AND v.visit_start_date < DATE_ADD(e.expo_date, INTERVAL 365 DAY)) AS visit_y1
FROM `%s.visit_occurrence` v JOIN %s e USING(person_id)
WHERE v.visit_concept_id IN (9202, 581477)      -- Outpatient Visit, Office Visit
GROUP BY v.person_id", cdr, expo_sql))
for (x in list(cnd, vis)) x$person_id <- as.character(x$person_id)
cnd$person_id <- as.character(cnd$person_id); vis$person_id <- as.character(vis$person_id)
d <- merge(d, cnd, by = "person_id", all.x = TRUE); d <- merge(d, vis, by = "person_id", all.x = TRUE)
for (v in c("cond_post","cond_y1","visit_pre","visit_post","visit_y1")) d[[v]] <- ifelse(is.na(d[[v]]), 0, as.numeric(d[[v]]))
d$fu_y <- d$time_days / 365.25
d$cond_rate <- d$cond_post / d$fu_y; d$visit_rate <- d$visit_post / d$fu_y
cat(sprintf("participants with any outpatient/office visit before baseline: %.1f%%; after baseline: %.1f%%\n",
            100*mean(d$visit_pre > 0), 100*mean(d$visit_post > 0)))

## ---- (a) contact by cohesion, income, deprivation ----------------------------------------------------
hdr("(a) Post-baseline health-care contact by group (records per person-year; share with any contact in year 1 among those followed >= 1 year or with an event)")
d$coh_t <- cut(d$z_cohesion, quantile(d$z_cohesion, 0:3/3), include.lowest = TRUE, labels = c("Low","Middle","High"))
d$inc_g <- ifelse(is.na(d$income_n), "Not reported", ifelse(d$income_n <= 3, "<$35k", ifelse(d$income_n <= 5, "$35k-$75k", ">=$75k")))
d$dep_t <- cut(d$deprivation_index, quantile(d$deprivation_index, 0:3/3, na.rm = TRUE), include.lowest = TRUE, labels = c("Low","Middle","High"))
y1 <- d$time_days >= 365 | d$event == 1
grp <- function(g, lab){
  cat(sprintf("\n-- by %s --\n", lab))
  out <- do.call(rbind, lapply(split(seq_len(nrow(d)), d[[g]]), function(i) data.frame(
    n = mask(length(i)), cond_per_y_median = sprintf("%.1f", median(d$cond_rate[i])), cond_per_y_mean = sprintf("%.1f", mean(d$cond_rate[i])),
    visits_per_y_median = sprintf("%.1f", median(d$visit_rate[i])), visits_per_y_mean = sprintf("%.1f", mean(d$visit_rate[i])),
    any_cond_y1_pct = sprintf("%.1f", 100*mean(d$cond_y1[i & y1] > 0)), any_visit_y1_pct = sprintf("%.1f", 100*mean(d$visit_y1[i & y1] > 0)))))
  print(out)
}
grp("coh_t", "cohesion tertile"); grp("inc_g", "income group"); grp("dep_t", "area deprivation tertile")
cat("\nSpearman correlation of cohesion with post-baseline visit rate:", sprintf("%.3f", cor(d$z_cohesion, d$visit_rate, method = "spearman")), "\n")

## ---- (b) pre-baseline outpatient visits as the utilisation covariate --------------------------------
hdr("(b) Utilisation covariate: outpatient visits before baseline instead of condition records")
d$log_visit_pre <- log1p(d$visit_pre)
hr(fit(RHS, d), lab = "primary (log condition records before baseline)")
hr(fit(sub("log_util", "log_visit_pre", RHS), d), lab = "log outpatient visits before baseline instead")
hr(fit(paste(RHS, "+ log_visit_pre"), d), lab = "both")
hr(fit(paste(RHS, "+ log_visit_pre"), d), "log_visit_pre", "  (HR per log visit, for the record)")

## ---- (c) post-baseline visit rate (post-exposure; sensitivity only) ----------------------------------
hdr("(c) Adding the post-baseline outpatient visit rate (conditions on a post-exposure variable; sensitivity only)")
d$log_visit_rate <- log1p(d$visit_rate)
hr(fit(paste(RHS, "+ log_visit_rate"), d), lab = "primary + log post-baseline visits per year")

## ---- (d) negative control outcomes ----------------------------------------------------------------
hdr("(d) Negative control outcomes (same covariates; HR per SD cohesion expected near 1)")
for (nc in list(c(name = "cataract", parent = "Cataract"), c(name = "inguinal hernia", parent = "Inguinal hernia"),
                c(name = "fracture of bone", parent = "Fracture of bone"))){
  k <- n_concepts(nc[["parent"]]); fd <- first_date(nc[["parent"]])
  m <- match(d$person_id, fd$person_id); first <- fd$first_dt[m]
  dn <- d[is.na(first) | first > d$expo_date, ]                         # free of the condition at baseline
  f2 <- first[is.na(first) | first > d$expo_date]
  dn$event <- as.integer(!is.na(f2) & f2 <= dn$end_date)
  dn$time_days <- ifelse(dn$event == 1, as.numeric(f2 - dn$expo_date), dn$time_days)
  dn <- dn[dn$time_days > 0, ]
  cat(sprintf("%s: %s standard concepts; at risk %s; incident %s\n", nc[["name"]], mask(k), mask(nrow(dn)), mask(sum(dn$event))))
  hr(fit(RHS, dn), lab = sprintf("  cohesion -> incident %s", nc[["name"]]))
}
cat("\nNote for the text: if contact rates and the negative controls do not vary with cohesion, differential detection is an unlikely explanation of the main association; if contact does vary, report the direction.\n")
