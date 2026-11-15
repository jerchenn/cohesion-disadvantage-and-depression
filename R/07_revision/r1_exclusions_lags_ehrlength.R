# r1_exclusions_lags_ehrlength.R -- Reviewer 1 (psychiatric exclusions, EHR length) and Reviewer 2 (longer lags).
#  (a) exclude participants with a pre-baseline record of bipolar disorder, psychotic disorder, or a
#      substance use disorder (each alone and together); anxiety shown separately because it is common
#      and overlaps with depression in care pathways
#  (b) add length of prior EHR observation (prior_ehr_days) as a covariate
#  (c) lag periods of 365 and 730 days (events and follow-up before the lag excluded, as for the 180-day lag)
# Primary specification unchanged. Aggregate output only. Needs: survival, bigrquery.
source("R/07_revision/_common.R", local = TRUE); bq()

hdr("Reference: primary model")
d0 <- cox; m0 <- fit(RHS, d0); hr(m0, lab = "primary (Table 2)")

hdr("(a) Pre-baseline psychiatric history: concept sets and prevalence")
groups <- list(
  bipolar   = c("Bipolar disorder"),
  psychotic = c("Schizophrenia", "Psychotic disorder"),
  substance = c("Psychoactive substance use disorder", "Alcohol dependence", "Drug dependence", "Alcohol abuse", "Drug abuse", "Substance abuse"),
  anxiety   = c("Anxiety disorder"))
for (g in names(groups)){
  k  <- n_concepts(groups[[g]]); fd <- first_date(groups[[g]])
  d0[[paste0("pre_", g)]] <- pre_baseline(fd, d0)
  cat(sprintf("  %-10s parents: %-110s standard descendants: %s | cohort members with a record on/before baseline: %s (%.1f%%), events among them: %s\n",
              g, paste(groups[[g]], collapse = "; "), mask(k), mask(sum(d0[[paste0("pre_", g)]])), 100*mean(d0[[paste0("pre_", g)]]),
              mask(sum(d0$event[d0[[paste0("pre_", g)]]]))))
  if (k == 0) cat("     !! no concepts matched these parent names; check spelling against the concept table\n")
}
d0$pre_smi <- d0$pre_bipolar | d0$pre_psychotic | d0$pre_substance
cat(sprintf("  any of bipolar/psychotic/substance: %s (%.1f%%)\n", mask(sum(d0$pre_smi)), 100*mean(d0$pre_smi)))

hdr("(a) Primary model after exclusions")
hr(fit(RHS, d0[!d0$pre_bipolar, ]),   lab = "excluding prior bipolar disorder")
hr(fit(RHS, d0[!d0$pre_psychotic, ]), lab = "excluding prior psychotic disorder")
hr(fit(RHS, d0[!d0$pre_substance, ]), lab = "excluding prior substance use disorder")
hr(fit(RHS, d0[!d0$pre_smi, ]),       lab = "excluding all three (bipolar, psychotic, substance)")
hr(fit(RHS, d0[!d0$pre_smi & !d0$pre_anxiety, ]), lab = "excluding all three plus prior anxiety disorder")
hr(fit(paste(RHS, "+ pre_smi + pre_anxiety"), d0), lab = "full cohort, adjusting for the two history indicators")

hdr("(b) Length of prior EHR observation as a covariate")
cox$prior_ehr_y <- cox$prior_ehr_days / 365.25
cat(sprintf("prior EHR observation, years: median %.1f, IQR %.1f-%.1f\n", median(cox$prior_ehr_y), quantile(cox$prior_ehr_y, .25), quantile(cox$prior_ehr_y, .75)))
hr(fit(paste(RHS, "+ prior_ehr_y"), cox), lab = "primary + prior EHR years (linear)")
m_e <- fit(paste(RHS, "+ prior_ehr_y"), cox); hr(m_e, "prior_ehr_y", "  (HR per additional year of prior observation)")
hr(fit(paste(RHS, "+ log(prior_ehr_y)"), cox), lab = "primary + log prior EHR years")
hr(fit(paste(RHS, "+ prior_ehr_y + I(prior_ehr_y^2)"), cox), lab = "primary + prior EHR years (quadratic)")

hdr("(c) Longer lag periods (same construction as the 180-day lag in Table 2)")
for (L in c(180, 365, 730)){
  dl <- cox[cox$time_days > L, ]
  hr(fit(RHS, dl), lab = sprintf("lag %d days (events in first %d days excluded)", L, L))
}
cat("\nNote for the text: the lag analyses condition on remaining in follow-up past the lag; the time origin is unchanged.\n")
