# r2_ipcw.R -- Reviewer 2, third point: informative censoring. Follow-up ends at the last EHR record, which
# may depend on cohesion, income or deprivation. Inverse probability of censoring weights:
#  - "lost" = no event and follow-up ended more than 180 days before the end of the data (administrative end)
#  - censoring model: Cox for time to loss with the primary covariates (deprivation -> loss is a separate descriptive model)
#  - stabilised weight = S0(t_i) / S(t_i | X_i), evaluated at each person's own end of follow-up, truncated at the 99th pct
#  - primary refitted with the weights and robust SE
# Primary specification unchanged. Aggregate output only. Needs: survival.
source("R/07_revision/_common.R", local = TRUE)
d <- cox[!is.na(cox$z_cohesion), ]
d$end_date <- d$expo_date + d$time_days
admin_end <- max(d$end_date)                                    # used internally only; not printed
for (gap in c(180, 90)){
  hdr(sprintf("IPCW with loss defined as follow-up ending more than %d days before the end of the data", gap))
  d$lost <- as.integer(d$event == 0 & d$end_date < admin_end - gap)
  cat(sprintf("lost to follow-up: %s (%.1f%%); administratively censored: %s; events: %s\n",
              mask(sum(d$lost)), 100*mean(d$lost), mask(sum(d$event == 0 & d$lost == 0)), mask(sum(d$event))))
  ## censoring model
  cm <- coxph(as.formula(paste("Surv(time_days, lost) ~", RHS)), d)
  hr(cm, "z_cohesion", "censoring model: cohesion -> loss to follow-up (HR per SD)")
  hr(cm, "income_fNA", "censoring model: income not reported -> loss")
  cm0 <- coxph(Surv(time_days, lost) ~ 1, d)
  H  <- basehaz(cm,  centered = TRUE); H0 <- basehaz(cm0, centered = TRUE)
  Hi  <- approx(H$time,  H$hazard,  xout = d$time_days, method = "constant", rule = 2, f = 0)$y
  H0i <- approx(H0$time, H0$hazard, xout = d$time_days, method = "constant", rule = 2, f = 0)$y
  lp <- predict(cm, type = "lp")                                   # centred at sample means, matching centered=TRUE
  S_i <- exp(-Hi * exp(lp)); S_0 <- exp(-H0i)
  w <- S_0 / S_i
  cap <- quantile(w, .99); w[w > cap] <- cap
  d$w_ipcw <- w
  cat(sprintf("stabilised weights: median %.3f, IQR %.3f-%.3f, 99th pct %.2f (truncation point)\n", median(w), quantile(w,.25), quantile(w,.75), cap))
  hr(fit(RHS, d), lab = "primary, unweighted")
  mw <- coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), d, weights = w_ipcw, robust = TRUE)
  hr(mw, lab = "primary, IPCW-weighted (robust SE)")
  ## with deprivation in the censoring model
  dd <- d[!is.na(d$deprivation_index), ]; dd$dep_z <- as.numeric(scale(dd$deprivation_index))
  cmd <- coxph(as.formula(paste("Surv(time_days, lost) ~", RHS, "+ dep_z")), dd)
  hr(cmd, "dep_z", "censoring model: area deprivation -> loss (per SD)")
}
cat("\nNote for the text: weights near 1 and an unchanged HR mean the end of follow-up is not informative with respect to cohesion given the covariates.\n")
