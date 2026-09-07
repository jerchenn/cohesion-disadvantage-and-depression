# cox_absrisk3.R -- referee concern 10. Standardized (g-formula) 2-year cumulative incidence of depression
# at low vs high cohesion (-1 vs +1 SD), overall and within area-deprivation tertiles, with the absolute
# risk difference. Fits a Cox model, then averages model-implied risk over the covariate distribution under
# each counterfactual cohesion value (marginal standardization). Expresses the association, and its
# attenuation under disadvantage, on the absolute scale. Aggregate only. Needs: survival.

library(survival)
d <- readRDS("cox_dat3.rds")
d <- d[!is.na(d$deprivation_index), ]
d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                   ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
d$dep_t <- cut(d$deprivation_index, quantile(d$deprivation_index,0:3/3), include.lowest=TRUE,
               labels=c("Low disadvantage","Middle disadvantage","High disadvantage"))

absrisk <- function(sub, t=730){
  f <- coxph(Surv(time_days,event) ~ z_cohesion + age + sex_c + race_c + ethn_c +
               income_m + educ_m + log_util, sub)
  beta <- coef(f)[["z_cohesion"]]
  lp <- predict(f, type="lp")
  bh <- basehaz(f, centered=TRUE); H0t <- bh$hazard[max(which(bh$time <= t))]
  ci <- function(c) 1 - mean(exp(-H0t * exp(lp + beta*(c - sub$z_cohesion))))
  lo <- ci(-1); hi <- ci(1)
  c(HR=exp(beta), CI_low_coh=100*lo, CI_high_coh=100*hi, ARD_pct=100*(lo-hi))
}

cat("=== Standardized 2-year cumulative incidence (%), low (-1SD) vs high (+1SD) cohesion ===\n")
cat(sprintf("%-22s  HR    risk@lowCoh  risk@highCoh  abs.diff\n", "stratum"))
r <- absrisk(d)
cat(sprintf("%-22s  %.3f   %.2f%%        %.2f%%        %.2f pts\n","OVERALL", r["HR"],r["CI_low_coh"],r["CI_high_coh"],r["ARD_pct"]))
for (dt in levels(d$dep_t)){ s <- absrisk(d[d$dep_t==dt,])
  cat(sprintf("%-22s  %.3f   %.2f%%        %.2f%%        %.2f pts\n", dt, s["HR"],s["CI_low_coh"],s["CI_high_coh"],s["ARD_pct"])) }
cat("\nMarginal standardization over each stratum's covariate distribution; horizon 2 y.\n")
cat("abs.diff = excess 2-y incidence at low vs high cohesion (percentage points).\n")
