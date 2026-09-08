# cox_absrisk3.R -- absolute risk (eTable 7). Marginal (model-based) standardization of 2-year cumulative
# incidence of depression at low vs high cohesion (-1 vs +1 SD), overall and within area-deprivation
# tertiles, with the absolute risk difference and a nonparametric bootstrap 95% CI for that difference.
# Model-based standardization from the fitted Cox model, not an explicit counterfactual g-formula.
# Baseline cumulative hazard is evaluated explicitly at t = 730 days. Aggregate only. Needs: survival.

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

## standardized 2-year risk at cohesion = -1 and +1 SD, and their difference (single fit)
ard_point <- function(sub, t=730){
  f  <- coxph(Surv(time_days,event) ~ z_cohesion + age + sex_c + race_c + ethn_c +
                income_m + educ_m + log_util, sub)
  beta <- coef(f)[["z_cohesion"]]; lp <- predict(f, type="lp")
  bh <- basehaz(f, centered=TRUE)
  H0t <- approx(bh$time, bh$hazard, xout=t, method="constant", rule=2)$y   # explicit H0(730)
  risk <- function(cc) 1 - mean(exp(-H0t * exp(lp + beta*(cc - sub$z_cohesion))))
  lo <- risk(-1); hi <- risk(1)
  c(hr=exp(beta), lo=lo, hi=hi, ard=lo-hi)
}

## point estimate + nonparametric bootstrap CI for the absolute risk difference (slow: B fits per stratum)
absrisk <- function(sub, t=730, B=200){
  vars <- c("z_cohesion","age","sex_c","race_c","ethn_c","income_m","educ_m","log_util","time_days","event")
  sub <- sub[complete.cases(sub[,vars]), ]
  pt <- ard_point(sub, t)
  set.seed(1)
  bard <- replicate(B, { i <- sample.int(nrow(sub), nrow(sub), TRUE)
    tryCatch(ard_point(sub[i,], t)[["ard"]], error=function(e) NA_real_) })
  ci <- quantile(bard, c(.025,.975), na.rm=TRUE)
  c(HR=pt[["hr"]], low=100*pt[["lo"]], high=100*pt[["hi"]], ARD=100*pt[["ard"]],
    ARD_lo=100*ci[[1]], ARD_hi=100*ci[[2]])
}

cat("=== Standardized 2-year cumulative incidence (%), low (-1SD) vs high (+1SD) cohesion ===\n")
cat(sprintf("%-22s  HR     risk@low  risk@high  abs.diff (95%% CI)\n", "stratum"))
prnt <- function(lab,x) cat(sprintf("%-22s  %.3f  %.2f%%    %.2f%%    %.2f pts (%.2f to %.2f)\n",
  lab, x["HR"], x["low"], x["high"], x["ARD"], x["ARD_lo"], x["ARD_hi"]))
prnt("OVERALL", absrisk(d))
for (dt in levels(d$dep_t)) prnt(dt, absrisk(d[d$dep_t==dt,]))
cat("\nModel-based standardization over each stratum's covariate distribution; horizon 2 y;\n")
cat("abs.diff = excess 2-y incidence at low vs high cohesion; CI is a 200-rep nonparametric bootstrap.\n")
