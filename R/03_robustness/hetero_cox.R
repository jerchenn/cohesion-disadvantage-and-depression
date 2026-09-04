# hetero_cox.R -- heterogeneity + non-linearity of cohesion -> INCIDENT DEPRESSION
# (the prospective clinical outcome, not the earlier cross-sectional wellbeing one).
# Goal: a policy-actionable shape/subgroup signal (who benefits, where the curve bends).
library(survival); library(splines)
cox <- readRDS("cox_dat.rds")
cox <- cox[!is.na(cox$z_cohesion), ]
cox$sex_c  <- factor(ifelse(cox$sex %in% c("Female","Male"), cox$sex, "Other"))
cox$race_c <- factor(ifelse(cox$race %in% c("White","Black or African American"), cox$race, "Other"))
cox$ethn_c <- factor(ifelse(cox$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(cox$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
cox$income_m <- ifelse(is.na(cox$income_n), median(cox$income_n,na.rm=TRUE), cox$income_n)
cox$educ_m   <- ifelse(is.na(cox$educ_n),   median(cox$educ_n,  na.rm=TRUE), cox$educ_n)
covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"

## ---- 1. NON-LINEARITY / THRESHOLD -------------------------------------------
# quartiles of cohesion, reference = Q1 (lowest) -> shows where risk drops most
cox$coh_q <- cut(cox$z_cohesion, quantile(cox$z_cohesion, 0:4/4, na.rm=TRUE),
                 include.lowest=TRUE, labels=c("Q1low","Q2","Q3","Q4high"))
mq <- coxph(as.formula(paste("Surv(time_days,event) ~ coh_q +", covs)), cox)
cat("=== Cohesion quartiles (ref = Q1 lowest); HR<1 => lower risk than lowest-cohesion ===\n")
for (q in c("coh_qQ2","coh_qQ3","coh_qQ4high"))
  cat(sprintf("  %-12s HR %.3f (%.3f-%.3f)\n", q, exp(coef(mq)[[q]]),
      exp(confint(mq)[q,1]), exp(confint(mq)[q,2])))
# formal curvature: quadratic term + spline LR test
mlin <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion +", covs)), cox)
mqua <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion + I(z_cohesion^2) +", covs)), cox)
msp  <- coxph(as.formula(paste("Surv(time_days,event) ~ ns(z_cohesion,3) +", covs)), cox)
cat(sprintf("\n  quadratic term: %+.4f (p=%.3g)   [<0 concave/diminishing, >0 convex/accelerating]\n",
    coef(mqua)[["I(z_cohesion^2)"]], summary(mqua)$coefficients["I(z_cohesion^2)","Pr(>|z|)"]))
cat(sprintf("  non-linearity LR test (spline vs linear): p = %.3g\n", anova(mlin, msp)$`P(>|Chi|)`[2]))
# predicted HR across cohesion deciles (ref = median) for the shape figure
dec <- quantile(cox$z_cohesion, 1:9/10, na.rm=TRUE)
b <- coef(msp); X <- ns(dec, 3) - matrix(ns(median(cox$z_cohesion),3), nrow=length(dec), ncol=3, byrow=TRUE)
cat("  predicted HR vs median across deciles:\n   ",
    paste(sprintf("%.2f", exp(X %*% b[1:3])), collapse=" "), "\n")

## ---- 2. EFFECT MODIFICATION (who benefits most) -----------------------------
# drop the SES covariate that is collinear with each modifier
covs_noinc <- "age + sex_c + race_c + ethn_c + educ_m + log_util"
covs_noedu <- "age + sex_c + race_c + ethn_c + income_m + log_util"
covs_noage <- "sex_c + race_c + ethn_c + income_m + educ_m + log_util"
imod <- function(mvar, lab, cv){
  m <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion*", mvar, "+", cv)), cox)
  rn <- grep(":", names(coef(m)), value=TRUE)[1]
  cat(sprintf("  %-22s interaction HR %.3f (p=%.3g)\n", lab, exp(coef(m)[[rn]]),
      summary(m)$coefficients[rn,"Pr(>|z|)"]))
}
cox$low_income <- as.integer(cox$income_m <= 4)      # <= $50k band
cox$no_college <- ifelse(is.na(cox$educ_n), NA, as.integer(cox$educ_n < 4))  # raw educ, keep NA
cox$age65      <- as.integer(cox$age >= 65)
cat("\n=== Effect modification (interaction with z_cohesion) ===\n")
cat("(interaction HR >1 => cohesion LESS protective in the '1' group)\n")
imod("low_income", "cohesion x low-income", covs_noinc)
imod("no_college", "cohesion x no-college", covs_noedu)
imod("age65",      "cohesion x age>=65",    covs_noage)

## subgroup-specific cohesion HRs (interpretable), dropping the collinear covariate ----
sub <- function(cond, lab, cv){
  d <- cox[which(cond), ]
  if (nrow(d) < 50 || sum(d$event) < 10){ cat(sprintf("  %-22s (too few obs/events — skipped)\n", lab)); return(invisible()) }
  m <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion +", cv)), d)
  cat(sprintf("  %-22s HR %.3f (%.3f-%.3f)  [ev %d]\n", lab, exp(coef(m)[["z_cohesion"]]),
      exp(confint(m)["z_cohesion",1]), exp(confint(m)["z_cohesion",2]), m$nevent))
}
cat("\n=== Cohesion HR within subgroups ===\n")
sub(cox$low_income==1, "low income (<=50k)", covs_noinc); sub(cox$low_income==0, "higher income", covs_noinc)
sub(cox$no_college==1, "no college", covs_noedu);         sub(cox$no_college==0, "college+", covs_noedu)
sub(cox$age65==1,      "age >= 65", covs_noage);          sub(cox$age65==0,      "age < 65", covs_noage)
