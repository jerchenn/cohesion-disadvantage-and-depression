# mediation.R -- formal mediation of cohesion -> incident depression through
#  (a) objective physical activity (steps)  and  (b) self-reported affect (wellbeing).
# Difference method on the log-HR scale (Cox), reported as EXPLORATORY ATTENUATION after adjustment,
# NOT a causal proportion mediated. HR non-collapsibility (direct/total not cleanly separable in Cox)
# and cross-sectional mediators (measured near baseline, not strictly on the exposure->outcome path)
# preclude a causal decomposition.
library(survival)

pm <- function(dat, med, lab){
  covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
  # collapsed covariates
  dat$sex_c  <- factor(ifelse(dat$sex %in% c("Female","Male"), dat$sex, "Other"))
  dat$race_c <- factor(ifelse(dat$race %in% c("White","Black or African American"), dat$race, "Other"))
  dat$ethn_c <- factor(ifelse(dat$ethnicity=="Hispanic or Latino","Hispanic",
                       ifelse(dat$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  dat$income_m <- ifelse(is.na(dat$income_n), median(dat$income_n,na.rm=TRUE), dat$income_n)
  dat$educ_m   <- ifelse(is.na(dat$educ_n),   median(dat$educ_n,  na.rm=TRUE), dat$educ_n)
  keep <- !is.na(dat$z_cohesion) & !is.na(dat[[med]])
  d <- dat[keep, ]
  a  <- lm(as.formula(paste(med, "~ z_cohesion +", covs)), d)$coef[["z_cohesion"]]  # a-path
  mt <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion +", covs)), d)
  md <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion +", med, "+", covs)), d)
  bt <- coef(mt)[["z_cohesion"]]; bd <- coef(md)[["z_cohesion"]]
  cat(sprintf("\n[%s]  N=%d events=%d\n", lab, nrow(d), mt$nevent))
  cat(sprintf("  a-path cohesion->%s : %+.3f SD\n", med, a))
  cat(sprintf("  total  HR(cohesion)        : %.3f\n", exp(bt)))
  cat(sprintf("  direct HR(cohesion|%s): %.3f\n", med, exp(bd)))
  cat(sprintf("  log-HR attenuation after adjustment: %.1f%% (exploratory; not a causal proportion mediated)\n", 100*(bt-bd)/bt))
}

fb  <- readRDS("fitbit_dat.rds")
fb$z_steps <- scale(fb$steps)[,1]                 # higher = more steps
pm(fb, "z_steps", "Objective physical activity (steps)")

swb <- readRDS("swb_dat.rds")
pm(swb, "z_wb", "Self-reported affect (wellbeing index)")

cat("\nCaveats: Cox HRs are non-collapsible (direct/total shift partly for non-causal reasons);\n")
cat("mediators are cross-sectional; treat proportions as approximate, not identified NIE/NDE.\n")
