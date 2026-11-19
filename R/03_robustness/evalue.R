# evalue.R -- E-values for the key HRs (VanderWeele-Ding). Quantifies how strong an
# unmeasured confounder (e.g. residential self-selection) would have to be, on the
# risk-ratio scale, to explain away each association. ~5.5% incidence -> rare=TRUE ok.
# install.packages("EValue") if needed
# The primary HR is computed here from the cohort file, so the E-value follows the data. E-values for the
# Figure 3 measures (steps, happiness) are computed in the supplement from the logged Figure 3 estimates.
library(EValue); library(survival)
show <- function(lab, est, lo, hi){
  e <- evalues.HR(est, lo, hi, rare = TRUE)
  ci <- e["E-values", c("lower","upper")]
  ci <- ci[!is.na(ci)][1]                     # for protective HR the CI E-value is the non-NA bound
  cat(sprintf("%-28s HR %.4f (%.4f-%.4f)  E-value point %.2f | CI %.2f\n",
              lab, est, lo, hi, e["E-values","point"], ci))
}
cat("=== E-values (point estimate | CI limit nearest null) ===\n")
cox <- readRDS("cox_dat3.rds")
m <- coxph(Surv(time_days,event) ~ z_cohesion + age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util, cox)
ci <- exp(confint(m)["z_cohesion",])
show("Primary: cohesion->dep", exp(coef(m)[["z_cohesion"]]), ci[1], ci[2])
