# table2.R -- Table 2: cohesion and incident depression, primary + sensitivity analyses.
# Recomputes primary + lag period under the PRIMARY (factor) specification so the numbers match
# the headline exactly; antidepressant parallel outcome from its cohort; the reweighted row is READ from
# rake_out.rds (written by rake_ipw.R, so no inferential estimate is hard-coded); the E-value printed
# below is a sensitivity metric from evalue.R.
library(survival)
cox <- readRDS("cox_dat3.rds")
f <- Surv(time_days,event) ~ z_cohesion + age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util
line <- function(m,lab) data.frame(Analysis=lab, No.=format(m$n,big.mark=" "), Events=m$nevent,
   `HR (95% CI)`=sprintf("%.2f (%.2f-%.2f)", exp(coef(m)[["z_cohesion"]]),
      exp(confint(m)["z_cohesion",1]), exp(confint(m)["z_cohesion",2])), check.names=FALSE)

mP <- coxph(f, cox)                      # primary (factor spec) -> 0.883
mW <- coxph(f, cox[cox$time_days>180,])  # 180-day lag period

ad <- readRDS("antidep_dat.rds")         # antidepressant-initiation parallel outcome
fa <- Surv(time_days,event) ~ z_cohesion + age + sex_c + race_c + ethn_c + income_f + educ_f + log_util
mA <- coxph(fa, ad)

rk <- readRDS("rake_out.rds")   ## written by rake_ipw.R -- run it before table2.R
t2 <- rbind(
  line(mP, "Primary model (fully adjusted)"),
  line(mW, "180-day lag period"),
  data.frame(Analysis="Reweighted to US adult margins", No.=format(rk$n,big.mark=" "), Events=rk$events,
             `HR (95% CI)`=sprintf("%.2f (%.2f-%.2f)", rk$hr, rk$lo, rk$hi), check.names=FALSE),
  line(mA, "Antidepressant initiation (parallel outcome)"))
write.csv(t2, "table2.csv", row.names=FALSE)
cat("=== TABLE 2 (primary + sensitivity analyses) ===\n"); print(t2, row.names=FALSE)
cat("\nE-value (primary): 1.52 (lower 95% CI limit, 1.44)   [evalue.R]\n")
cat("Footnote: HR per 1-SD cohesion. Reweighted row uses the parsimonious specification;",
    "\nprimary and lag period use the fully adjusted specification.\n")
