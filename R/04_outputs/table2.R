# table2.R -- Table 2: cohesion and incident depression, primary + sensitivity analyses.
# Recomputes primary + lag period under the PRIMARY (factor) specification so the numbers match
# the headline exactly; antidepressant parallel outcome from its cohort; raked + E-value are
# stable single values carried from rake_ipw.R / evalue.R (commented provenance).
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

t2 <- rbind(
  line(mP, "Primary model (fully adjusted)"),
  line(mW, "180-day lag period"),
  data.frame(Analysis="Reweighted to US adult margins", No.="29 836", Events=2080,
             `HR (95% CI)`="0.88 (0.83-0.95)", check.names=FALSE),          # rake_ipw.R (bootstrap CI)
  line(mA, "Antidepressant initiation (parallel outcome)"))
write.csv(t2, "table2.csv", row.names=FALSE)
cat("=== TABLE 2 (primary + sensitivity analyses) ===\n"); print(t2, row.names=FALSE)
cat("\nE-value (primary): 1.52 (lower 95% CI limit, 1.44)   [evalue.R]\n")
cat("Footnote: HR per 1-SD cohesion. Reweighted row uses the parsimonious specification;",
    "\nprimary and lag period use the fully adjusted specification.\n")
