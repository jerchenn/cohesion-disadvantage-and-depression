# R/06_subjective_objective/07_models.R
# THE credibility test for B. Reads srh_analytic.rds (person-level, in Workbench). Nested regressions on
# ONE common complete-case sample: does the SES gradient in self-rated PHYSICAL health SURVIVE progressively
# richer OBJECTIVE-health adjustment (BMI -> +BP -> +labs)? Objective health modelled FLEXIBLY (poly, 2)
# so we don't under-adjust and leave residual health confounding.
#  - SES coef stable across models  => reporting heterogeneity (the thesis)
#  - SES coef collapses on adjustment => differential unmeasured health (null on the interesting claim)
# Labs (HbA1c/glucose) are UNOBSERVABLE to the person: SES persisting net of labs = not "they know they're sicker".
# Only aggregate coefficients leave the Workbench (DUCC).

d <- readRDS("srh_analytic.rds")
d$race <- relevel(factor(ifelse(is.na(d$race),"Other",d$race)), ref="White")

## common complete-case sample (BMI+BP+glucose; HbA1c DROPPED -- selectively measured, shrinks+biases)
cc <- complete.cases(d[,c("sr_phys","income","educ","race","age","female","bmi","sbp","glucose")])
b <- d[cc,]
cat("common complete-case n =", nrow(b), "\n\n")

ses <- function(m, lab){
  s <- summary(m)$coefficients
  rows <- grep("income|educ|raceBlack|raceAIAN|raceAsian", rownames(s), value=TRUE)
  cat("--", lab, " (R2=", round(summary(m)$r.squared,3), ") --\n", sep="")
  for (r in rows) cat(sprintf("   %-14s b=%+.3f (SE %.3f)\n", r, s[r,1], s[r,2]))
}

m0 <- lm(sr_phys ~ income + educ + race + age + female, b)
m1 <- lm(sr_phys ~ income + educ + race + age + female + poly(bmi,2), b)
m2 <- lm(sr_phys ~ income + educ + race + age + female + poly(bmi,2) + poly(sbp,2), b)
m3 <- lm(sr_phys ~ income + educ + race + age + female + poly(bmi,2) + poly(sbp,2) + poly(glucose,2), b)

cat("== SES gradient in self-rated physical health, across nested objective-health adjustment ==\n")
cat("   (income 1-3, educ 1-4 as linear per-level; race vs White; SAME n =", nrow(b), "throughout)\n\n")
ses(m0, "M0: SES + age + sex only")
ses(m1, "M1: + BMI (flexible)")
ses(m2, "M2: + BMI + systolic BP")
ses(m3, "M3: + BMI + BP + glucose (full objective)")

cat("\n== attenuation of the SES slopes from M0 -> M3 ==\n")
inc0<-coef(m0)["income"]; inc3<-coef(m3)["income"]; ed0<-coef(m0)["educ"]; ed3<-coef(m3)["educ"]
cat(sprintf("   income slope: %.3f -> %.3f  (%.0f%% retained)\n", inc0, inc3, 100*inc3/inc0))
cat(sprintf("   educ   slope: %.3f -> %.3f  (%.0f%% retained)\n", ed0, ed3, 100*ed3/ed0))
cat("   (high retention => reporting heterogeneity survives rich objective adjustment)\n")

cat("\nDONE\n")
