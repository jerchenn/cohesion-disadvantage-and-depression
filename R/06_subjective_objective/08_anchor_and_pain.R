# R/06_subjective_objective/08_anchor_and_pain.R
# Two tests of the reporting-heterogeneity interpretation (advisor). Reads srh_analytic.rds (in-Workbench).
#
# TEST 1 -- EXTERNAL ANCHOR: compare the SES gradient in SUBJECTIVE health to the SES gradient in an
#   OBJECTIVE health index of comparable scale (both standardized). If the SUBJECTIVE gradient is STEEPER,
#   self-report exaggerates the social health gradient => reporting heterogeneity, without claiming to have
#   measured all objective health. (income + educ reported SEPARATELY.)
# TEST 2 -- PAIN as unmeasured symptom burden: does the income/educ gradient in self-rated PHYSICAL health
#   shrink when self-reported PAIN (the domain BMI/BP/glucose miss) is added? Large shrink => part of the
#   "reporting" gap is real symptom burden. (Caveat: pain is also self-reported -> shared reporting style
#   makes this a conservative lower bound on the SES gradient.)
# Aggregate output only.

d <- readRDS("srh_analytic.rds")
d$race <- relevel(factor(ifelse(is.na(d$race),"Other",d$race)), ref="White")
d$age  <- 2024 - as.numeric(d$year_of_birth)

z <- function(x) as.numeric(scale(x))
cc <- complete.cases(d[,c("sr_phys","income","educ","race","age","bmi","sbp","glucose")])
b <- d[cc,]

## objective health index (higher = HEALTHIER), standardized: negative sum of standardized risk markers
b$obj_health <- z(-( z(b$bmi) + z(b$sbp) + z(b$glucose) ))
b$sr_phys_z  <- z(b$sr_phys)                      # higher = better self-rated health
cat("anchor-test n =", nrow(b), "\n")

cat("\n== TEST 1: SES gradient in SUBJECTIVE vs OBJECTIVE health (both standardized, same people) ==\n")
mS <- lm(sr_phys_z  ~ income + educ + race + age, b)   # social gradient in SUBJECTIVE health
mO <- lm(obj_health ~ income + educ + race + age, b)   # social gradient in OBJECTIVE health
csS <- summary(mS)$coefficients; csO <- summary(mO)$coefficients
for (v in c("income","educ")) cat(sprintf("   %-7s: subjective slope %+.3f (SE %.3f)  vs  objective slope %+.3f (SE %.3f)  -> subj/obj = %.2f\n",
   v, csS[v,1], csS[v,2], csO[v,1], csO[v,2], csS[v,1]/csO[v,1]))
cat("   (subj slope >> obj slope => self-report exaggerates the social health gradient = reporting heterogeneity)\n")

## robustness: objective index as COUNT of abnormal clinical markers (BMI>=30, SBP>=140, glucose>=126)
b$n_abn <- (b$bmi>=30) + (b$sbp>=140) + (b$glucose>=126)
b$obj_health2 <- z(-b$n_abn)
mO2 <- lm(obj_health2 ~ income + educ + race + age, b); csO2 <- summary(mO2)$coefficients
cat("\n   [robustness, objective = -count(abnormal markers)]:\n")
for (v in c("income","educ")) cat(sprintf("     %-7s subj %+.3f vs obj2 %+.3f -> ratio %.2f\n", v, csS[v,1], csO2[v,1], csS[v,1]/csO2[v,1]))

cat("\n== TEST 2: does self-reported PAIN absorb the SES gradient in self-rated physical health? ==\n")
bp <- b[!is.na(b$sr_pain),]
a0 <- lm(sr_phys ~ income + educ + race + age + poly(bmi,2)+poly(sbp,2)+poly(glucose,2), bp)
a1 <- lm(sr_phys ~ income + educ + race + age + poly(bmi,2)+poly(sbp,2)+poly(glucose,2) + poly(sr_pain,2), bp)
for (v in c("income","educ")){ b0<-coef(a0)[v]; b1<-coef(a1)[v]
  cat(sprintf("   %-7s: %.3f (no pain) -> %.3f (+pain)  (%.0f%% retained)\n", v, b0, b1, 100*b1/b0)) }
cat("   (large drop => part of the 'reporting' gap is unmeasured symptom burden; small drop => not pain)\n")

cat("\nDONE\n")
