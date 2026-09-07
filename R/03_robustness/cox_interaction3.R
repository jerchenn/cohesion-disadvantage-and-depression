# cox_interaction3.R -- referee concern 6 (+ PH diagnostics). On the corrected cohort (cox_dat3.rds):
#  - cohesion x area-deprivation and cohesion x income interactions, each with the PRIMARY adjustment set
#    (not the parsimonious spec the reviewer flagged), ZIP3 cluster-robust
#  - both interactions in ONE model to test whether each survives adjustment for the other (are the two
#    disadvantage axes independent modifiers?)
#  - full Schoenfeld PH table to identify which covariate drives the global test
# Aggregate only. Needs: survival.

library(survival)
d <- readRDS("cox_dat3.rds")
d <- d[!is.na(d$deprivation_index), ]
z <- function(v){ v[is.na(v)] <- median(v,na.rm=TRUE); as.numeric(scale(v)) }
d$dep_z <- z(d$deprivation_index); d$income_z <- z(d$income_n); d$clus <- factor(d$deprivation_index)
prim <- "age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
prim_noinc <- "age + sex + race + ethnicity + educ_f + emp_f + log_util"
ix <- function(m, term, lab){ b<-summary(m)$coefficients[term,]
  cat(sprintf("%-40s HR %.3f (%.3f-%.3f)  P=%s\n", lab, exp(b[["coef"]]),
      exp(b[["coef"]]-1.96*b[["robust se"]]), exp(b[["coef"]]+1.96*b[["robust se"]]),
      signif(b[["Pr(>|z|)"]],2))) }

cat("=== Interaction, primary adjustment set, ZIP3 cluster-robust ===\n")
ma <- coxph(as.formula(sprintf("Surv(time_days,event) ~ z_cohesion*dep_z + %s + cluster(clus)", prim)), d)
ix(ma, "z_cohesion:dep_z", "cohesion x area-deprivation (alone)")
mi <- coxph(as.formula(sprintf("Surv(time_days,event) ~ z_cohesion*income_z + %s + cluster(clus)", prim_noinc)), d)
ix(mi, "z_cohesion:income_z", "cohesion x income (alone)")

cat("\n=== Both axes in one model (do they survive mutual adjustment?) ===\n")
mb <- coxph(as.formula(sprintf("Surv(time_days,event) ~ z_cohesion*dep_z + z_cohesion*income_z + %s + cluster(clus)", prim_noinc)), d)
ix(mb, "z_cohesion:dep_z", "cohesion x area-deprivation (mutual)")
ix(mb, "z_cohesion:income_z", "cohesion x income (mutual)")
cat("  (dep_z higher = more deprived; income_z higher = higher income. Positive dep interaction or\n")
cat("   negative income interaction => protection weakens as disadvantage rises.)\n")

cat("\n=== Full Schoenfeld PH table (identify the global-test driver) ===\n")
print(round(cox.zph(coxph(as.formula(sprintf("Surv(time_days,event) ~ z_cohesion + %s", prim)), d))$table, 4))
