# strata_income_deprivation.R -- cohesion HR within income and area-deprivation strata on the
# rebuilt cohort (replaces an ad hoc notebook cell run before the 13 Nov rebuild, whose "<75k"
# group also contained participants with income not reported). Primary specification; within an
# income stratum the income factor is kept, so it adjusts for brackets inside the stratum.
source("R/07_revision/_common.R", local = TRUE)
prim <- "age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
fit_s <- function(s, lab){ m <- coxph(as.formula(sprintf("Surv(time_days,event) ~ z_cohesion + %s", prim)), s)
  ci <- exp(confint(m)["z_cohesion", ])
  cat(sprintf("%-40s HR %.4f (%.4f-%.4f)  n=%s ev=%s\n", lab, exp(coef(m)[["z_cohesion"]]), ci[1], ci[2], mask(m$n), mask(m$nevent))) }
d <- cox[!is.na(cox$z_cohesion), ]
hdr("Income strata (income_n: 1-9 brackets; >=75k is brackets 6-9)")
fit_s(d[!is.na(d$income_n) & d$income_n >= 6, ], "income >=$75k")
fit_s(d[!is.na(d$income_n) & d$income_n <= 5, ], "income <$75k (reported)")
fit_s(d[is.na(d$income_n), ], "income not reported")
fit_s(d[is.na(d$income_n) | d$income_n <= 5, ], "income <$75k or not reported (old definition)")
hdr("Area deprivation tertiles (ZIP3)")
dd <- d[!is.na(d$deprivation_index), ]
dd$dep_t <- cut(dd$deprivation_index, quantile(dd$deprivation_index, 0:3/3), include.lowest = TRUE, labels = c("T1 least deprived","T2","T3 most deprived"))
for (t in levels(dd$dep_t)) fit_s(dd[dd$dep_t == t, ], paste("deprivation", t))
