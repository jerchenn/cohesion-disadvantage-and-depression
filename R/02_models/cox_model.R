# cox_model.R -- Study A: cohesion -> incident depression (Cox PH)
library(survival)
cox_dat <- readRDS("cox_dat.rds")

## main model: HR per SD cohesion, adjusted (incl. healthcare utilization for ascertainment bias)
m <- coxph(Surv(time_days, event) ~ z_cohesion + age + sex + race + ethnicity +
             income_f + educ_f + emp_f + log_util, data = cox_dat)
cat("=== Main model (N =", m$n, ", events =", m$nevent, ") ===\n")
print(summary(m)$coefficients["z_cohesion", ])
cat("HR per SD cohesion:", round(exp(coef(m)[["z_cohesion"]]), 4),
    " 95% CI:", paste(round(exp(confint(m)["z_cohesion", ]), 4), collapse=" - "), "\n")

## sensitivity 1: 6-month washout (drop events in first 180 days -> guards against prevalent-detection)
d2 <- cox_dat[cox_dat$time_days > 180, ]
m2 <- coxph(Surv(time_days, event) ~ z_cohesion + age + sex + race + ethnicity +
              income_f + educ_f + emp_f + log_util, data = d2)
cat("\n=== 6-month washout (events =", m2$nevent, ") ===\n")
cat("HR per SD cohesion:", round(exp(coef(m2)[["z_cohesion"]]), 4),
    " 95% CI:", paste(round(exp(confint(m2)["z_cohesion", ]), 4), collapse=" - "), "\n")

## proportional-hazards check for the exposure
cat("\n=== PH test (Schoenfeld) for z_cohesion ===\n")
print(cox.zph(m)["z_cohesion"]$table)
