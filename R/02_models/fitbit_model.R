# fitbit_model.R -- two directions on the device layer:
#  (1) cohesion -> objective behaviour (activity/sleep)  [behavioural embedding]
#  (2) objective behaviour -> incident clinical depression [non-self-report replication]
library(survival); library(sandwich); library(lmtest)
fb <- readRDS("fitbit_dat.rds")

metrics <- c(steps="steps", sed_min="sed_min", mvpa_min="mvpa_min",
             asleep_min="asleep_min", efficiency="efficiency", waso_min="waso_min")
# direction: +1 = healthier (steps/mvpa/efficiency/asleep up good; sedentary/waso up bad)
healthy_sign <- c(steps=1, sed_min=-1, mvpa_min=1, asleep_min=1, efficiency=1, waso_min=-1)
for (v in metrics) fb[[paste0("z_",v)]] <- healthy_sign[v] * scale(fb[[v]])[,1]

## collapsed covariates (kill rare-level separation in the Cox nuisance terms)
fb$sex_c  <- factor(ifelse(fb$sex %in% c("Female","Male"), fb$sex, "Other"))
fb$race_c <- factor(ifelse(fb$race %in% c("White","Black or African American"), fb$race, "Other"))
fb$ethn_c <- factor(ifelse(fb$ethnicity=="Hispanic or Latino","Hispanic",
                    ifelse(fb$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
fb$income_m <- ifelse(is.na(fb$income_n), median(fb$income_n, na.rm=TRUE), fb$income_n)
fb$educ_m   <- ifelse(is.na(fb$educ_n),   median(fb$educ_n,   na.rm=TRUE), fb$educ_n)
covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"

## (1) cohesion -> behaviour (native units per SD cohesion, HC3) ----
cat("=== (1) Cohesion -> device behaviour (native units per +1 SD cohesion) ===\n")
for (v in metrics){
  m <- lm(as.formula(paste(v, "~ z_cohesion +", covs)), fb)
  ct <- coeftest(m, vcov=vcovHC(m,"HC3"))["z_cohesion",]
  cat(sprintf("  %-11s %+8.3f  (SE %.3f, p %.2g, N=%d)\n", v, ct[1], ct[2], ct[4], length(m$residuals)))
}

## (2) behaviour -> incident depression (HR per SD, oriented healthier => HR<1) ----
cat("\n=== (2) Device behaviour -> incident depression (HR per +1 SD healthier) ===\n")
hr1 <- function(v, dat){
  m <- coxph(as.formula(paste("Surv(time_days, event) ~", v, "+", covs)), dat)
  ci <- confint(m)
  sprintf("%.3f (%.3f-%.3f)  p=%.2g  [ev %d]", exp(coef(m)[[v]]), exp(ci[v,1]), exp(ci[v,2]),
          summary(m)$coefficients[v,"Pr(>|z|)"], m$nevent)
}
for (v in metrics){
  zv <- paste0("z_", v)
  cat(sprintf("  %-11s full: %s\n", v, hr1(zv, fb)))
  cat(sprintf("  %-11s 180d: %s\n", "", hr1(zv, fb[which(fb$time_days>180),])))
}

## (3) do activity + sleep predict independently? joint model ----
cat("\n=== (3) Joint: steps + sleep efficiency, mutually adjusted ===\n")
mj <- coxph(as.formula(paste("Surv(time_days, event) ~ z_steps + z_efficiency +", covs)), fb)
for (v in c("z_steps","z_efficiency")) cat(sprintf("  %-12s HR %.3f (%.3f-%.3f)  p=%.2g\n", v,
    exp(coef(mj)[[v]]), exp(confint(mj)[v,1]), exp(confint(mj)[v,2]),
    summary(mj)$coefficients[v,"Pr(>|z|)"]))
cat("  events:", mj$nevent, " N:", mj$n, "\n")
