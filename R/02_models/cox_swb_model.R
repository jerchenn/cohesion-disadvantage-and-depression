# cox_swb_model.R -- SWB convergence layer (robust covariates + 180-day washout)
library(survival)
swb_dat <- readRDS("swb_dat.rds")

# collapsed factors (>=2 levels guaranteed) + numeric median-filled SES
swb_dat$sex_c  <- factor(ifelse(swb_dat$sex %in% c("Female","Male"), swb_dat$sex, "Other"))
swb_dat$race_c <- factor(ifelse(swb_dat$race %in% c("White","Black or African American"),
                                swb_dat$race, "Other"))
swb_dat$ethn_c <- factor(ifelse(swb_dat$ethnicity == "Hispanic or Latino", "Hispanic",
                         ifelse(swb_dat$ethnicity == "Not Hispanic or Latino", "NotHispanic", "Other")))
swb_dat$income_m <- ifelse(is.na(swb_dat$income_n), median(swb_dat$income_n, na.rm=TRUE), swb_dat$income_n)
swb_dat$educ_m   <- ifelse(is.na(swb_dat$educ_n),   median(swb_dat$educ_n,   na.rm=TRUE), swb_dat$educ_n)

cat("levels -> sex_c:", nlevels(swb_dat$sex_c), " race_c:", nlevels(swb_dat$race_c),
    " ethn_c:", nlevels(swb_dat$ethn_c), "\n\n")

covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
hr <- function(m, v) sprintf("%.4f (%.4f - %.4f)",
                             exp(coef(m)[[v]]), exp(confint(m)[v,1]), exp(confint(m)[v,2]))

## Test 2: SWB -> incident depression, full + 180-day washout
m_wb  <- coxph(as.formula(paste("Surv(time_days, event) ~ z_wb +", covs)), swb_dat)
m_wb2 <- coxph(as.formula(paste("Surv(time_days, event) ~ z_wb +", covs)),
               swb_dat[swb_dat$time_days > 180, ])
cat("=== SWB -> incident depression ===\n")
cat("Full    (events", m_wb$nevent,  "): HR per SD wb =", hr(m_wb,  "z_wb"), "\n")
cat("Washout (events", m_wb2$nevent, "): HR per SD wb =", hr(m_wb2, "z_wb"), "\n")

## Test 1: cohesion -> SWB -> depression (difference method)
m_tot <- coxph(as.formula(paste("Surv(time_days, event) ~ z_cohesion +", covs)), swb_dat)
m_dir <- coxph(as.formula(paste("Surv(time_days, event) ~ z_cohesion + z_wb +", covs)), swb_dat)
lt <- coef(m_tot)[["z_cohesion"]]; ld <- coef(m_dir)[["z_cohesion"]]
cat("\n=== Cohesion -> incident depression: total vs direct (holding SWB) ===\n")
cat("Total  HR (cohesion):       ", hr(m_tot, "z_cohesion"), "\n")
cat("Direct HR (cohesion | SWB): ", hr(m_dir, "z_cohesion"), "\n")
cat("Approx proportion mediated by SWB:", round((lt - ld)/lt, 3),
    " (construct-overlap + HR non-collapsibility caveats)\n")
