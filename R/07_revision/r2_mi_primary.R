# r2_mi_primary.R -- Reviewer 2, fourth point: multiple imputation as the PRIMARY handling of missing income,
# education and employment, with the imputation model and auxiliary variables reported, and a sensitivity
# analysis under departures from missing at random (delta adjustment of the imputed values).
# Same base-R predictive mean matching as mi_sensitivity.R (no mice on the Workbench image); m = 20.
# Imputation model predictors: cohesion, age, gender, race, ethnicity, log utilisation, prior EHR years,
# area deprivation (with a missing indicator), employment, the Nelson-Aalen cumulative hazard and the event
# indicator (White & Royston 2009), plus income when imputing education and education when imputing income.
# Aggregate output only. Needs: survival.
source("R/07_revision/_common.R", local = TRUE)
set.seed(1); M <- 20

d <- data.frame(time_days = cox$time_days, event = cox$event, z_cohesion = cox$z_cohesion, age = cox$age,
                sex = factor(cox$sex), race = factor(cox$race), ethnicity = factor(cox$ethnicity),
                income = as.numeric(cox$income_n), education = as.numeric(cox$educ_n),
                employment = as.character(cox$employment), log_util = cox$log_util,
                prior_ehr_y = cox$prior_ehr_days / 365.25, dep = cox$deprivation_index)
d$employment[d$employment %in% c("PMI: Skip", "PMI: Prefer Not To Answer")] <- NA
d <- d[!is.na(d$z_cohesion), ]
d$dep_miss <- as.integer(is.na(d$dep)); d$dep[is.na(d$dep)] <- median(d$dep, na.rm = TRUE)
cat(sprintf("Analytic rows: %s | missing: income %.1f%%, education %.1f%%, employment %.1f%%, deprivation %.1f%% (median-filled with indicator as a predictor)\n",
            mask(nrow(d)), 100*mean(is.na(d$income)), 100*mean(is.na(d$education)), 100*mean(is.na(d$employment)), 100*mean(d$dep_miss)))
bh <- basehaz(coxph(Surv(time_days, event) ~ 1, d), centered = FALSE)
d$na_haz <- approx(bh$time, bh$hazard, xout = d$time_days, method = "constant", rule = 2)$y
emp_levels <- names(sort(table(d$employment), decreasing = TRUE))
X0 <- model.matrix(~ z_cohesion + age + sex + race + ethnicity + log_util + prior_ehr_y + dep + dep_miss + na_haz + event, d)

pmm_impute <- function(y, X, k = 5){
  obs <- which(!is.na(y)); mis <- which(is.na(y)); if (!length(mis)) return(y)
  bb <- sample(obs, length(obs), replace = TRUE)
  b  <- qr.coef(qr(X[bb, , drop = FALSE]), y[bb]); b[is.na(b)] <- 0
  yhat <- as.vector(X %*% b)
  o <- order(yhat[obs]); yo <- yhat[obs][o]; yv <- y[obs][o]
  pos <- findInterval(yhat[mis], yo, all.inside = TRUE)
  off <- sample(-((k-1) %/% 2):(k %/% 2), length(mis), replace = TRUE)
  y[mis] <- yv[pmin(pmax(pos + off, 1L), length(yo))]; y
}
## employment: impute as the most likely level by multinomial-free PMM on a numeric code ordered by event rate
emp_code <- rank(tapply(d$event, d$employment, mean))[d$employment]   # numeric code, higher = higher crude risk
code2lev <- setNames(names(rank(tapply(d$event, d$employment, mean))), rank(tapply(d$event, d$employment, mean)))

impute_once <- function(delta_inc = 0, delta_edu = 0){
  di <- d
  ## chained: income | X, education(obs or current); education | X, income; employment | X, income, education
  edu_fill <- ifelse(is.na(d$education), median(d$education, na.rm = TRUE), d$education)
  di$income <- pmm_impute(d$income, cbind(X0, edu_fill))
  di$education <- pmm_impute(d$education, cbind(X0, di$income))
  di$income <- pmm_impute(d$income, cbind(X0, di$education))          # second pass with imputed education
  ec <- pmm_impute(as.numeric(emp_code), cbind(X0, di$income, di$education))
  di$employment <- factor(code2lev[as.character(round(ec))], levels = names(code2lev))
  ## delta adjustment (MNAR): shift IMPUTED values only; bracket scale, floored/capped
  mi_inc <- is.na(d$income); di$income[mi_inc] <- pmin(pmax(di$income[mi_inc] + delta_inc, 1), 9)
  mi_edu <- is.na(d$education); di$education[mi_edu] <- pmin(pmax(di$education[mi_edu] + delta_edu, 1), 4)
  fm <- coxph(Surv(time_days, event) ~ z_cohesion + age + sex + race + ethnicity + factor(income) + factor(education) + employment + log_util, di)
  c(coef(fm)[["z_cohesion"]], vcov(fm)["z_cohesion", "z_cohesion"], mean(di$income[mi_inc]), mean(di$education[mi_edu]))
}
pool <- function(res, lab){
  Q <- mean(res[1, ]); U <- mean(res[2, ]); B <- var(res[1, ]); Tv <- U + (1 + 1/M) * B; SE <- sqrt(Tv)
  fmi <- (1 + 1/M) * B / Tv
  cat(sprintf("%-58s HR %.3f (%.3f-%.3f)  FMI %.3f | mean imputed income bracket %.2f, education level %.2f\n",
              lab, exp(Q), exp(Q - 1.96*SE), exp(Q + 1.96*SE), fmi, mean(res[3, ]), mean(res[4, ])))
}
run_mi <- function(delta_inc = 0, delta_edu = 0){ t0 <- Sys.time(); r <- sapply(1:M, function(m) impute_once(delta_inc, delta_edu))
  cat(sprintf("  (%d imputations, %.0f s)\n", M, as.numeric(difftime(Sys.time(), t0, units = "secs")))); r }

hdr("Reference: observed-data benchmarks")
cat(sprintf("observed income: mean bracket %.2f; observed education: mean level %.2f\n", mean(d$income, na.rm = TRUE), mean(d$education, na.rm = TRUE)))
hr(coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), cox), lab = "missing-indicator specification (current Table 2 primary)")

hdr("MI as primary (MAR), m = 20, categorical income/education/employment as in the primary model")
pool(run_mi(0, 0), "MI, MAR")
hdr("MNAR sensitivity: delta adjustment of the imputed values")
pool(run_mi(-1, 0), "imputed income one bracket LOWER than MAR predicts")
pool(run_mi(-2, 0), "imputed income two brackets lower")
pool(run_mi(+1, 0), "imputed income one bracket higher")
pool(run_mi(0, -1), "imputed education one level lower")
pool(run_mi(-1, -1), "income one bracket lower and education one level lower")
cat("\nNote for the text: the delta scenarios bound the plausible departures (non-responders poorer or less educated than MAR implies). Report the range of HRs.\n")
