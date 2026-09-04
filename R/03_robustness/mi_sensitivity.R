# mi_sensitivity.R -- sensitivity: multiple imputation for missing income/education, replacing the
# missing-indicator specification. BASE-R implementation (no mice/lme4/nloptr, which fail to compile
# in some Workbench images): predictive mean matching with a Bayesian-bootstrap donor model,
# m=20 imputations, Nelson-Aalen cumulative hazard + event as predictors (White & Royston,
# Stat Med 2009), pooled by Rubin's rules. Needs only 'survival'.
#
# INPUT : cox_dat.rds (from build_cox.R)     OUTPUT: printed pooled HR
# RUN AFTER: build_cox.R.

library(survival)
cox <- readRDS("cox_dat.rds")

d <- data.frame(
  time_days=cox$time_days, event=cox$event, z_cohesion=cox$z_cohesion, age=cox$age,
  sex=factor(cox$sex), race=factor(cox$race), ethnicity=factor(cox$ethnicity),
  income=as.numeric(cox$income_n), education=as.numeric(cox$educ_n),
  employment=cox$emp_f, log_util=cox$log_util)          # emp_f keeps NA as a level (addNA)

## restrict to rows complete on the NON-imputed covariates (income/education are imputed below).
## This matches the primary complete-covariate cohort and keeps d aligned with the design matrix.
keep <- stats::complete.cases(d[, c("time_days","event","z_cohesion","age","sex","race",
                                    "ethnicity","employment","log_util")])
d <- droplevels(d[keep, ])
cat(sprintf("Analytic rows: %d | to impute: income %d (%.1f%%), education %d (%.1f%%)\n",
    nrow(d), sum(is.na(d$income)), 100*mean(is.na(d$income)),
    sum(is.na(d$education)), 100*mean(is.na(d$education))))

## Nelson-Aalen cumulative hazard as an imputation predictor (base survival; no mice)
bh <- basehaz(coxph(Surv(time_days, event) ~ 1, d), centered=FALSE)
d$na_haz <- approx(bh$time, bh$hazard, xout=d$time_days, method="constant", rule=2)$y

## fully observed predictors for the imputation models (no NAs after the restriction above)
X <- model.matrix(~ z_cohesion + age + sex + race + ethnicity + log_util + na_haz + event, d)
stopifnot(nrow(X) == nrow(d))

## fast predictive mean matching for one incomplete numeric variable.
## Bayesian bootstrap of donors injects parameter uncertainty (proper MI); nearest-donor match on
## the linear predictor keeps imputed values on the observed (ordinal) scale.
pmm_impute <- function(y, X, k=5){
  obs <- which(!is.na(y)); mis <- which(is.na(y))
  if(!length(mis)) return(y)
  bb <- sample(obs, length(obs), replace=TRUE)              # Bayesian bootstrap
  b  <- qr.coef(qr(X[bb,,drop=FALSE]), y[bb]); b[is.na(b)] <- 0
  yhat <- as.vector(X %*% b)
  o  <- order(yhat[obs]); yo <- yhat[obs][o]; yv <- y[obs][o]   # donor yhat/values, sorted
  pos <- findInterval(yhat[mis], yo, all.inside=TRUE)          # nearest donor index (O(log n))
  off <- sample(-((k-1)%/%2):(k%/%2), length(mis), replace=TRUE)
  y[mis] <- yv[pmin(pmax(pos+off, 1L), length(yo))]           # random draw among k neighbours
  y
}

set.seed(1); M <- 20; est <- var_ <- numeric(M)   # set M <- 10 to halve the runtime if needed
t0 <- Sys.time()
for(m in 1:M){
  di <- d
  di$income    <- pmm_impute(d$income,    X)
  di$education <- pmm_impute(d$education,  X)
  fm <- coxph(Surv(time_days, event) ~ z_cohesion + age + sex + race + ethnicity +
                income + education + employment + log_util, di)
  est[m]  <- coef(fm)[["z_cohesion"]]
  var_[m] <- vcov(fm)["z_cohesion","z_cohesion"]
  cat(sprintf("  imputation %2d/%d done  (HR %.3f, elapsed %.0fs)\n",
      m, M, exp(est[m]), as.numeric(difftime(Sys.time(), t0, units="secs"))))
  flush.console()
}

## Rubin's rules
Q <- mean(est); U <- mean(var_); B <- var(est); Tv <- U + (1 + 1/M)*B; SE <- sqrt(Tv)
cat(sprintf("=== MI (m=20, PMM) primary z_cohesion HR: %.3f (%.3f-%.3f) ===\n",
    exp(Q), exp(Q - 1.96*SE), exp(Q + 1.96*SE)))
cat(sprintf("Missing before imputation: income %.1f%%, education %.1f%%\n",
    100*mean(is.na(d$income)), 100*mean(is.na(d$education))))
cat("Compare: missing-indicator primary HR 0.882 (0.857-0.908); parsimonious ~0.87 (0.85-0.90).\n")
cat("Interpretation: a pooled HR close to these indicates the estimate does not depend on the\n",
    "missing-indicator handling of income and education.\n", sep="")
