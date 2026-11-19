# r2_reri_subsamples_vif.R -- Reviewer 2, fifth point and the two final points.
#  (a) effect modification on the ADDITIVE scale: RERI, attributable proportion and synergy index for
#      low cohesion (bottom tertile) x high area deprivation (top tertile) and low cohesion x low income (<$35k),
#      delta-method CIs; multiplicative interaction P values with a Bonferroni correction for the two tests
#  (b) the cohesion -> incident depression association estimated WITHIN the SWB and Fitbit subsamples with the
#      primary specification, next to the main cohort, and how the subsamples differ
#  (c) collinearity diagnostics for the joint five-item model (VIF, condition number)
# Primary specification unchanged. Aggregate output only. Needs: survival.
source("R/07_revision/_common.R", local = TRUE)

## ---- (a) additive interaction -------------------------------------------------------------------
reri <- function(m, a, b, ab, lab){
  cf <- coef(m)[c(a, b, ab)]; V <- vcov(m)[c(a, b, ab), c(a, b, ab)]
  h11 <- exp(sum(cf)); h10 <- exp(cf[[1]]); h01 <- exp(cf[[2]])
  R <- h11 - h10 - h01 + 1
  g <- c(h11 - h10, h11 - h01, h11)                     # dR/d(b_a, b_b, b_ab)
  se <- sqrt(as.numeric(t(g) %*% V %*% g))
  AP <- R / h11; S <- (h11 - 1) / ((h10 - 1) + (h01 - 1))
  p_mult <- 2*pnorm(-abs(cf[[3]] / sqrt(V[3,3])))
  cat(sprintf("\n%s\n  HR exposed to A only %.2f, B only %.2f, both %.2f (reference: neither)\n", lab, h10, h01, h11))
  cat(sprintf("  RERI %.4f (95%% CI %.4f to %.4f); AP %.4f; synergy index %.4f\n", R, R - 1.96*se, R + 1.96*se, AP, S))
  cat(sprintf("  multiplicative interaction HR %.3f, P=%.3f (Bonferroni for 2 tests: P=%.3f)\n", exp(cf[[3]]), p_mult, min(1, 2*p_mult)))
  cat(sprintf("  n=%s events=%s\n", mask(m$n), mask(m$nevent)))
}
d <- cox[!is.na(cox$z_cohesion), ]
d$low_coh <- as.integer(d$z_cohesion <= quantile(d$z_cohesion, 1/3))
dd <- d[!is.na(d$deprivation_index), ]
dd$high_dep <- as.integer(dd$deprivation_index >= quantile(dd$deprivation_index, 2/3))
hdr("(a) Additive interaction")
cat(sprintf("A = low cohesion (bottom tertile): %s participants; B1 = high area deprivation (top tertile): %s; B2 = income <$35k: %s\n",
            mask(sum(d$low_coh)), mask(sum(dd$high_dep)), mask(sum(d$income_n <= 3, na.rm = TRUE))))
m1 <- fit(paste("low_coh * high_dep +", RHS_NOCOH), dd)
reri(m1, "low_coh", "high_dep", "low_coh:high_dep", "Low cohesion x high area deprivation (primary covariates)")
di <- d[!is.na(d$income_n), ]; di$low_inc <- as.integer(di$income_n <= 3)
m2 <- fit(paste("low_coh * low_inc +", sub("income_f \\+ ", "", RHS_NOCOH)), di)
reri(m2, "low_coh", "low_inc", "low_coh:low_inc", "Low cohesion x low income (<$35k) (primary covariates without the income factor)")
## continuous versions as in the paper (multiplicative, for the multiplicity note)
dd$dep_z <- as.numeric(scale(dd$deprivation_index)); di$income_z <- as.numeric(scale(di$income_n))
mc1 <- fit(paste("z_cohesion * dep_z +", RHS_NOCOH), dd); mc2 <- fit(paste("z_cohesion * income_z +", sub("income_f \\+ ", "", RHS_NOCOH)), di)
for (x in list(list(mc1, "z_cohesion:dep_z", "cohesion x deprivation, continuous"), list(mc2, "z_cohesion:income_z", "cohesion x income, continuous"))){
  b <- coef(x[[1]])[[x[[2]]]]; se <- sqrt(vcov(x[[1]])[x[[2]], x[[2]]]); p <- 2*pnorm(-abs(b/se))
  cat(sprintf("%-45s interaction HR %.3f (%.3f-%.3f) P=%.4f, Bonferroni P=%.4f\n", x[[3]], exp(b), exp(b-1.96*se), exp(b+1.96*se), p, min(1, 2*p)))
}

## ---- (b) cohesion -> depression within the subsamples ---------------------------------------------
hdr("(b) Cohesion -> incident depression within each subsample, primary specification")
hr(fit(RHS, cox), lab = "main cohort")
wb <- readRDS("wb_item_dat.rds"); wb_ids <- as.character(wb$person_id[!is.na(wb$event)])
fb <- readRDS("fitbit_dat.rds");  fb_ids <- as.character(fb$person_id[!is.na(fb$steps) & !is.na(fb$event)])
sub_wb <- cox[cox$person_id %in% wb_ids, ]; sub_fb <- cox[cox$person_id %in% fb_ids, ]
hr(fit(RHS, sub_wb), lab = "SWB subcohort (answered the affect survey)")
hr(fit(RHS, sub_fb), lab = "Fitbit subcohort (activity data)")
hr(fit(RHS, cox[!(cox$person_id %in% wb_ids), ]), lab = "main cohort WITHOUT the SWB subcohort")
hr(fit(RHS, cox[!(cox$person_id %in% fb_ids), ]), lab = "main cohort WITHOUT the Fitbit subcohort")
cat("\nhow the subsamples differ from the main cohort:\n")
desc <- function(x) c(n = mask(nrow(x)), event_pct = sprintf("%.2f", 100*mean(x$event)), age = sprintf("%.1f", mean(x$age)),
                      female_pct = sprintf("%.1f", 100*mean(x$sex == "Female")), white_pct = sprintf("%.1f", 100*mean(x$race == "White")),
                      income_ge75k_pct = sprintf("%.1f", 100*mean(x$income_n >= 6, na.rm = TRUE)), BAplus_pct = sprintf("%.1f", 100*mean(x$educ_n == 4, na.rm = TRUE)),
                      cohesion_mean = sprintf("%.2f", mean(x$z_cohesion, na.rm = TRUE)), fu_median_y = sprintf("%.2f", median(x$time_days)/365.25))
print(rbind(main = desc(cox), SWB = desc(sub_wb), Fitbit = desc(sub_fb)))

## ---- (c) collinearity in the joint item model --------------------------------------------------------
hdr("(c) Collinearity diagnostics for the joint five-item model")
items <- c("help","getalong","trust","values","watchout")
Z <- scale(cox[complete.cases(cox[, items]), items])
vif <- sapply(items, function(v) 1 / (1 - summary(lm(Z[, v] ~ Z[, setdiff(items, v)]))$r.squared))
print(round(vif, 2))
ev <- eigen(cor(Z))$values
cat(sprintf("condition number (sqrt of max/min eigenvalue of the item correlation matrix): %.2f; eigenvalues %s\n", sqrt(max(ev)/min(ev)), paste(round(ev, 2), collapse = ", ")))
cat("VIF below 5 and a condition number below 10 are the usual thresholds; report them in the eFigure legend and frame the joint model as hypothesis generating.\n")
