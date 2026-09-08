# rake_ipw.R -- representativeness sensitivity. Rake (iterative proportional fitting)
# the analytic sample to approximate US adult (18+) margins for age, sex, race/ethnicity,
# and education, then re-estimate the primary Cox with the weights. Tests whether the
# HR is an artifact of the cohort's advantaged, non-representative composition.
library(survival)
cox <- readRDS("cox_dat3.rds")

## --- raking variables ---------------------------------------------------------
cox$age_band <- cut(cox$age, c(-Inf,29,44,64,Inf), labels=c("18-29","30-44","45-64","65+"))
cox$sex_r    <- ifelse(cox$sex %in% c("Female","Male"), cox$sex, NA)
cox$race_eth <- with(cox, ifelse(ethnicity=="Hispanic or Latino","Hispanic",
                    ifelse(race=="White","WhiteNH",
                    ifelse(race=="Black or African American","BlackNH","Other"))))
cox$educ_band<- c("1"="<HS","2"="HS","3"="SomeCol","4"="BA+")[as.character(cox$educ_n)]

rk <- cox[!is.na(cox$age_band) & !is.na(cox$sex_r) & !is.na(cox$race_eth) & !is.na(cox$educ_band) &
          !is.na(cox$z_cohesion) & !is.na(cox$log_util), ]           # complete cases for the model
cat("Raking-complete N =", nrow(rk), "of", nrow(cox), "\n")

## --- target US adult margins (approximate ACS, 18+) ---------------------------
targets <- list(
  age_band = c("18-29"=.205,"30-44"=.255,"45-64"=.325,"65+"=.215),
  sex_r    = c("Female"=.515,"Male"=.485),
  race_eth = c("Hispanic"=.17,"WhiteNH"=.60,"BlackNH"=.12,"Other"=.11),
  educ_band= c("<HS"=.11,"HS"=.28,"SomeCol"=.29,"BA+"=.32))

## --- iterative proportional fitting -------------------------------------------
w <- rep(1, nrow(rk))
for (it in 1:30) for (v in names(targets)){
  g   <- rk[[v]]; tg <- targets[[v]]
  cur <- tapply(w, g, sum); cur <- cur[names(tg)]
  f   <- (tg * sum(w)) / cur                      # per-level rescale factor
  w   <- w * f[as.character(g)]
}
w <- w / mean(w)                                   # mean 1
cap <- quantile(w, .995); w[w > cap] <- cap        # trim extreme weights
rk$w <- w
cat(sprintf("weights: median %.2f, 95th %.2f, max %.2f; effective N = %.0f\n",
    median(w), quantile(w,.95), max(w), sum(w)^2/sum(w^2)))

## --- primary Cox: unweighted vs raked (same subset, so weighting is isolated) --
rk$sex_c<-factor(rk$sex_r); rk$race_c<-factor(rk$race_eth)
rk$ethn_c<-factor(ifelse(rk$ethnicity=="Hispanic or Latino","Hispanic",
           ifelse(rk$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
rk$income_m<-ifelse(is.na(rk$income_n),median(rk$income_n,na.rm=TRUE),rk$income_n)
rk$educ_m  <-ifelse(is.na(rk$educ_n),  median(rk$educ_n,  na.rm=TRUE),rk$educ_n)
rk <- droplevels(rk)
f  <- Surv(time_days,event) ~ z_cohesion + age + sex_c + race_c + ethn_c + income_m + educ_m + log_util
hr <- function(m) sprintf("%.3f (%.3f-%.3f)", exp(coef(m)[["z_cohesion"]]),
                          exp(confint(m)["z_cohesion",1]), exp(confint(m)["z_cohesion",2]))
mu <- coxph(f, rk)                                           # unweighted (same subset)
mw <- coxph(f, rk, weights=rk$w, robust=FALSE)              # force off auto-robust (avoids buggy path)

## bootstrap the raked HR CI (naive SE discarded)
set.seed(1); B <- 300; n <- nrow(rk); bc <- numeric(B)
for (b in 1:B){ i <- sample.int(n, n, TRUE)
  fit <- tryCatch(coxph(f, rk[i,], weights=rk$w[i], robust=FALSE), error=function(e) NULL)
  bc[b] <- if (is.null(fit)) NA else coef(fit)[["z_cohesion"]] }
bc <- bc[!is.na(bc)]; ci <- exp(quantile(bc, c(.025,.975)))
saveRDS(list(hr=exp(coef(mw)[["z_cohesion"]]), lo=ci[[1]], hi=ci[[2]],
             n=nrow(rk), events=sum(rk$event)), "rake_out.rds")
cat("\nPrimary HR on raking-complete subset:\n")
cat("  unweighted:", hr(mu), " (events", mu$nevent, ")\n")
cat(sprintf("  raked     : %.3f (%.3f-%.3f)  [%d bootstrap reps]\n",
    exp(coef(mw)[["z_cohesion"]]), ci[1], ci[2], length(bc)))
