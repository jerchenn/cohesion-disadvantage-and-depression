# sleep_headtohead.R -- GATING CHECK for the self-report/device sleep dissociation.
# The draft compares self-reported sleep (SWB cohort) with device sleep (Fitbit cohort),
# which are DIFFERENT samples. Test both in the SAME people (swb n fitbit) to confirm the
# dissociation is within-person, not a composition artifact. If it holds -> Figure panel;
# if not (or too few) -> downgrade to a descriptive footnote, drop the mechanistic claim.
library(survival)
wb <- readRDS("wb_item_dat.rds")     # self-report sleep item: z_sleep (higher = fewer problems)
fb <- readRDS("fitbit_dat.rds")      # device sleep: efficiency, waso_min, asleep_min

fb$z_eff  <-  scale(fb$efficiency)[,1]     # higher = better
fb$z_waso <- -scale(fb$waso_min)[,1]       # higher = better (less wake)
m <- merge(wb, fb[, c("person_id","z_eff","z_waso")], by="person_id")
m <- m[!is.na(m$event) & !is.na(m$z_sleep) & !is.na(m$z_eff), ]

m$sex_c  <- factor(ifelse(m$sex %in% c("Female","Male"), m$sex, "Other"))
m$race_c <- factor(ifelse(m$race %in% c("White","Black or African American"), m$race, "Other"))
m$ethn_c <- factor(ifelse(m$ethnicity=="Hispanic or Latino","Hispanic",
                   ifelse(m$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
m$income_m <- ifelse(is.na(m$income_n), median(m$income_n,na.rm=TRUE), m$income_n)
m$educ_m   <- ifelse(is.na(m$educ_n),   median(m$educ_n,  na.rm=TRUE), m$educ_n)
covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
cat("swb n fitbit intersection: n =", nrow(m), " events =", sum(m$event), "\n\n")

hr <- function(m, v){ ci<-confint(m); sprintf("%.3f (%.3f-%.3f)  p=%.2g",
     exp(coef(m)[[v]]), exp(ci[v,1]), exp(ci[v,2]), summary(m)$coefficients[v,"Pr(>|z|)"]) }
f <- function(rhs) coxph(as.formula(paste("Surv(time_days,event) ~", rhs, "+", covs)), m)

cat("Self-reported sleep alone : ", hr(f("z_sleep"), "z_sleep"), "\n")
cat("Device sleep eff alone    : ", hr(f("z_eff"),   "z_eff"),   "\n")
cat("Device WASO alone         : ", hr(f("z_waso"),  "z_waso"),  "\n")
mm <- f("z_sleep + z_eff")
cat("\nHEAD-TO-HEAD (same people, both in model):\n")
cat("  self-report sleep : ", hr(mm, "z_sleep"), "\n")
cat("  device efficiency : ", hr(mm, "z_eff"),   "\n")
cat("\nInterpretation: if self-report z_sleep stays clearly <1 and device z_eff ~1 in the\n")
cat("joint model, the dissociation is within-person (keep as a figure panel + claim).\n")
