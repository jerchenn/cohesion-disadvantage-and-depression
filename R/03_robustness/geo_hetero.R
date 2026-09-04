# geo_hetero.R -- does cohesion's protection vary by AREA deprivation? (place-based
# targeting question, mirroring the individual-income anti-buffering finding).
# Uses ds_zip_code_socioeconomic (person-level, ZIP3 resolution) DEPRIVATION_INDEX.
library(tidyverse); library(bigrquery); library(survival)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")
cox <- readRDS("cox_dat.rds")

geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index)
geo <- geo[!duplicated(geo$person_id), ]
d <- merge(cox, geo, by="person_id")
d <- d[!is.na(d$z_cohesion) & !is.na(d$deprivation_index), ]
d$z_dep <- scale(d$deprivation_index)[,1]
d$sex_c  <- factor(ifelse(d$sex %in% c("Female","Male"), d$sex, "Other"))
d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                   ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
covs <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util + z_dep"
cat("Geo-linked cohort n =", nrow(d), " events =", sum(d$event), "\n")

## interaction cohesion x area deprivation (>1 => cohesion LESS protective in deprived areas)
m <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion*z_dep +", covs)), d)
rn <- "z_cohesion:z_dep"
cat(sprintf("cohesion x area-deprivation interaction HR %.3f (p=%.3g)\n",
    exp(coef(m)[[rn]]), summary(m)$coefficients[rn,"Pr(>|z|)"]))

## cohesion HR by area-deprivation tertile
d$dep_t <- cut(d$deprivation_index, quantile(d$deprivation_index,0:3/3,na.rm=TRUE),
               include.lowest=TRUE, labels=c("low_dep","mid_dep","high_dep"))
for (t in levels(d$dep_t)){
  mm <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion +", covs)), d[d$dep_t==t,])
  cat(sprintf("  %-9s cohesion HR %.3f (%.3f-%.3f)  [ev %d]\n", t, exp(coef(mm)[["z_cohesion"]]),
      exp(confint(mm)["z_cohesion",1]), exp(confint(mm)["z_cohesion",2]), mm$nevent))
}
