# build_al.R -- exploratory expansion of the biology layer: cohesion -> ALLOSTATIC LOAD.
# AL = multi-system count of biomarkers in the high-risk range (Seeman/McEwen style),
# across cardiovascular, metabolic, and inflammatory systems. CROSS-SECTIONAL, lab-selected,
# more selected than HbA1c (needs many panels) -> exploratory supplement only.
library(tidyverse); library(bigrquery); library(sandwich); library(lmtest)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

## biomarker set: name pattern, plausibility bounds, direction (dir=-1 => low value is high-risk)
mk <- list(
  sbp =list(pat="systolic blood pressure", lo=70, hi=250, dir= 1),
  dbp =list(pat="diastolic blood pressure",lo=40, hi=150, dir= 1),
  hr  =list(pat="heart rate",              lo=30, hi=200, dir= 1),
  a1c =list(pat="hemoglobin a1c",          lo=3,  hi=20,  dir= 1),
  chol=list(pat="cholesterol [mass/volume] in serum or plasma", lo=50, hi=500, dir=1),
  hdl =list(pat="cholesterol in hdl",      lo=10, hi=150, dir=-1),
  trig=list(pat="triglyceride",            lo=20, hi=2000,dir= 1),
  crp =list(pat="c reactive protein",      lo=0,  hi=100, dir= 1),
  bmi =list(pat="body mass index",         lo=12, hi=80,  dir= 1))

vals <- data.frame(person_id=cox$person_id)
cat("=== biomarker concepts chosen (top standard concept by coverage) ===\n")
for (nm in names(mk)){
  m <- mk[[nm]]
  top <- run_sql(sprintf(
    "SELECT m.measurement_concept_id cid, ANY_VALUE(c.concept_name) nm, COUNT(DISTINCT m.person_id) n
     FROM `%s.measurement` m JOIN `%s.concept` c ON c.concept_id=m.measurement_concept_id
     WHERE LOWER(c.concept_name) LIKE '%%%s%%' AND c.standard_concept='S' AND c.domain_id='Measurement'
     GROUP BY cid ORDER BY n DESC LIMIT 1", cdr, cdr, m$pat))
  if (nrow(top)==0){ cat(sprintf("  %-5s: NONE FOUND for '%s'\n", nm, m$pat)); next }
  cid <- top$cid[1]
  cat(sprintf("  %-5s concept %s  %-45s n=%s\n", nm, cid, substr(top$nm[1],1,45), top$n[1]))
  v <- run_sql(sprintf(
    "SELECT person_id, AVG(value_as_number) val FROM `%s.measurement`
     WHERE measurement_concept_id=%s AND value_as_number BETWEEN %f AND %f GROUP BY person_id",
    cdr, cid, m$lo, m$hi))
  v$val <- as.numeric(v$val); names(v)[2] <- nm
  vals <- merge(vals, v, by="person_id", all.x=TRUE)
}

## high-risk flags: top quartile (bottom quartile for dir=-1, e.g. HDL)
flag <- data.frame(person_id=vals$person_id)
for (nm in names(mk)){
  if (!nm %in% names(vals)) next
  x <- vals[[nm]]; d <- mk[[nm]]$dir
  cut <- if (d==1) quantile(x, .75, na.rm=TRUE) else quantile(x, .25, na.rm=TRUE)
  flag[[nm]] <- if (d==1) as.integer(x >= cut) else as.integer(x <= cut)
}
bm <- names(mk)[names(mk) %in% names(flag)]
navail <- rowSums(!is.na(flag[bm]))
al <- rowSums(flag[bm], na.rm=TRUE)
al[navail < 6] <- NA                          # require >=6 systems measured
d <- merge(cox, data.frame(person_id=vals$person_id, al=al, n_bm=navail,
                           bmi_val=vals$bmi), by="person_id", all.x=TRUE)
cat(sprintf("\nAllostatic load computable for %d / %d (%.1f%%); median markers/person = %.0f\n",
    sum(!is.na(d$al)), nrow(d), 100*mean(!is.na(d$al)), median(navail[navail>0])))
cat("AL distribution (0-9 systems):\n"); print(table(d$al))

## selection check + models
d$has_al <- !is.na(d$al)
cat("\n=== Selection: AL-havers vs not ===\n")
for (v in c("z_cohesion","age","log_util")) cat(sprintf("  %-11s havers %.3f | non %.3f\n",
    v, mean(d[[v]][d$has_al],na.rm=TRUE), mean(d[[v]][!d$has_al],na.rm=TRUE)))

covs <- "age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"
fitp <- function(m,lab){ct<-coeftest(m,vcov=vcovHC(m,"HC3"))["z_cohesion",]
  cat(sprintf("  %-24s beta %+.4f  SE %.4f  p %.2g  (N=%d)\n",lab,ct[1],ct[2],ct[4],length(m$residuals)))}
cat("\n=== Cohesion -> allostatic load (systems dysregulated per SD cohesion) ===\n")
fitp(lm(as.formula(paste("al ~ z_cohesion +",covs)), d), "base (demog+SES+util)")
fitp(lm(as.formula(paste("al ~ z_cohesion + bmi_val +",covs)), d), "+ BMI")
saveRDS(d[c("person_id","al","n_bm")], "al_dat.rds")
