# item_cohesion_hetero.R -- which of the 5 cohesion items drive the association with incident
# depression, and which show resource-dependence (protection fading under disadvantage).
# Reports, per item: (1) marginal HR per SD; (2) independent HR in a joint 5-item model;
# (3) item x area-deprivation and item x income interactions (>1 = protection weakens under
# disadvantage). Aggregate output only.
#
# INPUT: cox_dat.rds + ds_survey + ds_zip_code_socioeconomic.  RUN AFTER build_cox.R.
# Needs: tidyverse, survival, bigrquery.

library(tidyverse); library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

## 1. pull the 5 cohesion items, recode (scale-aware), z-score each
ids <- c(help=40192463, getalong=40192411, trust=40192499, values=40192417, watchout=40192400)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(ids)[match(as.numeric(sv$question_concept_id), ids)]
wide <- tidyr::pivot_wider(sv[,c("person_id","item","answer")], names_from=item, values_from=answer)
map5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","Skip","Don't know","Prefer not to answer")
rc <- function(x,m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in c("help","getalong","trust","values")) wide[[nm]] <- rc(wide[[nm]], map5)
wide$watchout <- rc(wide$watchout, map4)
for (nm in names(ids)) wide[[paste0("z_",nm)]] <- as.numeric(scale(wide[[nm]]))
d <- merge(cox, wide[, c("person_id", paste0("z_", names(ids)))], by="person_id")

## parsimonious covariate block (matches the figures)
collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
d <- collapse(d)
COV <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
zitems <- paste0("z_", names(ids))
hrci <- function(m,v) sprintf("%.2f (%.2f-%.2f)", exp(coef(m)[[v]]),
                              exp(confint(m)[v,1]), exp(confint(m)[v,2]))

## 2. marginal HR per SD (each item alone)
cat("=== Marginal HR per SD, each cohesion item (adjusted) ===\n")
for (v in zitems){ m <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "+", COV)), d)
  cat(sprintf("  %-9s  %s\n", sub("z_","",v), hrci(m,v))) }

## 3. joint 5-item model -> independent contributions (note: items are correlated)
cat("\n=== Joint 5-item model (independent HRs; items intercorrelated) ===\n")
mj <- coxph(as.formula(paste("Surv(time_days,event) ~", paste(zitems,collapse=" + "), "+", COV)), d)
for (v in zitems) cat(sprintf("  %-9s  %s\n", sub("z_","",v), hrci(mj,v)))

## 4. item x disadvantage interactions (>1 = protection weakens as disadvantage rises)
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
gd <- merge(d, geo, by="person_id"); gd <- gd[!is.na(gd$deprivation_index),]; gd$z_dep <- scale(gd$deprivation_index)[,1]
d$z_inc <- scale(d$income_m)[,1]
cat("\n=== Item x disadvantage interaction (per SD; >1 = protection fades under disadvantage) ===\n")
cat(sprintf("  %-9s  %-24s  %-24s\n","item","x area deprivation","x lower income"))
for (v in zitems){
  md <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "* z_dep +", COV, "+ z_dep")), gd)
  id <- grep(":", names(coef(md)), value=TRUE)[1]
  mi <- coxph(as.formula(paste("Surv(time_days,event) ~", v,
              "* z_inc + age + sex_c + race_c + ethn_c + educ_m + log_util")), d)
  ii <- grep(":", names(coef(mi)), value=TRUE)[1]
  cat(sprintf("  %-9s  HR %.3f (p=%.2g)%s  HR %.3f /SD lower income\n",
      sub("z_","",v), exp(coef(md)[[id]]), summary(md)$coefficients[id,"Pr(>|z|)"],
      strrep(" ", 5), 1/exp(coef(mi)[[ii]])))
}
cat("\nInterpretation: larger marginal/joint HR below 1 = stronger driver; interaction HR farther\n",
    "above 1 = that item's protection fades most under disadvantage.\n", sep="")
