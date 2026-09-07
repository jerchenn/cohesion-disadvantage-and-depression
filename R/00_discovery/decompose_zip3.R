# decompose_zip3.R -- EXPLORATORY. Break the ZIP3 DEPRIVATION_INDEX into its component structural
# facets and ask, for each: (A) how it loads on the composite index, (B) its association with the
# two driving cohesion items (trust, help), (C) its own association with incident depression, and
# (D) whether it drives the cohesion resource-dependence (cohesion x facet interaction; HR > 1 =
# cohesion's protection fades as that facet worsens). Each facet is oriented so higher = more
# deprived (sign taken from its correlation with the composite index). Aggregate output only.
#
# INPUT: cox_dat.rds + ds_survey + ds_zip_code_socioeconomic.  RUN AFTER build_cox.R.
# Needs: tidyverse, survival, bigrquery.

library(tidyverse); library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}

## ---------- pull ZIP3 components + trust/help items ----------
comp <- c("FRACTION_POVERTY","FRACTION_VACANT_HOUSING","FRACTION_NO_HEALTH_INS",
          "FRACTION_ASSISTED_INCOME","FRACTION_HIGH_SCHOOL_EDU","MEDIAN_INCOME","DEPRIVATION_INDEX")
z <- run_sql(sprintf("SELECT person_id, %s FROM `%s.ds_zip_code_socioeconomic`",
                     paste(comp, collapse=","), cdr))
z <- z[!duplicated(z$person_id), ]
for (c in comp) z[[c]] <- as.numeric(z[[c]])

th <- c(trust=40192499, help=40192463)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(th, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(th)[match(as.numeric(sv$question_concept_id), th)]
wth <- tidyr::pivot_wider(sv[,c("person_id","item","answer")], names_from=item, values_from=answer)
map5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Dont Know","PMI: Prefer Not To Answer",
             "Skip","Don't know","Prefer not to answer")
rc <- function(x){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(map5[x]) }
wth$z_trust <- as.numeric(scale(rc(wth$trust))); wth$z_help <- as.numeric(scale(rc(wth$help)))

d <- collapse(merge(cox, z, by="person_id"))
d <- merge(d, wth[, c("person_id","z_trust","z_help")], by="person_id", all.x=TRUE)

## orient each facet so + = aligned with the composite (higher = more deprived); z-score
di <- "DEPRIVATION_INDEX"; zc <- setdiff(comp, di)
for (c in comp){ s <- sign(cor(d[[c]], d[[di]], use="complete.obs")); if (s==0) s <- 1
                 d[[paste0("z_",c)]] <- s * as.numeric(scale(d[[c]])) }
lab <- c(FRACTION_POVERTY="poverty", FRACTION_VACANT_HOUSING="vacant housing",
         FRACTION_NO_HEALTH_INS="no health insurance", FRACTION_ASSISTED_INCOME="assisted income",
         FRACTION_HIGH_SCHOOL_EDU="high-school education", MEDIAN_INCOME="median income (rev.)",
         DEPRIVATION_INDEX="composite index")

## ---------- (A) what the composite index is made of ----------
cat("=== (A) Facet correlation with the composite DEPRIVATION_INDEX (raw sign shows orientation) ===\n")
for (c in zc) cat(sprintf("  %-22s r=%+.2f\n", lab[c], cor(d[[c]], d[[di]], use="complete.obs")))

## ---------- (B) facet -> trust / help (std beta, adj age + person income) ----------
cat("\n=== (B) Facet -> trust / help (std beta per SD more-deprived, adj. age + person income) ===\n")
cat(sprintf("  %-22s  %8s  %8s\n","facet","->trust","->help"))
for (c in comp){ v <- paste0("z_",c)
  bt <- coef(lm(as.formula(paste("z_trust ~", v, "+ age + income_m")), d))[[v]]
  bh <- coef(lm(as.formula(paste("z_help  ~", v, "+ age + income_m")), d))[[v]]
  cat(sprintf("  %-22s  %+8.3f  %+8.3f\n", lab[c], bt, bh)) }

## ---------- (C) facet -> incident depression (HR per SD more-deprived, adjusted) ----------
COV <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
cat("\n=== (C) Facet -> incident depression (HR per SD more-deprived, adjusted) ===\n")
for (c in comp){ v <- paste0("z_",c)
  m <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "+", COV)), d[!is.na(d[[v]]),])
  cat(sprintf("  %-22s HR %.3f (%.3f-%.3f, p=%.2g)\n", lab[c], exp(coef(m)[[v]]),
      exp(confint(m)[v,1]), exp(confint(m)[v,2]), summary(m)$coefficients[v,"Pr(>|z|)"])) }

## ---------- (D) cohesion x facet interaction: which facet drives the resource-dependence ----------
cat("\n=== (D) Cohesion x facet interaction (HR>1 = cohesion protection fades as facet worsens) ===\n")
for (c in comp){ v <- paste0("z_",c); sub <- d[!is.na(d[[v]]) & !is.na(d$z_cohesion),]
  m <- coxph(as.formula(paste("Surv(time_days,event) ~ z_cohesion *", v, "+", COV, "+", v)), sub)
  iv <- grep(":", names(coef(m)), value=TRUE)[1]
  cat(sprintf("  %-22s int HR %.3f (%.3f-%.3f, p=%.2g)\n", lab[c], exp(coef(m)[[iv]]),
      exp(confint(m)[iv,1]), exp(confint(m)[iv,2]), summary(m)$coefficients[iv,"Pr(>|z|)"])) }

cat("\nNote: facets are collinear (they compose the index); read (D) marginally, not jointly.\n")
