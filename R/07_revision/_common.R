# _common.R -- shared pieces for the revision scripts (sourced with local=TRUE from each script).
# Aggregate-only printing: mask() hides any count below 20. The primary specification is unchanged
# from cox_model3.R; every script here is an ADDITIONAL sensitivity analysis.
suppressPackageStartupMessages(library(survival))
mask <- function(x){ x <- as.numeric(x); ifelse(!is.na(x) & x > 0 & x < 20, "<20", format(x, big.mark=",", trim=TRUE)) }
cox <- readRDS("cox_dat3.rds"); cox$person_id <- as.character(cox$person_id)
RHS <- "z_cohesion + age + sex + race + ethnicity + income_f + educ_f + emp_f + log_util"   # primary, as in cox_model3.R
RHS_NOCOH <- sub("z_cohesion \\+ ", "", RHS)
fit <- function(rhs, d, ...) coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", rhs)), d, ...)
hr <- function(m, term = "z_cohesion", lab = ""){
  b <- coef(m)[[term]]; se <- sqrt(vcov(m)[term, term])
  cat(sprintf("%-62s HR %.3f (%.3f-%.3f)  P=%.2g  n=%s ev=%s\n", lab, exp(b), exp(b - 1.96*se), exp(b + 1.96*se),
              2*pnorm(-abs(b/se)), mask(m$n), mask(m$nevent)))
  invisible(c(b = b, se = se))
}
hdr <- function(s) cat("\n==========", s, "==========\n")
## BigQuery, only for the scripts that need it
## assigns cdr, bq_proj and run_sql in the calling script's environment (not via <<-, which can hit
## a locked binding of the same name exported by an attached package)
bq <- function(env = parent.frame()){
  suppressPackageStartupMessages(library(bigrquery))
  p <- Sys.getenv("GOOGLE_PROJECT")
  assign("cdr", Sys.getenv("WORKSPACE_CDR"), envir = env)
  assign("bq_proj", p, envir = env)
  assign("run_sql", function(q) bq_table_download(bq_project_query(p, q), bigint = "character"), envir = env)
}
## all standard SNOMED Condition descendants of the named parent concepts (same construction as the outcome set)
desc_sql <- function(parent_names) sprintf(
  "SELECT DISTINCT ca.descendant_concept_id AS concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` ca ON ca.ancestor_concept_id = p.concept_id JOIN `%s.concept` c ON c.concept_id = ca.descendant_concept_id WHERE p.vocabulary_id = 'SNOMED' AND p.standard_concept = 'S' AND p.domain_id = 'Condition' AND LOWER(p.concept_name) IN (%s) AND c.standard_concept = 'S'",
  cdr, cdr, cdr, paste(sprintf("'%s'", tolower(parent_names)), collapse = ","))
n_concepts <- function(parent_names) as.numeric(run_sql(sprintf("SELECT COUNT(*) AS n FROM (%s)", desc_sql(parent_names)))$n)
first_date <- function(parent_names){
  x <- run_sql(sprintf("SELECT person_id, MIN(DATE(condition_start_date)) AS first_dt FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id", cdr, desc_sql(parent_names)))
  x$person_id <- as.character(x$person_id); x$first_dt <- as.Date(x$first_dt); x
}
## flag cohort members whose first record of a condition is on or before baseline
pre_baseline <- function(fd, d = cox){ m <- match(d$person_id, fd$person_id); !is.na(m) & fd$first_dt[m] <= d$expo_date }
