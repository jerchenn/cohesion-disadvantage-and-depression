# run_all.R -- replication script for "Neighborhood Social Cohesion, Socioeconomic Disadvantage, and Incident
# Depression" (JAMA Network Open JNO26-18095, revision 1). Rebuilds every analysis dataset from the All of Us
# Controlled Tier (C2025Q4R6) and regenerates every estimate in the manuscript and Supplement 1, in dependency
# order. Run from a fresh R notebook in an All of Us Workbench workspace:
#
#   if (!dir.exists("~/cohesion-disadvantage-and-depression"))
#     system("git clone https://github.com/jerchenn/cohesion-disadvantage-and-depression.git ~/cohesion-disadvantage-and-depression")
#   setwd("~/cohesion-disadvantage-and-depression"); system("git pull")
#   source("R/run_all.R")
#
# Each script runs in its own environment. A script that fails is logged and the run continues; the list of
# failures prints at the end. To resume after a failure, set START_AT to that script's path and source again,
# eg START_AT <- "R/07_revision/r2_ipcw.R"; source("R/run_all.R"). The full run takes several hours; the
# BigQuery pulls and the 20-imputation models are the slow steps. All printed output is aggregate (counts under
# 20 masked) and is copied to replication_log.txt in the repo root.

## 0. packages
need <- c("survival", "bigrquery", "dplyr", "tidyr", "tidyverse", "readr", "ggplot2", "sandwich", "lmtest", "gridExtra", "EValue")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss)
stopifnot(nzchar(Sys.getenv("WORKSPACE_CDR")), nzchar(Sys.getenv("GOOGLE_PROJECT")))

steps <- c(
  ## 1. analysis datasets (cox_dat3 first; the others read it)
  "R/01_build/build_cox3.R",            # primary cohort -> cox_dat3.rds
  "R/01_build/build_cox2.R",            # re-anchored cohort for the adversity models -> cox_dat2.rds
  "R/01_build/build_cox_swb.R",         # self-reported wellbeing subcohort -> swb_dat.rds
  "R/01_build/build_wb_items.R",        # item-level affect subcohort -> wb_item_dat.rds
  "R/01_build/build_fitbit.R",          # wearable subcohort -> fitbit_dat.rds
  "R/01_build/build_antidep.R",         # antidepressant parallel outcome -> antidep_dat.rds
  ## 2. checks on the cohort: recodes, integrity, alpha, depression code list (eTable 2)
  "R/00_discovery/verify_recodes.R",
  "R/00_discovery/sanity_checks.R",
  ## 3. primary model and original sensitivity analyses (Table 2, eTables 3, 4, 10, 15)
  "R/02_models/cox_model3.R",
  "R/03_robustness/cox_interaction3.R",
  "R/04_outputs/cox_absrisk3.R",
  "R/03_robustness/rake_ipw.R",
  "R/03_robustness/baseline_history.R",
  "R/03_robustness/mi_sensitivity.R",
  "R/03_robustness/cox_mdd3.R",
  "R/03_robustness/evalue.R",
  "R/03_robustness/etable_robustness.R",
  "R/03_robustness/exposome_adjustment_bySES.R",
  ## 4. subcohort analyses (Figure 3, eTables 4, 13, 14, mediation)
  "R/02_models/cox_swb_model.R",
  "R/02_models/wb_items_cox.R",
  "R/02_models/fitbit_model.R",
  "R/02_models/sleep_headtohead.R",
  "R/02_models/env_wb_linkage.R",
  "R/03_robustness/mediation.R",
  ## 5. analyses added at revision (eTables 1, 5-9, 11, 12; Table 2 rows)
  "R/07_revision/r1_responders.R",
  "R/07_revision/r1_exclusions_lags_ehrlength.R",
  "R/07_revision/r2_ascertainment.R",
  "R/07_revision/r2_ipcw.R",
  "R/07_revision/r2_mi_primary.R",
  "R/07_revision/r2_multilevel.R",
  "R/07_revision/r2_reri_subsamples_vif.R",
  "R/07_revision/strata_income_deprivation.R",
  "R/07_revision/r3_raking_subset.R",
  "R/07_revision/r3_final_checks.R",
  ## 6. tables and figures
  "R/04_outputs/table1.R",
  "R/04_outputs/table2.R",
  "R/04_outputs/figs_data.R",           # -> fig_estimates.csv
  "R/04_outputs/fig_items.R",           # -> fig_items_data.csv
  "R/04_outputs/fig_stratified.R",      # Figure 1
  "R/04_outputs/fig_forest.R",          # Figures 2 and 3 (log-scale forest plots) -> figures/
  "R/04_outputs/fig_doseresponse.R",    # eFigure 2
  "R/04_outputs/fig_gradient.R"         # eFigure 3
)

resume <- exists("START_AT", envir = globalenv())
if (resume) {
  start <- get("START_AT", envir = globalenv()); rm("START_AT", envir = globalenv())   # used once, so a later full run starts from the beginning
  if (!file.exists("cox_dat3.rds")) stop("START_AT was set but this folder has no cox_dat3.rds; run the full script without START_AT")
  if (start %in% steps) steps <- steps[match(start, steps):length(steps)]
}
failed <- character(0)
logcon <- file("replication_log.txt", open = if (resume) "a" else "w")
sink(logcon, split = TRUE)
cat("Replication run started", format(Sys.time()), "|", R.version.string, "| CDR", Sys.getenv("WORKSPACE_CDR"), "\n")
for (f in steps) {
  cat("\n>>>", f, format(Sys.time(), "%H:%M"), "\n")
  ok <- tryCatch({ source(f, local = new.env()); TRUE },
                 error = function(e) { cat("!!! FAILED:", f, "--", conditionMessage(e), "\n"); FALSE })
  if (!ok) failed <- c(failed, f)
}
cat("\nReplication run finished", format(Sys.time()), "\n")
if (length(failed)) cat("Scripts that failed (rerun with START_AT):", paste(failed, collapse = ", "), "\n") else cat("All scripts completed.\n")
sink(); close(logcon)
