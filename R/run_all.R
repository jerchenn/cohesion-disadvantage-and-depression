# run_all.R -- master reproducibility script. Regenerates every reported estimate from a single set of
# versioned cached datasets, in dependency order. Run inside the All of Us Workbench from the repo root.
#
# DEPENDENCY GRAPH (-> = "produces", indent = "reads"):
#
#   build_cox3.R      -> cox_dat3.rds      PRIMARY cohort (baseline age, pre-baseline utilization, ZIP3)
#   build_cox2.R      -> cox_dat2.rds      re-anchored cohort for the adversity-exposome robustness (eTable 5)
#                          reads cox_dat3.rds
#   build_cox_swb.R   -> swb_dat.rds       self-reported wellbeing subcohort   (baseline age from cox_dat3)
#   build_wb_items.R  -> wb_item_dat.rds   item-level affect subcohort         (baseline age from cox_dat3)
#   build_fitbit.R    -> fitbit_dat.rds    wearable subcohort                  (reads cox_dat3.rds)
#   build_antidep.R   -> antidep_dat.rds   antidepressant parallel outcome     (reads cox_dat3.rds)
#
#   Models / robustness / outputs below all read the cached .rds above; none re-queries the cohort.
#   Every primary/Table 2/interaction/absolute-risk estimate reads cox_dat3.rds.
#   Secondary (wellbeing/item/wearable) read the subcohort .rds, which derive covariates from cox_dat3.
#
# NOTE: cox_dat.rds (from the older build_cox.R) is NOT used by any reported analysis; it is retained only
#       for historical reference. The single analysis cohort for the manuscript is cox_dat3.rds.

run <- function(f){ cat("\n>>>", f, "\n"); source(f, local = new.env()) }

## 1. datasets (order matters: cox_dat3 first; subcohorts read it)
run("R/01_build/build_cox3.R")
run("R/01_build/build_cox2.R")
run("R/01_build/build_cox_swb.R")
run("R/01_build/build_wb_items.R")
run("R/01_build/build_fitbit.R")
run("R/01_build/build_antidep.R")

## 2. primary + robustness (Table 2, eTables 4-7)
run("R/02_models/cox_model3.R")
run("R/03_robustness/cox_interaction3.R")
run("R/04_outputs/cox_absrisk3.R")
run("R/03_robustness/rake_ipw.R")
run("R/03_robustness/baseline_history.R")
run("R/03_robustness/mi_sensitivity.R")
run("R/03_robustness/cox_mdd3.R")
run("R/03_robustness/evalue.R")
run("R/03_robustness/etable_robustness.R")
run("R/03_robustness/exposome_adjustment_bySES.R")

## 3. secondary + mechanistic (eTables 2-3, Figure 3, mediation)
run("R/02_models/cox_swb_model.R")
run("R/02_models/wb_items_cox.R")
run("R/02_models/fitbit_model.R")
run("R/02_models/sleep_headtohead.R")
run("R/02_models/env_wb_linkage.R")
run("R/03_robustness/mediation.R")

## 4. tables and figures
run("R/04_outputs/table1.R")
run("R/04_outputs/table2.R")
run("R/04_outputs/figs_data.R")
run("R/04_outputs/figs_plot.R")
run("R/04_outputs/fig_stratified.R")
run("R/04_outputs/fig_items.R")
run("R/04_outputs/fig_doseresponse.R")
run("R/04_outputs/fig_gradient.R")

cat("\n=== run_all complete: every reported estimate regenerated from cox_dat3.rds (+ subcohorts) ===\n")
