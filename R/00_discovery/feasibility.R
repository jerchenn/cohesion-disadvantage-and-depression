# feasibility.R -- decides whether the cross-modal design is viable BEFORE we build it.
# Reports, aggregate only:
#   (1) cox_dat.rds structure -- so the Cox re-run of block D can be wired to the right columns
#   (2) cohort sizes and overlaps: EHHWB symptom responders, SDOH responders, EHR-cohort (cox_dat),
#       and Fitbit participants with adequate wear -- the load-bearing n for the device arm
#   (3) which Fitbit/activity tables exist in this CDR (names vary), so the wear query targets the right one
#
# INPUT: cox_dat.rds + ds_survey + CDR fitbit tables.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
dataset <- sub("\\..*$", "", cdr)   # project.dataset -> dataset (for INFORMATION_SCHEMA)

## ---------- (1) cox_dat structure ----------
cox <- readRDS("cox_dat.rds")
cat("========== (1) cox_dat.rds ==========\n")
cat(sprintf("rows: %d | cols: %d\n", nrow(cox), ncol(cox)))
cat("columns:\n"); print(names(cox))
cat("\nclasses:\n"); print(sapply(cox, function(x) class(x)[1]))
## try to surface the event/time columns for the Cox
ev <- names(cox)[sapply(cox, function(x) is.logical(x) ||
        (is.numeric(x) && all(na.omit(unique(x)) %in% c(0,1))))]
tm <- names(cox)[sapply(cox, function(x) is.numeric(x) && any(na.omit(x) > 30))]
cat(sprintf("\ncandidate event (0/1) cols: %s\n", paste(ev, collapse=", ")))
cat(sprintf("candidate time cols (numeric, values>30): %s\n", paste(tm, collapse=", ")))
for (c in intersect(ev, names(cox)))
  cat(sprintf("  %-22s events(1)=%d  n=%d\n", c, sum(cox[[c]]==1,na.rm=TRUE), sum(!is.na(cox[[c]]))))

## ---------- (2) cohort sizes ----------
sx_ids  <- c(1704026,1704024,1703983,1704039,1704004,1704041,1704038,1703996,1703977,
             1703984,1703995,1704000,1703987,1704028,1703920,1704042)   # 16 EHHWB PHQ/GAD items
sdoh_probe <- c(40192463,40192499,40192384,40192517,40192507)           # a few SDOH anchors
n1 <- run_sql(sprintf("SELECT COUNT(DISTINCT person_id) n FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(sx_ids, collapse=",")))$n
n2 <- run_sql(sprintf("SELECT COUNT(DISTINCT person_id) n FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(sdoh_probe, collapse=",")))$n
n12 <- run_sql(sprintf(
  "SELECT COUNT(*) n FROM (
     SELECT person_id FROM `%s.ds_survey` WHERE question_concept_id IN (%s) GROUP BY person_id
     INTERSECT DISTINCT
     SELECT person_id FROM `%s.ds_survey` WHERE question_concept_id IN (%s) GROUP BY person_id)",
  cdr, paste(sx_ids, collapse=","), cdr, paste(sdoh_probe, collapse=",")))$n
cat("\n========== (2) cohort sizes ==========\n")
cat(sprintf("EHHWB symptom responders (any of 16): %s\n", n1))
cat(sprintf("SDOH responders (any of 5 anchors):   %s\n", n2))
cat(sprintf("EHHWB AND SDOH:                        %s\n", n12))
cat(sprintf("EHR cohort (cox_dat rows):             %d\n", nrow(cox)))

## ---------- (3) fitbit tables + wear ----------
cat("\n========== (3) Fitbit/activity tables in CDR ==========\n")
tabs <- run_sql(sprintf(
  "SELECT table_name FROM `%s.%s.INFORMATION_SCHEMA.TABLES`
   WHERE LOWER(table_name) LIKE '%%activit%%' OR LOWER(table_name) LIKE '%%fitbit%%'
      OR LOWER(table_name) LIKE '%%step%%'  OR LOWER(table_name) LIKE '%%heart%%'
      OR LOWER(table_name) LIKE '%%sleep%%' OR LOWER(table_name) LIKE '%%device%%'
   ORDER BY table_name", proj, dataset))
print(tabs)
## if an activity_summary table exists, count participants by wear-days
if ("activity_summary" %in% tabs$table_name){
  w <- run_sql(sprintf(
    "WITH d AS (SELECT person_id, COUNT(DISTINCT date) days FROM `%s.activity_summary`
                WHERE steps IS NOT NULL GROUP BY person_id)
     SELECT SUM(days>=1) any_wear, SUM(days>=30) w30, SUM(days>=90) w90 FROM d", cdr))
  cat(sprintf("\nactivity_summary wear: any=%s | >=30 days=%s | >=90 days=%s\n",
              w$any_wear, w$w30, w$w90))
  ov <- run_sql(sprintf(
    "WITH fit AS (SELECT person_id FROM `%s.activity_summary` WHERE steps IS NOT NULL
                  GROUP BY person_id HAVING COUNT(DISTINCT date) >= 30),
          ehh AS (SELECT person_id FROM `%s.ds_survey` WHERE question_concept_id IN (%s) GROUP BY person_id),
          sdo AS (SELECT person_id FROM `%s.ds_survey` WHERE question_concept_id IN (%s) GROUP BY person_id)
     SELECT COUNT(*) n FROM (SELECT person_id FROM fit
       INTERSECT DISTINCT SELECT person_id FROM ehh
       INTERSECT DISTINCT SELECT person_id FROM sdo)",
    cdr, cdr, paste(sx_ids, collapse=","), cdr, paste(sdoh_probe, collapse=",")))$n
  cat(sprintf("Fitbit(>=30d) AND EHHWB AND SDOH:      %s   <-- device-arm n\n", ov))
} else {
  cat("\nNo 'activity_summary' table by that name -- pick the right one from the list above and I'll target it.\n")
}
cat("\nDecision rule: device-arm n under ~2000 -> Fitbit is a sensitivity analysis, not a pillar.\n")
