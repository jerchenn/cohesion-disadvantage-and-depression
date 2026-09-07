# ehhwb_timing.R -- BLOCKING check. build_cox.R anchors the risk window at expo_date = the SDOH cohesion
# item's survey date. But ace/trauma (and the symptom outcomes) are EHHWB items on a separate survey date.
# If EHHWB postdates expo_date -- or postdates first_dep -- then ACE endorsement can be inflated by having
# already been diagnosed (mood/diagnosis-congruent recall), inflating the coefficient that carries the thesis.
# This reports, per person in the analytic-style cohort:
#   expo_date (SDOH), ehhwb_date (EHHWB), first_dep; the gap distribution; and the share of EHHWB-after-window
#   and EHHWB-after-diagnosis cases. Aggregate only.
#
# INPUT: ds_survey + condition_occurrence + concept.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

ehhwb_ids <- c(1332946,1332947,1332948,1332950,1333208,1332857,1333210,1333212,1332951,1333202,1333203,  # ACE
               1703988,1704021,1703976,1704055,1704035,1703994,1704056,1704011,1704029,1704020,1703989,  # trauma
               1704026,1704024,1703983,1704039,1704004,1704041,1704038,1703996,1703977,                 # PHQ
               1703984,1703995,1704000,1703987,1704028,1703920,1704042)                                  # GAD
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)

q <- sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey`
              WHERE question_concept_id=40192463 GROUP BY person_id),
     ehh AS (SELECT person_id, MIN(DATE(survey_datetime)) AS ehhwb_date FROM `%s.ds_survey`
              WHERE question_concept_id IN (%s) GROUP BY person_id),
     dep AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence`
              WHERE condition_concept_id IN (%s) GROUP BY person_id)
SELECT expo.person_id, expo.expo_date, ehh.ehhwb_date, dep.first_dep
FROM expo JOIN ehh USING(person_id) LEFT JOIN dep USING(person_id)",
  cdr, cdr, paste(ehhwb_ids, collapse=","), cdr, dep_ids)
d <- run_sql(q)
d$expo_date <- as.Date(d$expo_date); d$ehhwb_date <- as.Date(d$ehhwb_date); d$first_dep <- as.Date(d$first_dep)
d$gap <- as.numeric(d$ehhwb_date - d$expo_date)   # + = EHHWB after the SDOH-anchored window start

cat(sprintf("People with both surveys: %d\n\n", nrow(d)))
cat("=== EHHWB date minus expo_date (days); + = EHHWB AFTER window start ===\n")
print(summary(d$gap))
cat(sprintf("\nEHHWB same day as expo (|gap|<=1): %.1f%%\n", 100*mean(abs(d$gap)<=1, na.rm=TRUE)))
cat(sprintf("EHHWB after window start (gap>1):   %.1f%%\n", 100*mean(d$gap>1, na.rm=TRUE)))
cat(sprintf("EHHWB before window start (gap<-1): %.1f%%\n", 100*mean(d$gap< -1, na.rm=TRUE)))

## the dangerous cases: incident depression, and EHHWB (where ACE is reported) completed AFTER the diagnosis
ev <- d[!is.na(d$first_dep) & d$first_dep > d$expo_date, ]   # incident cases (matches build_cox rule)
cat(sprintf("\n=== incident cases (first_dep after expo_date): %d ===\n", nrow(ev)))
cat(sprintf("  ACE reported AFTER diagnosis (ehhwb_date > first_dep): %.1f%%\n",
            100*mean(ev$ehhwb_date > ev$first_dep, na.rm=TRUE)))
cat(sprintf("  ACE reported BEFORE diagnosis (ehhwb_date <= first_dep): %.1f%%\n",
            100*mean(ev$ehhwb_date <= ev$first_dep, na.rm=TRUE)))

## clean subset for a re-run: EHHWB recorded at/before the window start
clean <- d[!is.na(d$gap) & d$gap <= 1, ]
cat(sprintf("\nEHHWB-precedes-window clean subset: %d (%.0f%% of both-survey cohort)\n",
            nrow(clean), 100*nrow(clean)/nrow(d)))
cat("\nDecision: if most EHHWB postdates expo_date, re-anchor expo_date at MAX(SDOH, EHHWB) date OR restrict\n",
    "to gap<=1, and re-run ACE. ACE HR surviving that restriction is the number the paper needs.\n", sep="")
