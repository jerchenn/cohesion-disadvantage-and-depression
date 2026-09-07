# build_cox2.R -- re-anchored survival dataset fixing the EHHWB timing problem.
# build_cox.R anchored the window at the SDOH cohesion date (expo_date); but ace/trauma/symptoms are EHHWB,
# completed a median 263 days later (47% of incident ACE reports postdate the diagnosis). Re-anchor at
# anchor = MAX(SDOH date, EHHWB date) so EVERY exposure precedes follow-up, and recompute event/time from
# anchor with prevalent-at-anchor excluded. Reuses all covariates from cox_dat.rds; only the survival
# columns change. Saves cox_dat2.rds with the same schema so exposome_cox.R can read it directly.
#
# INPUT: cox_dat.rds + ds_survey + condition_occurrence + concept.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

ehhwb_ids <- c(1332946,1332947,1332948,1332950,1333208,1332857,1333210,1333212,1332951,1333202,1333203,
               1703988,1704021,1703976,1704055,1704035,1703994,1704056,1704011,1704029,1704020,1703989,
               1704026,1704024,1703983,1704039,1704004,1704041,1704038,1703996,1703977,
               1703984,1703995,1704000,1703987,1704028,1703920,1704042)
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)

surv <- run_sql(sprintf("
WITH expo AS (SELECT person_id, MIN(DATE(survey_datetime)) AS expo_date FROM `%s.ds_survey` WHERE question_concept_id=40192463 GROUP BY person_id),
     ehh  AS (SELECT person_id, MIN(DATE(survey_datetime)) AS ehhwb_date FROM `%s.ds_survey` WHERE question_concept_id IN (%s) GROUP BY person_id),
     anc  AS (SELECT person_id, GREATEST(expo_date, ehhwb_date) AS anchor FROM expo JOIN ehh USING(person_id)),
     dep  AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
     ehr  AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT anc.person_id, anc.anchor, dep.first_dep, ehr.last_ehr, ehr.n_cond
FROM anc JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)
WHERE ehr.last_ehr >= anc.anchor AND ehr.first_ehr <= DATE_SUB(anc.anchor, INTERVAL 365 DAY)
  AND (dep.first_dep IS NULL OR dep.first_dep > anc.anchor)
", cdr, cdr, paste(ehhwb_ids, collapse=","), cdr, dep_ids, cdr))
surv$anchor <- as.Date(surv$anchor); surv$first_dep <- as.Date(surv$first_dep); surv$last_ehr <- as.Date(surv$last_ehr)
surv$n_cond <- as.numeric(surv$n_cond)
surv$event <- as.integer(!is.na(surv$first_dep))
surv$time_days <- ifelse(surv$event==1, as.numeric(surv$first_dep - surv$anchor),
                                         as.numeric(surv$last_ehr  - surv$anchor))
surv <- surv[surv$time_days > 0, ]

## reuse covariates from cox_dat; swap in the re-anchored survival + refreshed utilization
base <- cox[, setdiff(names(cox), c("event","time_days","n_cond","log_util"))]
cox2 <- merge(base, surv[c("person_id","event","time_days","n_cond")], by="person_id")
cox2$log_util <- log1p(cox2$n_cond)
saveRDS(cox2, "cox_dat2.rds")
cat(sprintf("cox_dat2 (re-anchored at MAX(SDOH,EHHWB)): n=%d  events=%d  median FU(yr)=%.2f\n",
            nrow(cox2), sum(cox2$event), median(cox2$time_days)/365.25))
cat(sprintf("  (was: n=%d events=%d under SDOH-only anchor)\n", nrow(cox), sum(cox$event)))
