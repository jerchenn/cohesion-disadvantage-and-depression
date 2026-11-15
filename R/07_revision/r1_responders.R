# r1_responders.R -- Reviewer 1, first point: how do Social Determinants of Health survey responders differ
# from the broader EHR-linked All of Us population, and does weighting for survey response change the estimate?
#  (a) demographic comparison: responders to the cohesion item vs all participants with EHR condition data
#  (b) response propensity (logistic: age band, gender, race/ethnicity, EHR record length) -> inverse probability
#      of response weights for the at-risk cohort -> primary refitted with the weights (robust SE)
# Person-level rows stay inside the Workbench; only percentages, masked counts and model estimates are printed.
# Needs: survival, bigrquery.
source("R/07_revision/_common.R", local = TRUE); bq()

pop <- run_sql(sprintf("
WITH ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond
             FROM `%s.condition_occurrence` GROUP BY person_id),
     resp AS (SELECT DISTINCT person_id FROM `%s.ds_survey` WHERE question_concept_id = 40192463)
SELECT p.person_id, EXTRACT(YEAR FROM p.birth_datetime) AS birth_year, g.concept_name AS gender, r.concept_name AS race, e.concept_name AS ethnicity,
       ehr.first_ehr, ehr.last_ehr, ehr.n_cond, IF(resp.person_id IS NULL, 0, 1) AS responder
FROM `%s.person` p JOIN ehr USING(person_id)
LEFT JOIN resp USING(person_id)
LEFT JOIN `%s.concept` g ON p.gender_concept_id = g.concept_id
LEFT JOIN `%s.concept` r ON p.race_concept_id = r.concept_id
LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id = e.concept_id", cdr, cdr, cdr, cdr, cdr, cdr))
pop$person_id <- as.character(pop$person_id); pop$responder <- as.integer(pop$responder)
pop$age <- 2023 - as.numeric(pop$birth_year)
pop$age_band <- cut(pop$age, c(-Inf, 29, 44, 64, Inf), labels = c("18-29","30-44","45-64","65+"))
pop$gender_c <- ifelse(pop$gender %in% c("Female","Male"), pop$gender, "Other/not reported")
pop$race_eth <- with(pop, ifelse(ethnicity == "Hispanic or Latino", "Hispanic",
                     ifelse(race == "White", "White", ifelse(race == "Black or African American", "Black", ifelse(race == "Asian", "Asian", "Other/not reported")))))
pop$ehr_years <- as.numeric(as.Date(pop$last_ehr) - as.Date(pop$first_ehr)) / 365.25
pop$log_ncond <- log1p(as.numeric(pop$n_cond))

hdr("(a) Responders to the cohesion item vs all EHR-linked participants")
cat(sprintf("EHR-linked participants: %s; responders among them: %s (%.1f%%)\n", mask(nrow(pop)), mask(sum(pop$responder)), 100*mean(pop$responder)))
cmp <- function(v){ t <- prop.table(table(pop[[v]], pop$responder), 2) * 100
  out <- data.frame(level = rownames(t), non_responders_pct = sprintf("%.1f", t[, "0"]), responders_pct = sprintf("%.1f", t[, "1"]))
  cat(sprintf("\n-- %s --\n", v)); print(out, row.names = FALSE) }
cmp("age_band"); cmp("gender_c"); cmp("race_eth")
cat(sprintf("\nEHR record length, years: median %.1f non-responders vs %.1f responders; condition records: median %.0f vs %.0f\n",
            median(pop$ehr_years[pop$responder == 0], na.rm = TRUE), median(pop$ehr_years[pop$responder == 1], na.rm = TRUE),
            median(as.numeric(pop$n_cond[pop$responder == 0])), median(as.numeric(pop$n_cond[pop$responder == 1]))))

hdr("(b) Response propensity weights and the weighted primary model")
pm <- glm(responder ~ age_band + gender_c + race_eth + ehr_years + log_ncond, data = pop, family = binomial)
cf <- summary(pm)$coefficients
print(data.frame(term = rownames(cf), OR = sprintf("%.2f", exp(cf[, 1])), p = signif(cf[, 4], 2)), row.names = FALSE)
pop$p_resp <- fitted(pm)
cat(sprintf("response propensity among responders: median %.3f, IQR %.3f-%.3f; pseudo-R2 (McFadden) %.3f\n",
            median(pop$p_resp[pop$responder == 1]), quantile(pop$p_resp[pop$responder == 1], .25), quantile(pop$p_resp[pop$responder == 1], .75),
            1 - pm$deviance / pm$null.deviance))
d <- cox; d$w_resp <- 1 / pop$p_resp[match(d$person_id, pop$person_id)]
cat(sprintf("cohort members with a weight: %s of %s\n", mask(sum(!is.na(d$w_resp))), mask(nrow(d))))
d <- d[!is.na(d$w_resp), ]; d$w_resp <- d$w_resp / mean(d$w_resp)
cap <- quantile(d$w_resp, .99); d$w_resp[d$w_resp > cap] <- cap
cat(sprintf("normalised weights: median %.2f, 95th pct %.2f, truncated at %.2f; effective n %s\n", median(d$w_resp), quantile(d$w_resp, .95), cap, mask(round(sum(d$w_resp)^2 / sum(d$w_resp^2)))))
hr(fit(RHS, d), lab = "primary, unweighted (same rows)")
hr(coxph(as.formula(paste("Surv(time_days, event) ~", RHS)), d, weights = w_resp, robust = TRUE), lab = "primary, weighted for survey response (robust SE)")
cat("\nNote for the text: the weights correct for differential response by age, gender, race and ethnicity and EHR depth among EHR-linked participants; they cannot correct for selection into All of Us itself.\n")
