# build_wb_items.R -- item-level well-being/affect cohort for the combined paper.
# Pulls ALL individual environment items (SDOH) + a curated EHHWB affect item set,
# scores every item oriented "higher = better mental state", and attaches the
# prospective incident-depression survival outcome (baseline = EHHWB happy-item date).
library(tidyverse); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## ---- item dictionaries -------------------------------------------------------
env_ids <- c(help=40192463, getalong=40192411, trust=40192499, values=40192417, watchout=40192400,
             graffiti=40192420, noisy=40192522, vandalism=40192412, abandoned=40192469,
             hangaround=40192500, crime=40192493, drug=40192457, alcohol=40192476, trouble=40192404,
             clean=40192456, upkeep=40192386, safe=40192384, unsafe_night=40192492, unsafe_day=40192414,
             shops=40192436, transit=40192440, sidewalks=40192437, bike=40192431, recreation=40192410)
# affect items, tagged by construct-overlap tier with the depression Dx outcome
wb_ids  <- c(happy=1703980, meaning=1704001,                                   # positive wellbeing (LOW)
             cutoff=1703997,                                                   # social-affective (MOD)
             down=1704024, anhedonia=1704026, worthless=1704041, fatigue=1704039,
             sleep=1703983, appetite=1704004, concen=1704038, psychomotor=1703996,
             restless=1704028, suicidal=1703977,                              # PHQ-9 (HIGH)
             nervous=1703984, worrystop=1703995, worrymuch=1704000, relax=1703987,
             afraid=1704042, irritable=1703920)                               # GAD-7 (MOD-HIGH)
ses_ids <- c(income=1585375, education=1585940)
tier <- c(happy="wellbeing", meaning="wellbeing", cutoff="social",
          down="phq", anhedonia="phq", worthless="phq", fatigue="phq", sleep="phq",
          appetite="phq", concen="phq", psychomotor="phq", restless="phq", suicidal="phq",
          nervous="gad", worrystop="gad", worrymuch="gad", relax="gad", afraid="gad", irritable="gad")

all_ids <- c(env_ids, wb_ids, ses_ids)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey` WHERE question_concept_id IN (%s)",
                      cdr, paste(all_ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(all_ids)[match(as.numeric(sv$question_concept_id), all_ids)]
wide <- tidyr::pivot_wider(sv[!is.na(sv$item), c("person_id","item","answer")],
                           names_from=item, values_from=answer)

demo <- run_sql(sprintf(paste("SELECT p.person_id, DATE_DIFF(CURRENT_DATE, DATE(p.birth_datetime), YEAR) AS age,",
  "g.concept_name AS sex, r.concept_name AS race, e.concept_name AS ethnicity FROM `%s.person` p",
  "LEFT JOIN `%s.concept` g ON p.gender_concept_id=g.concept_id",
  "LEFT JOIN `%s.concept` r ON p.race_concept_id=r.concept_id",
  "LEFT JOIN `%s.concept` e ON p.ethnicity_concept_id=e.concept_id"), cdr,cdr,cdr,cdr))
demo$age <- as.numeric(demo$age)
d <- merge(wide, demo, by="person_id")

## ---- recode maps -------------------------------------------------------------
na_vals <- c("Don't know/Not sure","Does not apply to my neighborhood","PMI: Skip",
             "PMI: Prefer Not To Answer","Skip","Don't know","Prefer not to answer",
             "Response removed due to invalid value")
rc <- function(x,m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
map5_agree <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4_agree <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
map4_sw    <- c("Strongly disagree"=1,"Somewhat disagree"=2,"Somewhat agree"=3,"Strongly agree"=4)
phq_map    <- c("Not at all"=0,"Several days"=1,"More than half the days"=2,"Nearly every day"=3)
pcl_map    <- c("Not at all"=0,"A little bit"=1,"Moderately"=2,"Quite a bit"=3,"Extremely"=4)
happy_map  <- c("Extremely happy"=6,"Very happy"=5,"Moderately happy"=4,"Moderately unhappy"=3,"Very unhappy"=2,"Extremely unhappy"=1)
meaning_map<- c("Not at all"=1,"A little"=2,"A moderate amount"=3,"Very much"=4,"An extreme amount"=5)

## ---- environment items: orient higher = better, keep raw numeric -------------
for (nm in c("help","getalong","trust","values")) d[[nm]] <- rc(d[[nm]], map5_agree)
d$watchout <- rc(d$watchout, map4_agree)
for (nm in c("graffiti","noisy","vandalism","abandoned","hangaround","crime","drug","alcohol","trouble"))
  d[[nm]] <- 5 - rc(d[[nm]], map4_agree)                       # disorder present -> reverse
for (nm in c("clean","upkeep","safe")) d[[nm]] <- rc(d[[nm]], map4_agree)
for (nm in c("unsafe_night","unsafe_day")) d[[nm]] <- 5 - rc(d[[nm]], map4_sw)
for (nm in c("shops","transit","sidewalks","bike","recreation")) d[[nm]] <- rc(d[[nm]], map4_sw)
# cohesion composite (5 items: z each, mean, re-z) for the mediation/bridge tests
Zc <- scale(as.matrix(d[c("help","getalong","trust","values","watchout")]))
coh <- rowMeans(Zc, na.rm=TRUE); coh[rowMeans(!is.na(Zc)) < 0.6] <- NA
d$z_cohesion <- scale(coh)[,1]

## ---- affect items: orient higher = better mental state (reverse distress) ----
d$happy_n   <- rc(d$happy,   happy_map)
d$meaning_n <- rc(d$meaning, meaning_map)
d$cutoff_n  <- rc(d$cutoff,  pcl_map)
phq_gad <- names(tier)[tier %in% c("phq","gad")]
for (nm in phq_gad) d[[paste0(nm,"_n")]] <- rc(d[[nm]], phq_map)
# z-scores, all oriented so higher = better wellbeing / fewer symptoms
d$z_happy   <- scale(d$happy_n)[,1]
d$z_meaning <- scale(d$meaning_n)[,1]
d$z_cutoff  <- -scale(d$cutoff_n)[,1]
for (nm in phq_gad) d[[paste0("z_",nm)]] <- -scale(d[[paste0(nm,"_n")]])[,1]
wb_names <- names(tier)                                        # happy, meaning, cutoff, phq..., gad...

## ---- SES numeric -------------------------------------------------------------
income_map <- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
educ_map <- c("Less than a high school degree or equivalent"=1,"Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"College graduate or advanced degree"=4)
d$income_n <- rc(d$income, income_map); d$educ_n <- rc(d$education, educ_map)

## ---- prospective incident-depression outcome (baseline = happy-item date) ----
dep_ids <- sprintf("SELECT a.descendant_concept_id FROM `%s.concept` p JOIN `%s.concept_ancestor` a ON a.ancestor_concept_id=p.concept_id JOIN `%s.concept` c ON c.concept_id=a.descendant_concept_id WHERE p.vocabulary_id='SNOMED' AND p.standard_concept='S' AND p.concept_name IN ('Depressive disorder','Major depressive disorder','Recurrent depressive disorder','Dysthymia') AND c.standard_concept='S' AND c.domain_id='Condition' AND LOWER(c.concept_name) NOT LIKE '%%bipolar%%'", cdr,cdr,cdr)
base <- run_sql(sprintf("
WITH bl AS (SELECT person_id, MIN(DATE(survey_datetime)) AS base_date FROM `%s.ds_survey` WHERE question_concept_id=1703980 GROUP BY person_id),
dep AS (SELECT person_id, MIN(condition_start_date) AS first_dep FROM `%s.condition_occurrence` WHERE condition_concept_id IN (%s) GROUP BY person_id),
ehr AS (SELECT person_id, MIN(condition_start_date) AS first_ehr, MAX(condition_start_date) AS last_ehr, COUNT(*) AS n_cond FROM `%s.condition_occurrence` GROUP BY person_id)
SELECT bl.person_id, bl.base_date, dep.first_dep, ehr.last_ehr, ehr.n_cond
FROM bl JOIN ehr USING(person_id) LEFT JOIN dep USING(person_id)
WHERE ehr.last_ehr >= bl.base_date AND ehr.first_ehr <= DATE_SUB(bl.base_date, INTERVAL 365 DAY)
  AND (dep.first_dep IS NULL OR dep.first_dep > bl.base_date)
", cdr, cdr, dep_ids, cdr))
base$base_date <- as.Date(base$base_date); base$first_dep <- as.Date(base$first_dep); base$last_ehr <- as.Date(base$last_ehr)
base$n_cond <- as.numeric(base$n_cond)
base$event <- as.integer(!is.na(base$first_dep))
base$time_days <- ifelse(base$event==1, as.numeric(base$first_dep - base$base_date),
                                        as.numeric(base$last_ehr  - base$base_date))
base <- base[base$time_days > 0, ]

wb_item_dat <- merge(d, base[c("person_id","event","time_days","n_cond")], by="person_id", all.x=TRUE)
wb_item_dat$log_util <- log1p(wb_item_dat$n_cond)
attr(wb_item_dat, "wb_names") <- wb_names
attr(wb_item_dat, "tier")     <- tier
saveRDS(wb_item_dat, "wb_item_dat.rds")
cat("wb_item_dat: rows =", nrow(wb_item_dat),
    "| at-risk (has survival) =", sum(!is.na(wb_item_dat$event)),
    "| incident events =", sum(wb_item_dat$event, na.rm=TRUE), "\n")
cat("affect-item coverage (non-NA z):\n")
print(sapply(paste0("z_", wb_names), function(v) sum(!is.na(wb_item_dat[[v]]))))
