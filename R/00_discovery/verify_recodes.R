# verify_recodes.R -- audit EVERY hard-coded answer-string recode in the pipeline against the
# answer strings actually present in ds_survey. Motivated by the education bug (educ_map in the
# build scripts matched only 2 of the 8 "Highest Grade:" answers, so 69% was coded missing).
# For each question concept used anywhere in run_all.R: prints the answer distribution, which
# strings the pipeline's map covers, and the share of non-skip answers that would silently become
# NA. Output is AGGREGATE ONLY (cells < 20 are printed as "<20"). Read-only; writes nothing.
# INPUT: ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint = "character")

## ---- maps exactly as they appear in the build / robustness scripts ------------------------------
map5      <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4      <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
map4_sw   <- c("Strongly disagree"=1,"Somewhat disagree"=2,"Somewhat agree"=3,"Strongly agree"=4)
income_map<- c("Annual Income: less 10k"=1,"Annual Income: 10k 25k"=2,"Annual Income: 25k 35k"=3,"Annual Income: 35k 50k"=4,"Annual Income: 50k 75k"=5,"Annual Income: 75k 100k"=6,"Annual Income: 100k 150k"=7,"Annual Income: 150k 200k"=8,"Annual Income: more 200k"=9)
educ_map_current <- c("Less than a high school degree or equivalent"=1,"Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"College graduate or advanced degree"=4)
educ_map_fixed   <- c("Highest Grade: Never Attended"=1,"Highest Grade: One Through Four"=1,"Highest Grade: Five Through Eight"=1,"Highest Grade: Nine Through Eleven"=1,
                      "Highest Grade: Twelve Or GED"=2,"Highest Grade: College One to Three"=3,"Highest Grade: College Graduate"=4,"Highest Grade: Advanced Degree"=4)
happy_map   <- c("Extremely happy"=6,"Very happy"=5,"Moderately happy"=4,"Moderately unhappy"=3,"Very unhappy"=2,"Extremely unhappy"=1)
meaning_map <- c("Not at all"=1,"A little"=2,"A moderate amount"=3,"Very much"=4,"An extreme amount"=5)
phq_map     <- c("Not at all"=0,"Several days"=1,"More than half the days"=2,"Nearly every day"=3)
pcl_map     <- c("Not at all"=0,"A little bit"=1,"Moderately"=2,"Quite a bit"=3,"Extremely"=4)
disc6       <- c("Never"=1,"Less than once a year"=2,"A few times a year"=3,"A few times a month"=4,"At least once a week"=5,"Almost everyday"=6)
food3       <- c("Never true"=1,"Sometimes true"=2,"Often true"=3)
# union of the na_vals lists used across scripts (any of these is an intended missing)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","PMI: Dont Know","Skip","Don't know",
             "Prefer not to answer","NA","Response removed due to invalid value")

## ---- which map each concept is run through -------------------------------------------------------
spec <- list(
  map5    = c(help=40192463, getalong=40192411, trust=40192499, values=40192417),
  map4    = c(watchout=40192400, graffiti=40192420, noisy=40192522, vandalism=40192412, abandoned=40192469,
              hangaround=40192500, crime=40192493, drug=40192457, alcohol=40192476, trouble=40192404,
              clean=40192456, upkeep=40192386, safe=40192384),
  map4_sw = c(unsafe_night=40192492, unsafe_day=40192414, shops=40192436, transit=40192440,
              sidewalks=40192437, bike=40192431, recreation=40192410),
  income_map       = c(income=1585375),
  educ_map_current = c(education=1585940),
  educ_map_fixed   = c(education=1585940),
  employment_raw   = c(employment=1585952),           # used as a raw factor; listed for completeness
  happy_map   = c(happy=1703980),
  meaning_map = c(meaning=1704001),
  pcl_map     = c(cutoff=1703997),
  phq_map     = c(down=1704024, anhedonia=1704026, worthless=1704041, fatigue=1704039, sleep=1703983,
                  appetite=1704004, concen=1704038, psychomotor=1703996, restless=1704028, suicidal=1703977,
                  nervous=1703984, worrystop=1703995, worrymuch=1704000, relax=1703987, afraid=1704042, irritable=1703920),
  disc6 = c(disc1=40192466, disc2=40192489, disc3=40192416, disc4=40192490, disc5=40192380,
            disc6=40192395, disc7=40192496, disc8=40192519, disc9=40192451),
  food3 = c(food1=40192517, food2=40192426),
  cnt_regex = c(disab1=903573, disab2=903574, disab3=903575, disab4=903576, disab5=903577, disab6=903578,
                ace1=1332946, ace2=1332947, ace3=1332948, ace4=1332950, ace5=1333208, ace6=1332857,
                ace7=1333210, ace8=1333212, ace9=1332951, ace10=1333202, ace11=1333203,
                tr1=1703988, tr2=1704021, tr3=1703976, tr4=1704055, tr5=1704035, tr6=1703994,
                tr7=1704056, tr8=1704011, tr9=1704029, tr10=1704020, tr11=1703989)
)
maps <- list(map5=map5, map4=map4, map4_sw=map4_sw, income_map=income_map,
             educ_map_current=educ_map_current, educ_map_fixed=educ_map_fixed,
             happy_map=happy_map, meaning_map=meaning_map, pcl_map=pcl_map, phq_map=phq_map,
             disc6=disc6, food3=food3)
cnt_pos <- "Yes|Once|More than once"; cnt_neg <- "^No$|: No$|^Never$"       # as in etable_robustness.R cnt()

## ---- pull the answer distributions (one query, aggregate) ---------------------------------------
all_ids <- unique(unlist(spec))
a <- run_sql(sprintf("SELECT question_concept_id, ANY_VALUE(question) AS question, answer, COUNT(*) AS n
                      FROM `%s.ds_survey` WHERE question_concept_id IN (%s)
                      GROUP BY question_concept_id, answer", cdr, paste(all_ids, collapse=",")))
a$n <- as.numeric(a$n); a$qid <- as.numeric(a$question_concept_id)
fmt_n <- function(n) ifelse(n < 20, "<20", format(n, big.mark=","))

## ---- audit ---------------------------------------------------------------------------------------
summary_rows <- list()
for (mname in names(spec)) {
  for (lab in names(spec[[mname]])) {
    id <- spec[[mname]][[lab]]
    sub <- a[a$qid == id, ]; sub <- sub[order(-sub$n), ]
    cat(sprintf("\n==== %s  (%s, concept %d)  map = %s ====\n", lab, substr(sub$question[1], 1, 70), id, mname))
    if (nrow(sub) == 0) { cat("  NO ROWS IN ds_survey\n"); next }
    if (mname == "cnt_regex") {
      status <- ifelse(grepl(cnt_pos, sub$answer), "pos=1",
                ifelse(grepl(cnt_neg, sub$answer), "neg=0",
                ifelse(sub$answer %in% na_vals | grepl("Prefer Not To Answer", sub$answer), "skip", "UNMAPPED->NA")))
    } else if (mname == "employment_raw") {
      status <- ifelse(sub$answer %in% na_vals, "skip", "level")
    } else {
      status <- ifelse(sub$answer %in% names(maps[[mname]]), "mapped",
                ifelse(sub$answer %in% na_vals, "skip", "UNMAPPED->NA"))
    }
    print(data.frame(answer = sub$answer, n = fmt_n(sub$n), status = status), row.names = FALSE)
    substantive <- sub$n[status != "skip"]
    lost <- sum(sub$n[status == "UNMAPPED->NA"])
    pct <- if (sum(substantive) > 0) 100 * lost / sum(substantive) else NA
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      item = lab, concept = id, map = mname, n_rows = sum(sub$n),
      pct_lost = round(pct, 1), flag = ifelse(!is.na(pct) & pct > 0, "CHECK", "ok"))
  }
}
S <- do.call(rbind, summary_rows)
cat("\n\n######## SUMMARY: share of substantive (non-skip) answers the pipeline silently sets to NA ########\n")
print(S[order(-S$pct_lost), ], row.names = FALSE)
cat("\nAny row with pct_lost > 0 is a recode-map mismatch of the education kind and must be fixed before rerunning.\n")
