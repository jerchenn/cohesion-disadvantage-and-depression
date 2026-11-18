# r3_raking_subset.R -- Reviewer 2 point 4: compare the raking-complete subset (rake_ipw.R definition)
# with the full at-risk cohort and with the participants it leaves out. Aggregate only; counts under 20 masked.
source("R/07_revision/_common.R", local = TRUE)
cox$sex_r    <- ifelse(cox$sex %in% c("Female", "Male"), cox$sex, NA)
cox$educ_band <- c("1" = "<HS", "2" = "HS", "3" = "SomeCol", "4" = "BA+")[as.character(cox$educ_n)]
inrk <- !is.na(cox$age) & !is.na(cox$sex_r) & !is.na(cox$educ_band) & !is.na(cox$z_cohesion) & !is.na(cox$log_util)
re <- with(cox, ifelse(ethnicity == "Hispanic or Latino", "Hispanic", ifelse(race == "White", "White", ifelse(race == "Black or African American", "Black",
           ifelse(race == "Asian", "Asian", "Other or not reported")))))
pc <- function(x) { n <- sum(x, na.rm = TRUE); sprintf("%s (%.1f)", mask(n), 100 * n / length(x)) }
col <- function(i) { d <- cox[i, ]; r <- re[i]
  c(No. = mask(nrow(d)), `Incident depression` = pc(d$event == 1), `Age, mean (SD)` = sprintf("%.1f (%.1f)", mean(d$age), sd(d$age)),
    Women = pc(d$sex == "Female"), Asian = pc(r == "Asian"), Black = pc(r == "Black"), Hispanic = pc(r == "Hispanic"), White = pc(r == "White"),
    `Other or not reported` = pc(r == "Other or not reported"), `Income >=$75k` = pc(!is.na(d$income_n) & d$income_n >= 6),
    `Income not reported` = pc(is.na(d$income_n)), `College graduate or higher` = pc(!is.na(d$educ_n) & d$educ_n == 4),
    `Education not reported` = pc(is.na(d$educ_n)), `Cohesion score, mean (SD)` = sprintf("%.2f (%.2f)", mean(d$z_cohesion, na.rm = TRUE), sd(d$z_cohesion, na.rm = TRUE))) }
hdr("Raking-complete subset vs full at-risk cohort")
print(cbind(`Full at-risk cohort` = col(rep(TRUE, nrow(cox))), `Raking subset` = col(inrk), `Not in raking subset` = col(!inrk)), quote = FALSE)
cat("raking subset share of at-risk cohort:", sprintf("%.1f%%", 100 * mean(inrk)), "\n")
