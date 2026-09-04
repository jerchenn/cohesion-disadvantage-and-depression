# table1.R -- Table 1 (cohort characteristics by cohesion tertile) + eTable 1
# (subcohort comparison: main vs SWB vs Fitbit -> the selection-transparency table).
# JAMA format: No. (%) for categorical, mean (SD) or median (IQR) for continuous.
cox <- readRDS("cox_dat.rds")

race_eth <- function(d) with(d, ifelse(ethnicity=="Hispanic or Latino","Hispanic",
   ifelse(race=="White","White, non-Hispanic",
   ifelse(race=="Black or African American","Black, non-Hispanic",
   ifelse(race=="Asian","Asian","Other/Not reported")))))
inc_cat <- function(d) ifelse(is.na(d$income_n),"Not reported",
   ifelse(d$income_n<=3,"<$35k", ifelse(d$income_n<=5,"$35k-$75k","≥$75k")))

npct <- function(x,v) sprintf("%d (%.1f)", sum(x==v,na.rm=TRUE), 100*mean(x==v,na.rm=TRUE))
msd  <- function(x)   sprintf("%.1f (%.1f)", mean(x,na.rm=TRUE), sd(x,na.rm=TRUE))
miqr <- function(x){q<-quantile(x,c(.5,.25,.75),na.rm=TRUE); sprintf("%.1f (%.1f-%.1f)",q[1],q[2],q[3])}

summ <- function(d){
  d$re <- race_eth(d); d$ic <- inc_cat(d); fu <- d$time_days/365.25
  c(
   "No."                                = format(nrow(d), big.mark=" "),
   "Incident depression, No. (%)"       = npct(d$event,1),
   "Follow-up, median (IQR), y"         = miqr(fu),
   "Age, mean (SD), y"                  = msd(d$age),
   "Female, No. (%)"                    = npct(d$sex,"Female"),
   "  White, non-Hispanic"              = npct(d$re,"White, non-Hispanic"),
   "  Black, non-Hispanic"              = npct(d$re,"Black, non-Hispanic"),
   "  Hispanic"                         = npct(d$re,"Hispanic"),
   "  Asian"                            = npct(d$re,"Asian"),
   "  Other/Not reported"               = npct(d$re,"Other/Not reported"),
   "  Income <$35k"                     = npct(d$ic,"<$35k"),
   "  Income $35k-$75k"                 = npct(d$ic,"$35k-$75k"),
   "  Income ≥$75k"                = npct(d$ic,"≥$75k"),
   "  Income not reported"              = npct(d$ic,"Not reported"),
   "EHR conditions, median (IQR)"       = miqr(d$n_cond),
   "Cohesion score, mean (SD)"          = msd(d$z_cohesion))
}

## ---- Table 1: by cohesion tertile ----
cox$coh_t <- cut(cox$z_cohesion, quantile(cox$z_cohesion,0:3/3,na.rm=TRUE),
                 include.lowest=TRUE, labels=c("Low","Middle","High"))
t1 <- data.frame(Characteristic=names(summ(cox)),
                 Overall = summ(cox),
                 `Low cohesion`    = summ(cox[which(cox$coh_t=="Low"),]),
                 `Middle cohesion` = summ(cox[which(cox$coh_t=="Middle"),]),
                 `High cohesion`   = summ(cox[which(cox$coh_t=="High"),]),
                 check.names=FALSE)
write.csv(t1, "table1.csv", row.names=FALSE)
cat("=== TABLE 1 (by cohesion tertile) ===\n"); print(t1, row.names=FALSE)

## ---- eTable 1: subcohort comparison ----
fb  <- readRDS("fitbit_dat.rds"); fb <- fb[!is.na(fb$steps) & !is.na(fb$event),]
wb  <- readRDS("wb_item_dat.rds"); wb <- wb[!is.na(wb$event),]
e1 <- data.frame(Characteristic=names(summ(cox)),
                 `Main cohort`   = summ(cox),
                 `SWB subcohort` = summ(wb),
                 `Fitbit subcohort` = summ(fb),
                 check.names=FALSE)
write.csv(e1, "etable1.csv", row.names=FALSE)
cat("\n=== eTABLE 1 (subcohort comparison) ===\n"); print(e1, row.names=FALSE)
cat("\nWrote table1.csv and etable1.csv\n")
cat("NOTE: check the education variable separately; educ_n==4 (college graduate) appeared\n")
cat("absent in the cohesion cohort (response-string mapping issue) -> income used for SES here.\n")
