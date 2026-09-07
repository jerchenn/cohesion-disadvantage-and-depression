# exposome_cox.R -- the cross-modal test. Same exposures, same distress block, SAME sample, two outcomes:
#   (I)  self-reported symptom severity (EHHWB composite)  -- lm
#   (II) EHR-diagnosed incident depression (cox_dat)        -- Cox, Surv(time_days,event)
# For each: total (exposures+covs) vs direct (+distress block), and the share absorbed by distress.
# If a mediated share collapses when the outcome is clinical rather than self-reported, that localizes
# common-method inflation with a number.
# PLUS the ascertainment fork: do discrim / discrim_hc predict LOWER utilization (log_util, n_cond)?
#   -- decides whether a self-report/EHR gap means "distress not disease" or "under-diagnosis".
# Cross-sectional exposures, prospective EHR outcome. Aggregate only. Covariates age/sex/income (matched
# across both outcomes so the modality is the only difference).
#
# INPUT: cox_dat.rds + ds_survey.  Needs: bigrquery, survival.

library(bigrquery); library(survival)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

## ---------- scales + item spec (same recodes as exposure_model.R) ----------
S <- list(
 agree5 = c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5),
 agree4 = c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4),
 sw4    = c("Strongly disagree"=1,"Somewhat disagree"=2,"Somewhat agree"=3,"Strongly agree"=4),
 time5  = c("None of the time"=1,"A little of the time"=2,"Some of the time"=3,"Most of the time"=4,"All of the time"=5),
 lonely4= c("Never"=1,"Rarely"=2,"Sometimes"=3,"Often"=4),
 disc6  = c("Never"=1,"Less than once a year"=2,"A few times a year"=3,"A few times a month"=4,"At least once a week"=5,"Almost everyday"=6),
 hc5    = c("Never"=1,"Rarely"=2,"Sometimes"=3,"Most of the time"=4,"Always"=5),
 pss5   = c("Never"=1,"Almost Never"=2,"Sometimes"=3,"Fairly Often"=4,"Very Often"=5),
 food3  = c("Never true"=1,"Sometimes true"=2,"Often true"=3),
 sx4    = c("Not at all"=1,"Several days"=2,"More than half the days"=3,"Nearly every day"=4),
 happy6 = c("Extremely unhappy"=1,"Very unhappy"=2,"Moderately unhappy"=3,"Moderately happy"=4,"Very happy"=5,"Extremely happy"=6),
 mean5  = c("Not at all"=1,"A little"=2,"A moderate amount"=3,"Very much"=4,"An extreme amount"=5))
na_vals <- c("PMI: Skip","PMI: Dont Know","PMI: Prefer Not To Answer","NA","Skip","Prefer not to answer",
             "Don't know/Not sure","Does not apply to my neighborhood","Response removed due to invalid value")
spec <- function(id,scale,rev=FALSE) list(id=id,scale=scale,rev=rev)
fam <- list(
 lowcohesion = list(help=spec(40192463,"agree5",T), getalong=spec(40192411,"agree5",T),
                    trust=spec(40192499,"agree5",T), values=spec(40192417,"agree5",T)),
 disorder = list(safe=spec(40192384,"agree4",T), clean=spec(40192456,"agree4",T), upkeep=spec(40192386,"agree4",T),
                 watchout=spec(40192400,"agree4",T), graffiti=spec(40192420,"agree4"), noisy=spec(40192522,"agree4"),
                 vandalism=spec(40192412,"agree4"), abandoned=spec(40192469,"agree4"), hangaround=spec(40192500,"agree4"),
                 crime=spec(40192493,"agree4"), drug=spec(40192457,"agree4"), alcohol=spec(40192476,"agree4"),
                 trouble=spec(40192404,"agree4"), unsafe_night=spec(40192492,"sw4"), unsafe_day=spec(40192414,"sw4")),
 lowwalk = list(shops=spec(40192436,"sw4",T), transit=spec(40192440,"sw4",T), sidewalks=spec(40192437,"sw4",T),
                bike=spec(40192431,"sw4",T), recreation=spec(40192410,"sw4",T)),
 lowsupport = list(bed=spec(40192442,"time5",T), doctor=spec(40192480,"time5",T), meals=spec(40192388,"time5",T),
                   chores=spec(40192511,"time5",T), goodtime=spec(40192439,"time5",T), suggest=spec(40192528,"time5",T),
                   understand=spec(40192399,"time5",T), love=spec(40192446,"time5",T)),
 loneliness = list(lackcomp=spec(40192507,"lonely4"), noturnto=spec(40192397,"lonely4"), leftout=spec(40192398,"lonely4"),
                   isolated=spec(40192501,"lonely4"), withdrawn=spec(40192390,"lonely4"), aroundnot=spec(40192494,"lonely4"),
                   outgoing=spec(40192504,"lonely4",T), canfind=spec(40192516,"lonely4",T)),
 discrim = list(courtesy=spec(40192466,"disc6"), respect=spec(40192489,"disc6"), service=spec(40192416,"disc6"),
                notsmart=spec(40192490,"disc6"), afraid=spec(40192380,"disc6"), dishonest=spec(40192395,"disc6"),
                better=spec(40192496,"disc6"), names=spec(40192519,"disc6"), threatened=spec(40192451,"disc6")),
 discrim_hc = list(courtesy=spec(40192497,"hc5"), respect=spec(40192425,"hc5"), service=spec(40192503,"hc5"),
                   notsmart=spec(40192505,"hc5"), afraid=spec(40192423,"hc5"), better=spec(40192383,"hc5"),
                   notlisten=spec(40192394,"hc5")),
 stress = list(upset=spec(40192452,"pss5"), unable=spec(40192381,"pss5"), nervous=spec(40192491,"pss5"),
               cope=spec(40192506,"pss5"), angered=spec(40192396,"pss5"), piling=spec(40192462,"pss5"),
               confident=spec(40192419,"pss5",T), yourway=spec(40192525,"pss5",T), irritations=spec(40192449,"pss5",T),
               ontop=spec(40192445,"pss5",T)),
 foodinsec = list(food1=spec(40192517,"food3"), food2=spec(40192426,"food3")),
 lowwellbeing = list(happy=spec(1703980,"happy6",T), meaning=spec(1704001,"mean5",T)))
disab_ids <- c(903573,903574,903575,903576,903577,903578)
ace_ids   <- c(1332946,1332947,1332948,1332950,1333208,1332857,1333210,1333212,1332951,1333202,1333203)
trauma_ids<- c(1703988,1704021,1703976,1704055,1704035,1703994,1704056,1704011,1704029,1704020,1703989)
moves_id  <- 40192441
sx_ids <- c(1704026,1704024,1703983,1704039,1704004,1704041,1704038,1703996,1703977,
            1703984,1703995,1704000,1703987,1704028,1703920,1704042)

## ---------- pull + recode ----------
allids <- unique(c(unlist(lapply(fam, function(f) sapply(f, `[[`, "id"))),
                   disab_ids, ace_ids, trauma_ids, moves_id, sx_ids))
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(allids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$qid <- as.numeric(sv$question_concept_id)
W <- tidyr::pivot_wider(sv[,c("person_id","qid","answer")], names_from=qid, values_from=answer)
getcol <- function(id){ c <- W[[as.character(id)]]; if (is.null(c)) rep(NA,nrow(W)) else c }
rc <- function(x,scale){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(S[[scale]][x]) }
z  <- function(v) as.numeric(scale(v))
FAM <- data.frame(person_id=W$person_id)
for (fn in names(fam)){
  its <- fam[[fn]]; M <- sapply(names(its), function(nm){
    s<-its[[nm]]; zz<-z(rc(getcol(s$id), s$scale)); if (isTRUE(s$rev)) -zz else zz })
  FAM[[fn]] <- z(rowMeans(M, na.rm=TRUE)) }
cnt <- function(ids, pos="Yes|Once|More than once", neg="No|Never"){
  m <- sapply(ids, function(id){ a<-as.character(getcol(id)); ifelse(grepl(pos,a),1,ifelse(grepl(neg,a),0,NA)) })
  rowMeans(m, na.rm=TRUE) * ncol(m) }
FAM$disability <- z(cnt(disab_ids)); FAM$ace <- z(cnt(ace_ids)); FAM$trauma <- z(cnt(trauma_ids))
FAM$moves <- z(suppressWarnings(as.numeric(ifelse(getcol(moves_id) %in% na_vals, NA, getcol(moves_id)))))
SX <- sapply(sx_ids, function(id) z(rc(getcol(id), "sx4")))
FAM$overall <- z(rowMeans(SX, na.rm=TRUE))

## ---------- merge onto EHR cohort (same sample for both outcomes) ----------
cox$income_m <- ifelse(is.na(cox$income_n), median(cox$income_n,na.rm=TRUE), cox$income_n)
cox$sex_c <- factor(ifelse(cox$sex %in% c("Female","Male"), cox$sex, "Other"))
d <- merge(FAM, cox[,c("person_id","event","time_days","age","sex_c","income_m","log_util","n_cond")],
           by="person_id")
exposures <- c("discrim","discrim_hc","foodinsec","disability","ace","trauma",
               "lowcohesion","disorder","lowwalk","moves")
mediators <- c("loneliness","lowsupport","stress","lowwellbeing")
COV <- "age + sex_c + income_m"
Ex <- paste(exposures, collapse=" + "); Me <- paste(mediators, collapse=" + ")
cc <- complete.cases(d[, c(exposures, mediators, "overall", "event", "time_days", "age","sex_c","income_m")])
d <- d[cc, ]
cat(sprintf("Cross-modal sample (complete cases, both outcomes): n=%d | incident-dep events=%d\n",
            nrow(d), sum(d$event==1)))

## ---------- (I) self-report severity: lm total vs direct ----------
lm_t <- lm(as.formula(sprintf("overall ~ %s + %s", Ex, COV)), d)
lm_d <- lm(as.formula(sprintf("overall ~ %s + %s + %s", Ex, Me, COV)), d)
sr_t <- coef(lm_t)[exposures]; sr_d <- coef(lm_d)[exposures]

## ---------- (II) EHR incidence: Cox total vs direct (HR per SD) ----------
cx_t <- coxph(as.formula(sprintf("Surv(time_days, event) ~ %s + %s", Ex, COV)), d)
cx_d <- coxph(as.formula(sprintf("Surv(time_days, event) ~ %s + %s + %s", Ex, Me, COV)), d)
lhr_t <- coef(cx_t)[exposures]; lhr_d <- coef(cx_d)[exposures]

## ---------- cross-modal comparison ----------
pct <- function(a,b) ifelse(abs(a) < 0.02, NA, round(100*(a-b)/a))
CM <- data.frame(
  sr_total = round(sr_t,3),  sr_pct_distress = pct(sr_t, sr_d),
  ehr_HR   = round(exp(lhr_t),3), ehr_direct_HR = round(exp(lhr_d),3),
  ehr_pct_distress = pct(lhr_t, lhr_d))
cat("\n=== Cross-modal: same exposures, same distress block, same sample ===\n")
cat("  sr_total = self-report severity beta (per SD); ehr_HR = incident-depression hazard ratio (per SD)\n")
cat("  *_pct_distress = share absorbed by distress block (blank when effect ~0). Compare the two shares.\n")
print(CM)

## ---------- ascertainment fork: does discrimination predict LOWER utilization? ----------
u1 <- lm(as.formula(sprintf("log_util ~ %s + %s", Ex, COV)), d)
u2 <- lm(as.formula(sprintf("n_cond  ~ %s + %s", Ex, COV)), d)
U <- data.frame(log_util = round(coef(u1)[exposures],3), n_cond = round(coef(u2)[exposures],3))
cat("\n=== Ascertainment: exposure -> healthcare utilization (joint, adj covs) ===\n")
cat("  NEGATIVE log_util/n_cond for discrim or discrim_hc => under-utilization => EHR under-diagnosis risk\n")
print(U)
cat(sprintf("\n  n_cond and event share the EHR channel: cor(n_cond, event) = %.3f\n",
            cor(d$n_cond, d$event, use="complete.obs")))
cat("\nRead: if a discrimination effect on self-report severity does NOT appear in EHR incidence AND\n",
    "discrimination predicts lower utilization -> the gap is differential ascertainment, not 'distress only'.\n", sep="")
