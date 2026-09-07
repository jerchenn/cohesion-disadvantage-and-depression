# exposome_map.R -- EXPLORATORY. Broad family-level map: how each predictor FAMILY associates with
# each depression/anxiety SYMPTOM dimension. Every predictor family is oriented to its ADVERSE pole
# (higher = more risk); symptoms are oriented higher = more severe. So a positive beta = risk factor.
# Cross-sectional, mostly self-report (common-method) -> read as a STRUCTURE map, not causal.
# Adjusted (age, sex, income). Aggregate only.
#
# INPUT: cox_dat.rds (covariates) + ds_survey.  Needs: bigrquery.

library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

## ---------- answer scales ----------
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

## ---------- item spec: name -> (id, scale, reverse?) ; reverse=TRUE flips z so higher = adverse ----------
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

## count families (binarize -> count of adversities)
disab_ids <- c(903573,903574,903575,903576,903577,903578)
ace_ids   <- c(1332946,1332947,1332948,1332950,1333208,1332857,1333210,1333212,1332951,1333202,1333203)
trauma_ids<- c(1703988,1704021,1703976,1704055,1704035,1703994,1704056,1704011,1704029,1704020,1703989)
moves_id  <- 40192441

## symptom outcomes (higher = more severe)
phq <- c(anhedonia=1704026,down=1704024,sleep=1703983,fatigue=1704039,appetite=1704004,
         worthless=1704041,concen=1704038,psychomotor=1703996,suicidal=1703977)
gad <- c(nervous=1703984,worrystop=1703995,worrymuch=1704000,relax=1703987,restless=1704028,
         irritable=1703920,afraid=1704042)
symp <- c(phq,gad)

## ---------- pull everything ----------
allids <- unique(c(unlist(lapply(fam, function(f) sapply(f, `[[`, "id"))),
                   disab_ids, ace_ids, trauma_ids, moves_id, symp))
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(allids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$qid <- as.numeric(sv$question_concept_id)
W <- tidyr::pivot_wider(sv[,c("person_id","qid","answer")], names_from=qid, values_from=answer)
getcol <- function(id) { c <- W[[as.character(id)]]; if (is.null(c)) rep(NA,nrow(W)) else c }
rc <- function(x,scale){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(S[[scale]][x]) }
z  <- function(v) as.numeric(scale(v))

## scale-family indices (oriented to adverse pole)
FAM <- data.frame(person_id=W$person_id)
for (fn in names(fam)){
  its <- fam[[fn]]; M <- sapply(names(its), function(nm){
    s<-its[[nm]]; zz<-z(rc(getcol(s$id), s$scale)); if (isTRUE(s$rev)) -zz else zz })
  FAM[[fn]] <- z(rowMeans(M, na.rm=TRUE))
}
## count families
cnt <- function(ids, pos_regex="Yes|Once|More than once", neg_regex="^No$|^Never$"){
  m <- sapply(ids, function(id){ a<-as.character(getcol(id))
    ifelse(grepl(pos_regex,a), 1, ifelse(grepl(neg_regex,a), 0, NA)) })
  rowMeans(m, na.rm=TRUE) * ncol(m) }               # mean*k ~ count, NA-robust
FAM$disability <- z(cnt(disab_ids))
FAM$ace        <- z(cnt(ace_ids))
FAM$trauma     <- z(cnt(trauma_ids))
FAM$moves      <- z(suppressWarnings(as.numeric(ifelse(getcol(moves_id) %in% na_vals, NA, getcol(moves_id)))))

## symptom outcomes
SY <- data.frame(person_id=W$person_id)
for (nm in names(symp)) SY[[nm]] <- z(rc(getcol(symp[nm]), "sx4"))

## ---------- assemble + covariates ----------
d <- merge(FAM, SY, by="person_id")
cv <- cox[, c("person_id","age","sex","income_n")]
cv$income_m <- ifelse(is.na(cv$income_n), median(cv$income_n,na.rm=TRUE), cv$income_n)
cv$sex_c <- factor(ifelse(cv$sex %in% c("Female","Male"), cv$sex, "Other"))
d <- merge(d, cv[,c("person_id","age","sex_c","income_m")], by="person_id")
families <- c(names(fam),"disability","ace","trauma","moves")
cat(sprintf("Analytic n (any symptom + covariates): %d\n", nrow(d)))
cat("\n=== Family coverage (non-missing index) ===\n")
print(sapply(families, function(f) sum(is.finite(d[[f]]))))

## ---------- family x symptom association matrix (adjusted beta; + = risk) ----------
beta <- function(y,x){ m<-lm(as.formula(sprintf("%s ~ %s + age + sex_c + income_m", y, x)), d)
                       coef(m)[[x]] }
M <- outer(families, names(symp), Vectorize(function(f,s) round(beta(s,f),2)))
dimnames(M) <- list(families, names(symp))
cat("\n=== Family (adverse pole) x symptom dimension (adj beta; + = more risk -> more symptom) ===\n")
print(M); write.csv(M, "exposome_symptom_matrix.csv")

## ---------- structure: family intercorrelations ----------
cat("\n=== Family intercorrelations (Pearson r) ===\n")
print(round(cor(d[,families], use="pairwise.complete.obs"), 2))
cat("\nOrientation: predictors higher = more adverse; symptoms higher = more severe; + beta = risk.\n")
