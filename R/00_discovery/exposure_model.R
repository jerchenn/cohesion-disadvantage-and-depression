# exposure_model.R -- EXPLORATORY. Keeps two kinds of variable strictly apart:
#   EXPOSURES  = things measured independently of mood (discrimination, food insecurity, disability,
#                ACEs, trauma, low cohesion, disorder, low walkability, low support, residential moves)
#   MEDIATORS  = self-reported distress that sits right next to the outcome (loneliness, stress, low
#                wellbeing) -- candidate mediators / near-outcomes, NOT put on the same footing as exposures.
# Outcomes: 16 PHQ/GAD symptoms collapsed into 3 clinical composites (somatic, mood, anxiety) + overall.
# Blocks:
#   A  exposure -> composite, UNIVARIATE (each exposure alone, adj covariates)
#   B  exposure -> composite, JOINT     (all exposures mutually adjusted) -> which survive
#   C  exposure -> each mediator, JOINT (a-paths: which exposures feed loneliness/stress/low wellbeing)
#   D  descriptive mediation on overall severity: total (exposures+covs) vs direct (+mediators), % via block
# Cross-sectional, self-report -> STRUCTURE map, not causal. Mediated share is an UPPER bound (mediators
# overlap the outcome). Adjusted age/sex/income. Aggregate only.
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

## ---------- item spec ----------
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

## ---------- symptom composites (higher = more severe) ----------
mood    <- c(anhedonia=1704026,down=1704024,worthless=1704041,suicidal=1703977,concen=1704038)
somatic <- c(sleep=1703983,fatigue=1704039,appetite=1704004,psychomotor=1703996)
anxiety <- c(nervous=1703984,worrystop=1703995,worrymuch=1704000,relax=1703987,
             restless=1704028,irritable=1703920,afraid=1704042)
symp <- c(mood,somatic,anxiety)

## ---------- pull ----------
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

## family indices (adverse pole)
FAM <- data.frame(person_id=W$person_id)
for (fn in names(fam)){
  its <- fam[[fn]]; M <- sapply(names(its), function(nm){
    s<-its[[nm]]; zz<-z(rc(getcol(s$id), s$scale)); if (isTRUE(s$rev)) -zz else zz })
  FAM[[fn]] <- z(rowMeans(M, na.rm=TRUE))
}
cnt <- function(ids, pos_regex="Yes|Once|More than once", neg_regex="No|Never"){
  m <- sapply(ids, function(id){ a<-as.character(getcol(id))
    ifelse(grepl(pos_regex,a), 1, ifelse(grepl(neg_regex,a), 0, NA)) })
  rowMeans(m, na.rm=TRUE) * ncol(m) }
FAM$disability <- z(cnt(disab_ids))
FAM$ace        <- z(cnt(ace_ids))
FAM$trauma     <- z(cnt(trauma_ids))
FAM$moves      <- z(suppressWarnings(as.numeric(ifelse(getcol(moves_id) %in% na_vals, NA, getcol(moves_id)))))

## symptom composites
zitems <- function(ids) sapply(names(ids), function(nm) z(rc(getcol(ids[nm]), "sx4")))
mk_comp <- function(ids){ Z <- zitems(ids); list(score=z(rowMeans(Z, na.rm=TRUE)), Z=Z) }
Cm <- mk_comp(mood); Cs <- mk_comp(somatic); Ca <- mk_comp(anxiety)
allZ <- cbind(Cm$Z, Cs$Z, Ca$Z)
COMP <- data.frame(person_id=W$person_id,
                   mood=Cm$score, somatic=Cs$score, anxiety=Ca$score,
                   overall=z(rowMeans(allZ, na.rm=TRUE)))

## ---------- assemble + covariates ----------
d <- merge(FAM, COMP, by="person_id")
cv <- cox[, c("person_id","age","sex","income_n")]
cv$income_m <- ifelse(is.na(cv$income_n), median(cv$income_n,na.rm=TRUE), cv$income_n)
cv$sex_c <- factor(ifelse(cv$sex %in% c("Female","Male"), cv$sex, "Other"))
d <- merge(d, cv[,c("person_id","age","sex_c","income_m")], by="person_id")

exposures <- c("discrim","discrim_hc","foodinsec","disability","ace","trauma",
               "lowcohesion","disorder","lowwalk","lowsupport","moves")
mediators <- c("loneliness","stress","lowwellbeing")
comps     <- c("mood","somatic","anxiety")
COV <- "age + sex_c + income_m"

cat(sprintf("Analytic n: %d\n", nrow(d)))
mir <- function(Z){ r <- cor(Z, use="pairwise.complete.obs"); mean(r[upper.tri(r)]) }
cat(sprintf("Composite internal consistency (mean inter-item r): mood %.2f | somatic %.2f | anxiety %.2f\n",
            mir(Cm$Z), mir(Cs$Z), mir(Ca$Z)))
cat(sprintf("Composite intercorrelations: mood-somatic %.2f | mood-anxiety %.2f | somatic-anxiety %.2f\n",
            cor(d$mood,d$somatic,use="complete.obs"), cor(d$mood,d$anxiety,use="complete.obs"),
            cor(d$somatic,d$anxiety,use="complete.obs")))

## ---------- A. univariate exposure -> composite ----------
uni <- function(y,x){ coef(lm(as.formula(sprintf("%s ~ %s + %s", y, x, COV)), d))[[x]] }
A <- outer(exposures, comps, Vectorize(function(x,y) round(uni(y,x),3)))
dimnames(A) <- list(exposures, comps)
cat("\n=== A. exposure -> composite, UNIVARIATE (adj covariates) ===\n"); print(A)

## ---------- B. joint exposure -> composite ----------
joint_block <- function(y, xs){
  m <- lm(as.formula(sprintf("%s ~ %s + %s", y, paste(xs,collapse=" + "), COV)), d)
  round(coef(m)[xs],3) }
B <- sapply(comps, function(y) joint_block(y, exposures))
cat("\n=== B. exposure -> composite, JOINT (all exposures mutually adjusted) ===\n"); print(B)

## ---------- C. a-paths: exposure -> each mediator, joint ----------
Cmat <- sapply(mediators, function(m) joint_block(m, exposures))
cat("\n=== C. exposure -> mediator, JOINT (a-paths; which exposures feed the distress block) ===\n"); print(Cmat)

## ---------- D. descriptive mediation on OVERALL severity ----------
tot <- lm(as.formula(sprintf("overall ~ %s + %s", paste(exposures,collapse=" + "), COV)), d)
dir <- lm(as.formula(sprintf("overall ~ %s + %s + %s",
          paste(exposures,collapse=" + "), paste(mediators,collapse=" + "), COV)), d)
bt <- coef(tot)[exposures]; bd <- coef(dir)[exposures]
D <- data.frame(total=round(bt,3), direct=round(bd,3),
                pct_via_distress=round(100*(bt-bd)/bt))
cat("\n=== D. overall severity: total (exposures+covs) vs direct (+distress block) ===\n")
cat("    total = joint effect holding other exposures fixed; direct = after adding loneliness+stress+lowwellbeing\n")
cat("    pct_via_distress = share absorbed by the distress block (UPPER bound; mediators overlap outcome)\n")
print(D)
cat(sprintf("\n    distress-block b-paths (in direct model): loneliness %.3f | stress %.3f | lowwellbeing %.3f\n",
            coef(dir)[["loneliness"]], coef(dir)[["stress"]], coef(dir)[["lowwellbeing"]]))
cat("\nOrientation: all predictors higher = more adverse; composites higher = more severe; + = risk.\n")
