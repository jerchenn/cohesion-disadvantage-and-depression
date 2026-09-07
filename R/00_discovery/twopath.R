# twopath.R -- Paper 2 core analysis: TWO neighborhood pathways to incident depression.
#   (A) SOCIAL pathway  : safety / low social disorder -> trust & willingness to help
#                         (cognitive social capital, self-reported) -> lower depression.
#   (B) AMENITY pathway : walkability / amenities -> objectively-measured physical activity
#                         (wearable steps) -> lower depression, INDEPENDENT of trust/help.
# The device-measured steps node anchors pathway B against common-method bias (both indices and
# trust/help are self-reported on the same SDOH survey).
#
# Three built-environment indices (social-safety, physical-decay, amenity) -- NOT pre-collapsed.
# Primary evidence = a-paths and b-paths reported separately; proportion-mediated is secondary
# (cross-sectional mediators; non-collapsible Cox -> descriptive decomposition, not identified NIE/NDE).
# Aggregate output only (All of Us DUA).
#
# INPUT: cox_dat.rds + fitbit_dat.rds + ds_survey + ds_zip_code_socioeconomic.  RUN AFTER build_cox.R.
# Needs: tidyverse, survival, bigrquery.

library(tidyverse); library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds"); fb <- readRDS("fitbit_dat.rds")

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
COV <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"

## ---------- pull built-environment items + trust/help, scale-aware recode ----------
ids <- c(safe=40192384, clean=40192456, upkeep=40192386, sidewalks=40192437, shops=40192436,
         transit=40192440, recreation=40192410, bike=40192431, crime=40192493,
         unsafe_night=40192492, unsafe_day=40192414, graffiti=40192420, vandalism=40192412,
         abandoned=40192469, drug=40192457, alcohol=40192476, trouble=40192404, noisy=40192522,
         hangaround=40192500, trust=40192499, help=40192463)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(ids)[match(as.numeric(sv$question_concept_id), ids)]
wide <- tidyr::pivot_wider(sv[,c("person_id","item","answer")], names_from=item, values_from=answer)

map_sw <- c("Strongly disagree"=1,"Somewhat disagree"=2,"Somewhat agree"=3,"Strongly agree"=4)
map_ad <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
map_5  <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
sw_items <- c("sidewalks","shops","transit","recreation","bike","unsafe_night","unsafe_day")
d5_items <- c("trust","help")
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Dont Know","PMI: Prefer Not To Answer",
             "Skip","Don't know","Prefer not to answer","Does not apply to my neighborhood")
scale_of <- function(nm) if (nm %in% sw_items) map_sw else if (nm %in% d5_items) map_5 else map_ad
rc <- function(x, m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in names(ids)) wide[[nm]] <- rc(wide[[nm]], scale_of(nm))

## z-score; flip disorder/danger ("bad") items so higher z = better environment on every item
bad <- c("crime","unsafe_night","unsafe_day","graffiti","vandalism","abandoned",
         "drug","alcohol","trouble","noisy","hangaround")
for (nm in names(ids)) wide[[paste0("z_",nm)]] <- as.numeric(scale(wide[[nm]]))
for (nm in bad)       wide[[paste0("z_",nm)]] <- -wide[[paste0("z_",nm)]]

## ---------- THREE indices (mean of member z-scores; higher = better environment) ----------
social_items <- c("safe","crime","unsafe_night","unsafe_day","drug","alcohol","trouble","noisy","hangaround")
decay_items  <- c("clean","upkeep","abandoned","graffiti","vandalism")
amenity_items<- c("recreation","shops","transit","sidewalks","bike")
zmean <- function(df, items){ m <- as.matrix(df[, paste0("z_", items)]); rowMeans(m, na.rm=TRUE) }
wide$social  <- zmean(wide, social_items)
wide$decay   <- zmean(wide, decay_items)
wide$amenity <- zmean(wide, amenity_items)
wide$cohesion_th <- rowMeans(as.matrix(wide[, c("z_trust","z_help")]), na.rm=TRUE)  # trust+help mediator
for (v in c("social","decay","amenity","cohesion_th")) wide[[v]] <- as.numeric(scale(wide[[v]]))

## ---------- coverage: per-index N and pairwise-complete N (guards the amenity contrast) ----------
cov_n <- sapply(c("social","decay","amenity","cohesion_th"),
                function(v) sum(!is.na(wide[[v]]) & is.finite(wide[[v]])))
cat("=== Per-index N (non-missing) ===\n"); print(cov_n)
cc <- complete.cases(wide[, c("social","decay","amenity","cohesion_th")])
cat("Pairwise-complete N (all indices + mediator):", sum(cc), "\n")
cat("\n=== Index separability (Pearson r) ===\n")
print(round(cor(wide[cc, c("social","decay","amenity","cohesion_th")]), 2))

## ---------- assemble analysis frames ----------
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
W <- wide[, c("person_id","social","decay","amenity","cohesion_th")]
d  <- collapse(merge(cox, W, by="person_id"))                    # full cohort (self-report mediator)
d  <- merge(d, geo, by="person_id", all.x=TRUE); d$z_dep <- as.numeric(scale(d$deprivation_index))
df <- collapse(merge(fb,  W, by="person_id"))                    # fitbit subcohort (objective mediator)
df <- merge(df, geo, by="person_id", all.x=TRUE)
df <- df[!is.na(df$steps), ]                                     # restrict to Fitbit wearers -> one consistent sample
df$z_dep   <- as.numeric(scale(df$deprivation_index))
df$z_steps <- as.numeric(scale(df$steps))
df$z_mvpa  <- as.numeric(scale(df$mvpa_min))
cat(sprintf("\nAnalytic N: full cohort %d (events %d) | fitbit subcohort %d (events %d)\n",
            nrow(d), sum(d$event), nrow(df), sum(df$event)))

std_beta <- function(y, x, data, cov="age + income_m + z_dep"){
  m <- lm(as.formula(paste(y,"~",x,"+",cov)), data); summary(m)$coefficients[x,"Estimate"] }
hrci <- function(m,v) sprintf("%.3f (%.3f-%.3f)", exp(coef(m)[[v]]), exp(confint(m)[v,1]), exp(confint(m)[v,2]))

## ---------- A-PATHS (exposure -> mediator), reported separately (PRIMARY evidence) ----------
cat("\n=== A-paths: built-env index -> MEDIATORS (std beta) ===\n")
cat(sprintf("  %-8s  ->trust/help    ->steps     ->MVPA\n",""))
for (v in c("social","decay","amenity")){
  b_th <- std_beta("cohesion_th", v, d)
  b_st <- std_beta("z_steps",     v, df)
  b_mv <- std_beta("z_mvpa",      v, df)
  cat(sprintf("  %-8s  %+.3f          %+.3f      %+.3f\n", v, b_th, b_st, b_mv)) }

## ---------- B-PATHS (mediator -> depression) ----------
cat("\n=== B-paths: MEDIATOR -> incident depression (HR per SD healthier) ===\n")
mB1 <- coxph(as.formula(paste("Surv(time_days,event) ~ cohesion_th +", COV)), d)
cat("  trust/help (full cohort) HR", hrci(mB1,"cohesion_th"), "\n")
mB2 <- coxph(as.formula(paste("Surv(time_days,event) ~ z_steps +", COV)), df)
cat("  steps (wearer subcohort) HR", hrci(mB2,"z_steps"), "\n")
mB3 <- coxph(as.formula(paste("Surv(time_days,event) ~ z_mvpa +", COV)), df)
cat("  MVPA  (wearer subcohort) HR", hrci(mB3,"z_mvpa"), "\n")

## ---------- Index -> depression: total, joint, and pathway-specific direct effects ----------
cat("\n=== Index -> incident depression (HR per SD better environment) ===\n")
for (v in c("social","decay","amenity")){
  m <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "+", COV)), d)
  cat(sprintf("  %-8s total   HR %s\n", v, hrci(m,v))) }
mJ <- coxph(as.formula(paste("Surv(time_days,event) ~ social + decay + amenity +", COV)), d)
cat("  --- joint (mutually adjusted) ---\n")
for (v in c("social","decay","amenity")) cat(sprintf("  %-8s joint   HR %s\n", v, hrci(mJ,v)))

## ---------- SOCIAL pathway: mediation through trust/help (full cohort) ----------
cat("\n=== SOCIAL pathway: social-safety -> depression, direct after trust/help ===\n")
mt <- coxph(as.formula(paste("Surv(time_days,event) ~ social +", COV)), d)
md <- coxph(as.formula(paste("Surv(time_days,event) ~ social + cohesion_th +", COV)), d)
pm <- 1 - log(exp(coef(md)[["social"]])) / log(exp(coef(mt)[["social"]]))
cat(sprintf("  total  %s | direct %s | proportion via trust/help ~ %.0f%% (secondary)\n",
            hrci(mt,"social"), hrci(md,"social"), 100*pm))

## ---------- AMENITY pathway: mediation through device steps (fitbit subcohort) ----------
cat("\n=== AMENITY pathway: amenity -> depression, direct after device activity (wearer subcohort, one sample) ===\n")
at  <- coxph(as.formula(paste("Surv(time_days,event) ~ amenity +", COV)), df)
ads <- coxph(as.formula(paste("Surv(time_days,event) ~ amenity + z_steps +", COV)), df)
adm <- coxph(as.formula(paste("Surv(time_days,event) ~ amenity + z_mvpa +", COV)), df)
ps  <- 1 - log(exp(coef(ads)[["amenity"]])) / log(exp(coef(at)[["amenity"]]))
pmv <- 1 - log(exp(coef(adm)[["amenity"]])) / log(exp(coef(at)[["amenity"]]))
cat(sprintf("  total        %s\n  direct|steps %s  (via steps ~%.0f%%)\n  direct|MVPA  %s  (via MVPA ~%.0f%%)\n",
            hrci(at,"amenity"), hrci(ads,"amenity"), 100*ps, hrci(adm,"amenity"), 100*pmv))
cat("  cross-check: amenity -> trust/help (should be ~0):", sprintf("%+.3f", std_beta("cohesion_th","amenity", d)), "\n")

## ---------- add subjective wellbeing (SWB) as a third mediator (EHHWB subcohort) ----------
## SWB = mean of z_happy, z_meaning, z_cutoff (already recoded in wb_item_dat; higher = better).
## NB: SWB overlaps subclinical depression, so its b-path is inflated -- treat as an affective
## correlate, not a clean mechanism (same caveat as Paper 1's affect mediation).
wb  <- readRDS("wb_item_dat.rds")
wb$swb <- rowMeans(as.matrix(wb[, c("z_happy","z_meaning","z_cutoff")]), na.rm=TRUE)
dS  <- merge(d, wb[, c("person_id","swb")], by="person_id")
dS  <- dS[is.finite(dS$swb), ]; dS$swb <- as.numeric(scale(dS$swb))
cat(sprintf("\n=== Subjective wellbeing (SWB) mediator: EHHWB subcohort n %d (events %d) ===\n",
            nrow(dS), sum(dS$event)))
cat("A-paths (index -> SWB, std beta):\n")
for (v in c("social","decay","amenity")) cat(sprintf("  %-8s -> SWB  %+.3f\n", v, std_beta("swb", v, dS)))
cat("  trust/help -> SWB  ", sprintf("%+.3f", std_beta("swb","cohesion_th", dS)), "\n")
mS <- coxph(as.formula(paste("Surv(time_days,event) ~ swb +", COV)), dS)
cat("B-path: SWB -> incident depression HR", hrci(mS,"swb"), " (inflated by construct overlap)\n")
## social pathway: nested attenuation through trust/help then + SWB (mutually adjusted)
m0 <- coxph(as.formula(paste("Surv(time_days,event) ~ social +", COV)), dS)
m1 <- coxph(as.formula(paste("Surv(time_days,event) ~ social + cohesion_th +", COV)), dS)
m2 <- coxph(as.formula(paste("Surv(time_days,event) ~ social + cohesion_th + swb +", COV)), dS)
cat(sprintf("Social -> depression:  total %s | +trust/help %s | +trust/help+SWB %s\n",
            hrci(m0,"social"), hrci(m1,"social"), hrci(m2,"social")))

## ---------- robustness: 180-day lag on the two total effects ----------
cat("\n=== 180-day lag period (total effects) ===\n")
mtl <- coxph(as.formula(paste("Surv(time_days,event) ~ social +", COV)), d[d$time_days>180,])
atl <- coxph(as.formula(paste("Surv(time_days,event) ~ amenity +", COV)), df[df$time_days>180,])
cat("  social  HR", hrci(mtl,"social"), " | amenity HR", hrci(atl,"amenity"), "\n")

cat("\nCaveats: mediators are cross-sectional (same SDOH wave for trust/help); Cox HRs are\n",
    "non-collapsible; proportions are descriptive decompositions, not identified NIE/NDE.\n", sep="")
