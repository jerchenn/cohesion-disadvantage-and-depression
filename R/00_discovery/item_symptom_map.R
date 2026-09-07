# item_symptom_map.R -- EXPLORATORY. Item-level association map from three predictor families
#   (built-environment items, subjective-wellbeing items, device physical activity)
# to individual depression (PHQ) and anxiety (GAD) SYMPTOM items.
# Everything is oriented so higher = healthier, so a positive beta = a protective association.
# Cross-sectional; survey predictors and symptoms share method (common-method), so read this as a
# STRUCTURE map (which feature tracks which symptom cluster), not causal effect. Adjusted (age, sex,
# income), HC3 SE. Aggregate output only.
#
# INPUT: wb_item_dat.rds (+ fitbit_dat.rds for the activity block).  Needs: sandwich, lmtest.

library(sandwich); library(lmtest)
wb <- readRDS("wb_item_dat.rds")
collapse <- function(d){
  d$sex_c    <- factor(ifelse(d$sex %in% c("Female","Male"), d$sex, "Other"))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n, na.rm=TRUE), d$income_n); d }
wb <- collapse(wb)

## ---- predictor families (raw recoded env/cohesion items are oriented higher = better; z-score them) ----
benv <- c("safe","crime","unsafe_night","unsafe_day","drug","alcohol","trouble","noisy","hangaround",
          "graffiti","vandalism","abandoned","clean","upkeep",                 # safety / disorder / decay
          "shops","transit","sidewalks","bike","recreation")                   # amenity / walkability
cohesion <- c("help","getalong","trust","values","watchout")
swb <- c("happy","meaning","cutoff")                                           # already z_ in wb_item_dat
for (nm in c(benv, cohesion)) wb[[paste0("z_",nm)]] <- scale(wb[[nm]])[,1]

## ---- symptom outcomes (z_ already present; higher = healthier) ----
phq <- c("down","anhedonia","fatigue","concen","worthless","appetite","sleep","restless","psychomotor","suicidal")
gad <- c("nervous","worrymuch","worrystop","relax","irritable","afraid")
symptoms <- c(phq, gad)

beta <- function(y, x, data){
  m <- lm(as.formula(sprintf("z_%s ~ z_%s + age + sex_c + income_m", y, x)), data)
  tryCatch(coeftest(m, vcov=vcovHC(m,"HC3"))[paste0("z_",x),"Estimate"], error=function(e) NA) }
mk <- function(preds, data){
  M <- outer(preds, symptoms, Vectorize(function(p,s) beta(s,p,data)))
  dimnames(M) <- list(preds, symptoms); round(M, 2) }

## ---- (1) built-environment + cohesion items -> symptoms ----
cat(sprintf("=== Built-environment + cohesion items x symptoms (adj beta; + = protective; n=%d) ===\n",
            sum(!is.na(wb$z_safe))))
Mb <- mk(c(benv, cohesion), wb); print(Mb); write.csv(Mb, "map_builtenv_symptom.csv")

## ---- (2) subjective-wellbeing items -> symptoms ----
cat("\n=== Subjective-wellbeing items x symptoms ===\n")
Ms <- mk(swb, wb); print(Ms); write.csv(Ms, "map_swb_symptom.csv")

## ---- (3) device physical activity -> symptoms (Fitbit n EHHWB intersection) ----
fb <- readRDS("fitbit_dat.rds")
act <- collapse(merge(wb, fb[, c("person_id","steps","mvpa_min","sed_min","efficiency","waso_min")], by="person_id"))
act$z_steps <- scale(act$steps)[,1]; act$z_mvpa <- scale(act$mvpa_min)[,1]
act$z_sed <- -scale(act$sed_min)[,1]; act$z_eff <- scale(act$efficiency)[,1]; act$z_waso <- -scale(act$waso_min)[,1]
actv <- c("steps","mvpa","sed","eff","waso")                                   # oriented higher = healthier
cat(sprintf("\n=== Device activity x symptoms (Fitbit n EHHWB, n=%d; sed/waso reversed) ===\n",
            sum(!is.na(act$steps))))
Ma <- mk(actv, act); print(Ma); write.csv(Ma, "map_activity_symptom.csv")

## ---- (4) quick structure read: top 3 predictors per symptom, across all families ----
All <- rbind(Mb, Ms, Ma)
cat("\n=== Top 3 predictors per symptom (by |beta|) ===\n")
for (s in symptoms){
  o <- order(-abs(All[,s]))[1:3]
  cat(sprintf("  %-11s : %s\n", s,
      paste(sprintf("%s %+.2f", rownames(All)[o], All[o,s]), collapse="   "))) }
cat("\nOrientation: all items higher = healthier; a positive beta means the healthier/better level of\n",
    "the predictor goes with fewer symptoms. Survey predictors and symptoms share method.\n", sep="")
