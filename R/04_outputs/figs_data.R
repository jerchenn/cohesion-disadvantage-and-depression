# figs_data.R -- recompute EVERY figure estimate from saved .rds (no hardcoding;
# several estimates exist in 2 vintages this session, so recompute is the safe path).
# Writes one tidy CSV consumed by figs_plot.R. Also prints the primary under the
# collapsed covariate spec to reconcile with the factor-spec headline (0.882).
library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
COVS <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
row <- function(panel,label,m,v,extra=NA){
  ci<-confint(m); data.frame(panel,label,
    hr=exp(coef(m)[[v]]), lo=exp(ci[v,1]), hi=exp(ci[v,2]),
    p=summary(m)$coefficients[v,"Pr(>|z|)"], n=m$n, events=m$nevent, extra, stringsAsFactors=FALSE) }
cx <- function(rhs,d,cv=COVS) coxph(as.formula(paste("Surv(time_days,event) ~",rhs,"+",cv)), d)
out <- list()

## ---- primary (collapsed spec) full + lag period : reconcile with 0.882 ----
cox <- collapse(readRDS("cox_dat3.rds"))
mP <- cx("z_cohesion", cox)
cat(sprintf("RECONCILE primary (collapsed spec): HR %.4f (%.4f-%.4f)  [factor-spec headline was 0.882]\n",
    exp(coef(mP)[["z_cohesion"]]), exp(confint(mP)["z_cohesion",1]), exp(confint(mP)["z_cohesion",2])))
out[["p1"]] <- row("primary","Cohesion (primary)", mP, "z_cohesion")
out[["p2"]] <- row("primary","Cohesion (180-d lag period)", cx("z_cohesion",cox[cox$time_days>180,]), "z_cohesion")

## ---- Effect modification by income (reported in Results; drawn as gradient in eFigure 2) ----
cox$low_income <- as.integer(cox$income_m<=4)
cvNoInc <- "age + sex_c + race_c + ethn_c + educ_m + log_util"
mi <- cx("z_cohesion*low_income", cox, cvNoInc); pI <- summary(mi)$coefficients["z_cohesion:low_income","Pr(>|z|)"]
out[["i1"]] <- row("mod_income","Higher income (>$50k)", cx("z_cohesion",cox[cox$low_income==0,],cvNoInc),"z_cohesion", pI)
out[["i2"]] <- row("mod_income","Low income ($50k or less)", cx("z_cohesion",cox[cox$low_income==1,],cvNoInc),"z_cohesion", pI)

## ---- Effect modification by AREA deprivation, ZIP3 (Results; gradient in eFigure 2) ----
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
cx3 <- readRDS("cox_dat3.rds"); cx3$deprivation_index <- NULL
g <- collapse(merge(cx3, geo, by="person_id"))
g <- g[!is.na(g$z_cohesion) & !is.na(g$deprivation_index),]; g$z_dep <- scale(g$deprivation_index)[,1]
md <- cx("z_cohesion*z_dep", g, paste(COVS,"+ z_dep")); pD <- summary(md)$coefficients["z_cohesion:z_dep","Pr(>|z|)"]
g$dep_t <- cut(g$deprivation_index, quantile(g$deprivation_index,0:3/3,na.rm=TRUE),
               include.lowest=TRUE, labels=c("Low","Middle","High"))
for (t in c("Low","Middle","High"))
  out[[paste0("d_",t)]] <- row("mod_deprivation", paste0(t," deprivation"),
      cx("z_cohesion", g[g$dep_t==t,], paste(COVS,"+ z_dep")), "z_cohesion", pD)

## ---- Figure 3: multimodal concordance (each modality -> incident depression, per SD healthier) ----
wb <- collapse(readRDS("wb_item_dat.rds"))
fb <- collapse(readRDS("fitbit_dat.rds"))
fb$z_steps<-scale(fb$steps)[,1]; fb$z_mvpa<-scale(fb$mvpa_min)[,1]
fb$z_eff<-scale(fb$efficiency)[,1]; fb$z_waso<- -scale(fb$waso_min)[,1]
out[["c0"]] <- row("concordance","Cohesion (exposure)", mP, "z_cohesion")
out[["c1"]] <- row("concordance","Self-report: happiness", cx("z_happy",wb),"z_happy")
out[["c2"]] <- row("concordance","Self-report: meaning",   cx("z_meaning",wb),"z_meaning")
out[["c3"]] <- row("concordance","Self-report: social connection", cx("z_cutoff",wb),"z_cutoff")
out[["c4"]] <- row("concordance","Device: steps",  cx("z_steps",fb),"z_steps")
out[["c5"]] <- row("concordance","Device: MVPA",   cx("z_mvpa",fb),"z_mvpa")
out[["c6"]] <- row("concordance","Device: sleep efficiency", cx("z_eff",fb),"z_eff")
out[["c7"]] <- row("concordance","Device: sleep WASO",       cx("z_waso",fb),"z_waso")

est <- do.call(rbind, out)
names(est)[names(est)=="extra"] <- "interaction_p"
write.csv(est, "fig_estimates.csv", row.names=FALSE)
cat("\nWrote fig_estimates.csv\n"); print(est, row.names=FALSE, digits=3)
