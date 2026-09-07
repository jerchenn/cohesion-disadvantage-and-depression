# etable_robustness.R -- Paper 1 supplement eTable: cohesion -> incident depression, robustness to the
# individual/structural adversity exposome, pooled and by area-deprivation tertile.
#   Panel A: pooled re-anchored cohort -- base spec, each adversity added singly, then all five jointly.
#   Panel B: by area-deprivation tertile -- base vs + full exposome.
#   Footnote: published full-cohort primary estimate + cohesion x deprivation interaction.
# Covariate parameterization matches fig_stratified.R (collapsed race/ethnicity, continuous income/
# education) so every model converges cleanly (no sparse addNA factor cells -> no separation). Re-anchored
# cohort (cox_dat2.rds) so adversity measures precede follow-up. Aggregate only; writes etable_robustness.csv.
#
# INPUT: cox_dat.rds + cox_dat2.rds + ds_survey + ds_zip_code_socioeconomic.  Needs: bigrquery, survival.

library(bigrquery); library(survival)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

## ---------- adversity families ----------
S <- list(disc6=c("Never"=1,"Less than once a year"=2,"A few times a year"=3,"A few times a month"=4,"At least once a week"=5,"Almost everyday"=6),
          food3=c("Never true"=1,"Sometimes true"=2,"Often true"=3))
na_vals <- c("PMI: Skip","PMI: Dont Know","PMI: Prefer Not To Answer","NA","Skip","Prefer not to answer",
             "Don't know/Not sure","Response removed due to invalid value")
disc_ids<-c(40192466,40192489,40192416,40192490,40192380,40192395,40192496,40192519,40192451)
food_ids<-c(40192517,40192426); disab_ids<-c(903573,903574,903575,903576,903577,903578)
ace_ids <-c(1332946,1332947,1332948,1332950,1333208,1332857,1333210,1333212,1332951,1333202,1333203)
trauma_ids<-c(1703988,1704021,1703976,1704055,1704035,1703994,1704056,1704011,1704029,1704020,1703989)
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey` WHERE question_concept_id IN (%s)",
                      cdr, paste(unique(c(disc_ids,food_ids,disab_ids,ace_ids,trauma_ids)),collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]),]; sv$qid <- as.numeric(sv$question_concept_id)
W <- tidyr::pivot_wider(sv[,c("person_id","qid","answer")], names_from=qid, values_from=answer)
getcol<-function(id){c<-W[[as.character(id)]]; if(is.null(c)) rep(NA,nrow(W)) else c}
rc<-function(x,s){x<-as.character(x); x[x%in%na_vals]<-NA; as.numeric(S[[s]][x])}; z<-function(v) as.numeric(scale(v))
idx<-function(ids,s){M<-sapply(ids,function(id) z(rc(getcol(id),s))); z(rowMeans(M,na.rm=TRUE))}
cnt<-function(ids,pos="Yes|Once|More than once",neg="No|Never"){m<-sapply(ids,function(id){a<-as.character(getcol(id)); ifelse(grepl(pos,a),1,ifelse(grepl(neg,a),0,NA))}); z(rowMeans(m,na.rm=TRUE)*ncol(m))}
ADV <- data.frame(person_id=W$person_id, ace=cnt(ace_ids), trauma=cnt(trauma_ids),
                  disability=cnt(disab_ids), discrim=idx(disc_ids,"disc6"), foodinsec=idx(food_ids,"food3"))

## ---------- assemble ----------
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n); d }
d <- collapse(merge(merge(readRDS("cox_dat2.rds"), ADV, by="person_id"), geo, by="person_id"))
need <- c("z_cohesion","ace","trauma","disability","discrim","foodinsec","deprivation_index",
          "age","sex_c","race_c","ethn_c","income_m","educ_m","log_util","event","time_days")
d <- d[complete.cases(d[,need]), ]
d$dep_t <- cut(d$deprivation_index, quantile(d$deprivation_index,0:3/3), include.lowest=TRUE,
               labels=c("Low disadvantage","Middle disadvantage","High disadvantage"))

base_rhs <- "z_cohesion + age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
full_rhs <- paste(base_rhs, "+ ace + discrim + foodinsec + disability + trauma")
one_rhs  <- function(a) paste(base_rhs, "+", a)
row <- function(rhs, data, panel, label){
  m <- coxph(as.formula(sprintf("Surv(time_days,event) ~ %s", rhs)), data)
  b <- summary(m)$coefficients["z_cohesion",]; ci <- exp(confint(m)["z_cohesion",])
  data.frame(panel=panel, model=label, coh_HR=round(exp(b[["coef"]]),3),
             CI=sprintf("%.3f-%.3f", ci[1], ci[2]),
             p=signif(b[[ncol(summary(m)$coefficients)]],3), n=m$n, events=m$nevent) }

## Panel A -- pooled ladder
A <- row(base_rhs, d, "A. Pooled", "Base spec (cohesion + covariates)")
for (a in c("ace","discrim","foodinsec","disability","trauma"))
  A <- rbind(A, row(one_rhs(a), d, "A. Pooled", sprintf("  + %s", a)))
A <- rbind(A, row(full_rhs, d, "A. Pooled", "  + full adversity exposome"))
## Panel B -- by deprivation tertile
B <- do.call(rbind, lapply(levels(d$dep_t), function(dt){ sub<-d[d$dep_t==dt,]
  rbind(row(base_rhs, sub, "B. By area deprivation", paste0(dt, ": base")),
        row(full_rhs, sub, "B. By area deprivation", paste0(dt, ": + full exposome"))) }))
ET <- rbind(A, B)
write.csv(ET, "etable_robustness.csv", row.names=FALSE)

cat("========== eTable. Cohesion (HR per SD) robustness to the adversity exposome ==========\n")
for (pn in unique(ET$panel)){ cat(sprintf("\n-- %s --\n", pn)); s <- ET[ET$panel==pn,]
  for (i in seq_len(nrow(s)))
    cat(sprintf("  %-34s  HR %.3f (%s)  P=%s  [n=%d, events=%d]\n",
                s$model[i], s$coh_HR[i], s$CI[i], format(s$p[i]), s$n[i], s$events[i])) }

## footnote: published primary (full cohort, factor spec) + interaction
ref <- coxph(Surv(time_days,event) ~ z_cohesion + age + sex + race + ethnicity +
               income_f + educ_f + emp_f + log_util, readRDS("cox_dat.rds"))
rb <- summary(ref)$coefficients["z_cohesion",]; rc2 <- exp(confint(ref)["z_cohesion",])
d$dep_z <- z(d$deprivation_index)
mi <- coxph(Surv(time_days,event) ~ z_cohesion*dep_z + age + sex_c + race_c + ethn_c +
              income_m + educ_m + log_util, d)
ib <- summary(mi)$coefficients["z_cohesion:dep_z",]; ic <- exp(confint(mi)["z_cohesion:dep_z",])
cat(sprintf("\nFootnotes:\n  Published primary (full cohort, n=%d, events=%d): HR %.3f (%.3f-%.3f), P=%s\n",
            ref$n, ref$nevent, exp(rb[["coef"]]), rc2[1], rc2[2], signif(rb[[ncol(summary(ref)$coefficients)]],3)))
cat(sprintf("  Cohesion x deprivation interaction (per SD x SD): HR %.3f (%.3f-%.3f), P=%s\n",
            exp(ib[["coef"]]), ic[1], ic[2], signif(ib[[ncol(summary(mi)$coefficients)]],3)))
cat("  Exposure = z_cohesion per SD (higher=more cohesion; HR<1 protective). Re-anchored cohort.\n")
cat("  Covariates: age, sex, race, ethnicity (collapsed), income, education, healthcare utilization.\n")
