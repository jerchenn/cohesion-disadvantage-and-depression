# builtenv_trust.R -- EXPLORATORY. (1) List everything available at ZIP3 level. (2) Link the
# built-environment survey items (safety, disorder, walkability, upkeep) to the two cohesion items
# that drive the depression association (trust, help), and to incident depression itself.
# All built-environment items are oriented so higher = better environment. Aggregate output only.
#
# INPUT: cox_dat.rds + ds_survey + ds_zip_code_socioeconomic.  RUN AFTER build_cox.R.
# Needs: tidyverse, survival, bigrquery.

library(tidyverse); library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat.rds")

## ---------- (1) what is available at ZIP3 level ----------
z1 <- run_sql(sprintf("SELECT * FROM `%s.ds_zip_code_socioeconomic` LIMIT 5", cdr))
cat("=== ds_zip_code_socioeconomic columns (ZIP3-level variables available) ===\n")
print(names(z1))

## ---------- (2) pull built-environment + trust/help items ----------
good <- c(safe=40192384, clean=40192456, upkeep=40192386, sidewalks=40192437,
          shops=40192436, transit=40192440, recreation=40192410, bike=40192431)
bad  <- c(crime=40192493, unsafe_night=40192492, unsafe_day=40192414, graffiti=40192420,
          vandalism=40192412, abandoned=40192469, drug=40192457, alcohol=40192476,
          trouble=40192404, noisy=40192522, hangaround=40192500)
drv  <- c(trust=40192499, help=40192463)
ids  <- c(good, bad, drv)

sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(ids)[match(as.numeric(sv$question_concept_id), ids)]
wide <- tidyr::pivot_wider(sv[,c("person_id","item","answer")], names_from=item, values_from=answer)

## three distinct answer scales are used across these items (confirmed from answer distributions):
map_sw <- c("Strongly disagree"=1,"Somewhat disagree"=2,"Somewhat agree"=3,"Strongly agree"=4)  # walkability + unsafe day/night
map_ad <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)                     # safety/disorder/upkeep (no neutral)
map_5  <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)  # trust, help
sw_items <- c("sidewalks","shops","transit","recreation","bike","unsafe_night","unsafe_day")
d5_items <- c("trust","help")
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Dont Know","PMI: Prefer Not To Answer",
             "Skip","Don't know","Prefer not to answer","Does not apply to my neighborhood")
scale_of <- function(nm) if (nm %in% sw_items) map_sw else if (nm %in% d5_items) map_5 else map_ad
rc <- function(x, m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in names(ids)) wide[[nm]] <- rc(wide[[nm]], scale_of(nm))
## coverage check (should now be uniformly high for every item)
cov <- sapply(names(ids), function(nm) mean(!is.na(wide[[nm]])))
cat("\n=== Item coverage (non-missing after scale-aware recode; expect ~0.9+ for all) ===\n")
print(round(sort(cov), 2))
## z-score each item, then flip the "bad" (disorder/danger) items so higher z = better environment
for (nm in names(ids)) wide[[paste0("z_",nm)]] <- as.numeric(scale(wide[[nm]]))
for (nm in names(bad)) wide[[paste0("z_",nm)]] <- -wide[[paste0("z_",nm)]]

## merge with cohort + ZIP3 deprivation
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
d <- merge(cox, wide[, c("person_id", paste0("z_", names(ids)))], by="person_id")
d <- merge(d, geo, by="person_id", all.x=TRUE)
d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
d$z_dep <- as.numeric(scale(d$deprivation_index))
benv <- names(c(good,bad))

## ---------- built-environment -> trust and help (adjusted std. beta) ----------
beta <- function(y, x) { m <- lm(as.formula(paste(y,"~",x,"+ age + income_m + z_dep")), d)
                         summary(m)$coefficients[x,"Estimate"] }
tab <- data.frame(item=benv,
  to_trust=sapply(paste0("z_",benv), function(x) beta("z_trust", x)),
  to_help =sapply(paste0("z_",benv), function(x) beta("z_help",  x)))
tab <- tab[order(-rowMeans(tab[,c("to_trust","to_help")])), ]
cat("\n=== Built-environment item -> trust / help (std beta, adj. age+income+ZIP3 deprivation) ===\n")
print(cbind(item=tab$item, round(tab[,c("to_trust","to_help")],3)), row.names=FALSE)

## ---------- built-environment -> incident depression (HR per SD, oriented healthier) ----------
cat("\n=== Built-environment item -> incident depression (HR per SD, higher=better env) ===\n")
for (nm in benv){ m <- coxph(as.formula(paste("Surv(time_days,event) ~ z_",nm,
      " + age + income_m + z_dep", sep="")), d)
  cat(sprintf("  %-12s HR %.3f (p=%.2g)\n", nm, exp(coef(m)[[paste0("z_",nm)]]),
      summary(m)$coefficients[paste0("z_",nm),"Pr(>|z|)"])) }

## ---------- ZIP3 deprivation -> trust / help ----------
cat("\n=== ZIP3 area deprivation -> trust / help (std beta per SD deprivation, adj. age+income) ===\n")
for (y in c("z_trust","z_help")){ m<-lm(as.formula(paste(y,"~ z_dep + age + income_m")), d)
  cat(sprintf("  %-8s beta %.3f (p=%.2g)\n", sub("z_","",y),
      coef(m)[["z_dep"]], summary(m)$coefficients["z_dep","Pr(>|t|)"])) }
