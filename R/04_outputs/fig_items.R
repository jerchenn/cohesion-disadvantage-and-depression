# fig_items.R -- item-level decomposition of the cohesion-depression association and its
# resource-dependence. Panel A: each cohesion item's HR per SD, marginally and in a joint 5-item
# model (which items INDEPENDENTLY drive the effect). Panel B: each item's interaction with area
# deprivation and with lower income (HR > 1 = protection fades under disadvantage).
# Aggregate curve/point coordinates only.
#
# INPUT: cox_dat.rds + ds_survey + ds_zip_code_socioeconomic.  RUN AFTER build_cox.R.
# Needs: tidyverse, survival, bigrquery, ggplot2.

library(tidyverse); library(survival); library(bigrquery)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")
cox <- readRDS("cox_dat3.rds"); cox$deprivation_index <- NULL   # geo re-supplies deprivation_index

ids <- c(help=40192463, getalong=40192411, trust=40192499, values=40192417, watchout=40192400)
lab <- c(help="Willing to help", getalong="Get along", trust="Can be trusted",
         values="Share values", watchout="Watch out for each other")

## pull + recode + z-score the 5 items
sv <- run_sql(sprintf("SELECT person_id, question_concept_id, answer FROM `%s.ds_survey`
                       WHERE question_concept_id IN (%s)", cdr, paste(ids, collapse=",")))
sv <- sv[!duplicated(sv[c("person_id","question_concept_id")]), ]
sv$item <- names(ids)[match(as.numeric(sv$question_concept_id), ids)]
wide <- tidyr::pivot_wider(sv[,c("person_id","item","answer")], names_from=item, values_from=answer)
map5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral (neither agree nor disagree)"=3,"Agree"=4,"Strongly agree"=5)
map4 <- c("Strongly disagree"=1,"Disagree"=2,"Agree"=3,"Strongly agree"=4)
na_vals <- c("Don't know/Not sure","PMI: Skip","PMI: Prefer Not To Answer","Skip","Don't know","Prefer not to answer")
rc <- function(x,m){ x<-as.character(x); x[x %in% na_vals]<-NA; as.numeric(m[x]) }
for (nm in c("help","getalong","trust","values")) wide[[nm]] <- rc(wide[[nm]], map5)
wide$watchout <- rc(wide$watchout, map4)
for (nm in names(ids)) wide[[paste0("z_",nm)]] <- as.numeric(scale(wide[[nm]]))
d <- merge(cox, wide[, c("person_id", paste0("z_", names(ids)))], by="person_id")

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
d <- collapse(d)
COV <- "age + sex_c + race_c + ethn_c + income_m + educ_m + log_util"
zitems <- paste0("z_", names(ids))
row3 <- function(m,v){ ci<-confint(m); c(hr=exp(coef(m)[[v]]), lo=exp(ci[v,1]), hi=exp(ci[v,2])) }
mk <- function(nm,panel,series,x) data.frame(item=lab[nm], panel=panel, series=series,
                                             hr=x["hr"], lo=x["lo"], hi=x["hi"], row.names=NULL)

## Panel A: marginal + joint (independent) HRs
mj <- coxph(as.formula(paste("Surv(time_days,event) ~", paste(zitems,collapse=" + "), "+", COV)), d)
rows <- list()
for (v in zitems){
  mm <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "+", COV)), d)
  rows[[length(rows)+1]] <- mk(sub("z_","",v), "A  Association with depression (HR per SD)", "Marginal",            row3(mm,v))
  rows[[length(rows)+1]] <- mk(sub("z_","",v), "A  Association with depression (HR per SD)", "Joint (independent)", row3(mj,v))
}

## Panel B: item x disadvantage interaction (both oriented so HR > 1 = protection fades)
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
gd <- merge(d, geo, by="person_id"); gd <- gd[!is.na(gd$deprivation_index),]; gd$z_dep <- scale(gd$deprivation_index)[,1]
d$z_inc <- scale(d$income_m)[,1]
for (v in zitems){
  md <- coxph(as.formula(paste("Surv(time_days,event) ~", v, "* z_dep +", COV, "+ z_dep")), gd)
  idn <- grep(":", names(coef(md)), value=TRUE)[1]
  rows[[length(rows)+1]] <- mk(sub("z_","",v), "B  Modification by disadvantage (interaction HR)", "x area deprivation", row3(md,idn))
  mi <- coxph(as.formula(paste("Surv(time_days,event) ~", v,
              "* z_inc + age + sex_c + race_c + ethn_c + educ_m + log_util")), d)
  iin <- grep(":", names(coef(mi)), value=TRUE)[1]
  hci <- as.numeric(exp(confint(mi)[iin,]))                       # HR per SD higher income (unnamed)
  inv <- c(hr=1/exp(coef(mi)[[iin]]), lo=1/hci[2], hi=1/hci[1])  # -> per SD LOWER income
  rows[[length(rows)+1]] <- mk(sub("z_","",v), "B  Modification by disadvantage (interaction HR)", "x lower income", inv)
}

dat <- do.call(rbind, rows)
ord <- c("Can be trusted","Willing to help","Get along","Share values","Watch out for each other")
dat$item   <- factor(dat$item,   levels=rev(ord))
dat$series <- factor(dat$series, levels=c("Marginal","Joint (independent)","x area deprivation","x lower income"))
write.csv(dat, "fig_items_data.csv", row.names=FALSE)

p <- ggplot(dat, aes(hr, item, shape=series, colour=series)) +
  geom_vline(xintercept=1, linetype=2, colour="grey55") +
  geom_pointrange(aes(xmin=lo, xmax=hi), position=position_dodge(width=0.55), linewidth=0.5, size=0.4) +
  facet_wrap(~panel, scales="free_x") +
  scale_x_log10() +
  scale_shape_manual(values=c("Marginal"=16,"Joint (independent)"=1,
                              "x area deprivation"=15,"x lower income"=5), name=NULL) +
  scale_colour_manual(values=c("Marginal"="#0072B2","Joint (independent)"="#0072B2",
                               "x area deprivation"="#D55E00","x lower income"="#D55E00"), name=NULL) +
  labs(x="Hazard ratio (95% CI)", y=NULL,
       title="Trust and willingness to help drive the association, and fade most under disadvantage",
       subtitle="A: each item alone and in a joint 5-item model.  B: item x disadvantage interaction (HR>1 = protection fades)") +
  theme_bw(base_size=11) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        strip.text=element_text(hjust=0, face="bold"), plot.title.position="plot",
        plot.title=element_text(face="bold", size=11), plot.subtitle=element_text(size=8.8))
ggsave("figure_items.png", p, width=9.2, height=4.4, dpi=300)
ggsave("figure_items.pdf", p, width=9.2, height=4.4)
cat("Item-level figure written (Panel A drivers; Panel B resource-dependence).\n")
