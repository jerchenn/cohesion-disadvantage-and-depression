# fig_stratified.R -- Cohesion effect WITHIN area-deprivation strata. Adjusted cumulative incidence
# for LOW vs HIGH cohesion, faceted by deprivation tertile. Shows that cohesion separates depression
# risk within every disadvantage level (so it is not merely a disadvantage proxy) while the gap
# narrows as disadvantage rises (effect modification). Curves are covariate-adjusted (Cox model,
# fitted within each stratum, predicted at that stratum's mean/reference profile). Aggregate curve
# coordinates only; truncated where any cell's risk set falls below 20 (All of Us cell-size policy).
#
# INPUT: cox_dat.rds + ds_zip_code_socioeconomic   OUTPUT: fig_stratified_data.csv, figure_stratified.{png,pdf}
# RUN AFTER build_cox.R.   Needs: survival, bigrquery, ggplot2.

library(survival); library(bigrquery); library(ggplot2)
cdr <- Sys.getenv("WORKSPACE_CDR"); proj <- Sys.getenv("GOOGLE_PROJECT")
run_sql <- function(q) bq_table_download(bq_project_query(proj, q), bigint="character")

collapse <- function(d){
  d$sex_c  <- factor(ifelse(d$sex  %in% c("Female","Male"), d$sex, "Other"))
  d$race_c <- factor(ifelse(d$race %in% c("White","Black or African American"), d$race, "Other"))
  d$ethn_c <- factor(ifelse(d$ethnicity=="Hispanic or Latino","Hispanic",
                     ifelse(d$ethnicity=="Not Hispanic or Latino","NotHispanic","Other")))
  d$income_m <- ifelse(is.na(d$income_n), median(d$income_n,na.rm=TRUE), d$income_n)
  d$educ_m   <- ifelse(is.na(d$educ_n),   median(d$educ_n,  na.rm=TRUE), d$educ_n)
  d
}
mode_lvl <- function(f) names(sort(table(f), decreasing=TRUE))[1]

geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
g <- collapse(merge(readRDS("cox_dat.rds"), geo, by="person_id"))
g <- g[!is.na(g$z_cohesion) & !is.na(g$deprivation_index), ]

## cohesion tertiles (keep Low vs High for a clean 2-curve contrast); area-deprivation tertiles
g$coh_t <- cut(g$z_cohesion, quantile(g$z_cohesion, 0:3/3), include.lowest=TRUE,
               labels=c("Low cohesion","Middle","High cohesion"))
g$dep_t <- cut(g$deprivation_index, quantile(g$deprivation_index, 0:3/3), include.lowest=TRUE,
               labels=c("Low disadvantage","Middle disadvantage","High disadvantage"))

## within each deprivation tertile: adjusted Cox, predict Low vs High cohesion at the stratum mean profile
res <- list()
for (dt in levels(g$dep_t)){
  sub <- g[g$dep_t==dt, ]
  fit <- coxph(Surv(time_days,event) ~ coh_t + age + sex_c + race_c + ethn_c +
                 income_m + educ_m + log_util, sub)
  nd <- data.frame(coh_t=factor(c("Low cohesion","High cohesion"), levels=levels(g$coh_t)),
                   age=mean(sub$age), sex_c=mode_lvl(sub$sex_c), race_c=mode_lvl(sub$race_c),
                   ethn_c=mode_lvl(sub$ethn_c), income_m=mean(sub$income_m),
                   educ_m=mean(sub$educ_m), log_util=mean(sub$log_util))
  sf <- survfit(fit, newdata=nd)
  ## truncate where the Low- or High-cohesion risk set drops below 20
  sc <- sub[sub$coh_t %in% c("Low cohesion","High cohesion"), ]; sc$coh_t <- droplevels(sc$coh_t)
  km <- survfit(Surv(time_days,event) ~ coh_t, sc); sid <- rep(seq_along(km$strata), km$strata)
  tmax <- min(sapply(seq_along(km$strata), function(s){
    ix <- sid==s; ok <- km$time[ix][km$n.risk[ix] >= 20]; if(length(ok)) max(ok) else 0 }))
  keep <- sf$time <= tmax
  for (i in 1:2) res[[paste(dt,i)]] <- data.frame(
    years=sf$time[keep]/365.25, dep=dt, cohesion=c("Low cohesion","High cohesion")[i],
    cuminc=1-sf$surv[keep,i], lower=1-sf$upper[keep,i], upper=1-sf$lower[keep,i])
}
out <- do.call(rbind, res)
out$dep      <- factor(out$dep, levels=levels(g$dep_t))
out$cohesion <- factor(out$cohesion, levels=c("Low cohesion","High cohesion"))
write.csv(out, "fig_stratified_data.csv", row.names=FALSE)

pal <- c("Low cohesion"="#D55E00", "High cohesion"="#009E73")
p <- ggplot(out, aes(years, cuminc, colour=cohesion, fill=cohesion, linetype=cohesion)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, group=cohesion), alpha=0.15, colour=NA) +
  geom_step(linewidth=0.9) +
  facet_wrap(~dep, ncol=3) +
  scale_colour_manual(values=pal, name=NULL) +
  scale_fill_manual(values=pal, name=NULL, guide="none") +
  scale_linetype_manual(values=c("Low cohesion"="solid","High cohesion"="dotted"), name=NULL) +
  scale_y_continuous(labels=function(x) paste0(round(100*x), "%")) +
  labs(x="Years since baseline", y="Cumulative incidence of depression",
       title="Cohesion separates depression risk within every disadvantage stratum",
       subtitle="Low vs high neighborhood cohesion, by area-deprivation tertile (adjusted)") +
  theme_bw(base_size=11) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        strip.text=element_text(face="bold"), plot.title.position="plot",
        plot.title=element_text(face="bold", size=11.5), plot.subtitle=element_text(size=9.5))
ggsave("figure_stratified.png", p, width=8.6, height=4.2, dpi=300)
ggsave("figure_stratified.pdf", p, width=8.6, height=4.2)
cat("Stratified cumulative incidence written (low vs high cohesion x area-deprivation tertile).\n")
