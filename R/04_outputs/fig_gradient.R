# fig_gradient.R -- Resource-dependence as a continuous gradient. Plots the cohesion hazard ratio
# across the continua of (A) area deprivation and (B) household income, each from a Cox model with a
# linear cohesion x moderator interaction, with 95% CI bands. Both x-axes are oriented so that
# "increasing disadvantage" runs left -> right; the protective HR rising toward 1.0 on the right is
# the resource-dependence finding. Aggregate curve coordinates only.
#
# INPUT: cox_dat.rds + ds_zip_code_socioeconomic   OUTPUT: fig_gradient_data.csv, figure_gradient.{png,pdf}
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

## cohesion log-HR at moderator value z: b_coh + b_int*z ; SE from the 2x2 covariance
hr_curve <- function(fit, zseq){
  b <- coef(fit); V <- vcov(fit)
  ic <- "z_cohesion"; ii <- grep("z_cohesion:", names(b), value=TRUE)[1]
  est <- b[ic] + b[ii]*zseq
  se  <- sqrt(V[ic,ic] + zseq^2*V[ii,ii] + 2*zseq*V[ic,ii])
  data.frame(hr=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se))
}
pgrid <- 5:95

## ---- (A) area deprivation (higher = more disadvantage) ----
geo <- run_sql(sprintf("SELECT person_id, deprivation_index FROM `%s.ds_zip_code_socioeconomic`", cdr))
geo$deprivation_index <- as.numeric(geo$deprivation_index); geo <- geo[!duplicated(geo$person_id),]
cx <- readRDS("cox_dat3.rds"); cx$deprivation_index <- NULL   # geo re-supplies deprivation_index
g <- collapse(merge(cx, geo, by="person_id"))
g <- g[!is.na(g$z_cohesion) & !is.na(g$deprivation_index), ]; g$z_dep <- scale(g$deprivation_index)[,1]
md <- coxph(Surv(time_days,event) ~ z_cohesion*z_dep + age + sex_c + race_c + ethn_c +
              income_m + educ_m + log_util, g)
zA <- as.numeric(quantile(g$z_dep, pgrid/100))
A  <- cbind(hr_curve(md, zA), disadvantage=pgrid, panel="A  Area deprivation (higher deprivation at right)")

## ---- (B) household income (lower = more disadvantage; x reversed to 100 - income pct) ----
d <- collapse(readRDS("cox_dat3.rds")); d <- d[!is.na(d$z_cohesion), ]; d$z_inc <- scale(d$income_m)[,1]
mi <- coxph(Surv(time_days,event) ~ z_cohesion*z_inc + age + sex_c + race_c + ethn_c +
              educ_m + log_util, d)                    # income is the moderator, not also a covariate
zB <- as.numeric(quantile(d$z_inc, pgrid/100))
B  <- cbind(hr_curve(mi, zB), disadvantage=100-pgrid, panel="B  Household income (lower income at right)")

out <- rbind(A, B)
write.csv(out, "fig_gradient_data.csv", row.names=FALSE)

p <- ggplot(out, aes(disadvantage, hr)) +
  geom_hline(yintercept=1, linetype=2, colour="grey55") +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15) +
  geom_line(linewidth=0.8) +
  facet_wrap(~panel, ncol=2) +
  scale_y_log10(breaks=c(0.80,0.85,0.90,0.95,1.0)) +
  labs(x="Increasing socioeconomic disadvantage (percentile)",
       y="Cohesion hazard ratio per 1-SD (95% CI)",
       title="Protection by cohesion fades with socioeconomic disadvantage") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(hjust=0, face="bold"),
        plot.title.position="plot", plot.title=element_text(face="bold", size=11))
ggsave("figure_gradient.png", p, width=7.6, height=4.4, dpi=300)
ggsave("figure_gradient.pdf", p, width=7.6, height=4.4)
cat("Cohesion-HR gradients written. Deprivation interaction HR",
    sprintf("%.3f", exp(coef(md)[grep('z_cohesion:',names(coef(md)))])),
    "| income interaction HR", sprintf("%.3f", exp(coef(mi)[grep('z_cohesion:',names(coef(mi)))])), "\n")
